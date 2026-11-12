#!/usr/bin/env python3
"""Finish copy-labelled chromosomes with local assembly and supported bridges.

Every accepted change needs unique mappings to the combined copy reference,
local HiFi assembly, unambiguous flanking alignments and multiple raw molecules
spanning the complete replacement. Ambiguous loci remain unresolved.
"""
import argparse
import re
import shutil

# Direct script execution also works outside the repository.
import sys
import tempfile
from collections import Counter, OrderedDict
from contextlib import ExitStack
from pathlib import Path

import pysam
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
# Imports below follow the standalone CLI path bootstrap.
# pylint: disable=wrong-import-position
from assembly.lib.io import assembly_scratch, checksum, execute, load, publish, save
from assembly.lib.postassembly import (
    assembly_bundle,
    assembly_identity,
    inspect_sequences,
    map_reads,
    native_tidk,
    passthrough,
    read_bundle,
    rebase_symlinks,
    retain_files,
    reverse_complement,
    supported_alignments,
    terminal_repeats,
    write_bundle,
    write_sequences,
    write_table,
)

from assembly.lib.sequences import fasta_dict, gfa_to_fasta


def anchor_hits(candidates, anchors, work, threads, identity):
    write_sequences(work / "anchors.fa", anchors)
    with (work / "anchors.paf").open("w", encoding="utf-8") as out:
        execute(["minimap2", "-c", "-x", "asm5", "-t", threads, candidates, work / "anchors.fa"], stdout=out)
    hits = {name: [] for name in anchors}
    for line in (work / "anchors.paf").read_text(encoding="utf-8").splitlines():
        f = line.split("\t")
        if (
            int(f[10])
            and int(f[9]) / int(f[10]) >= identity
            and int(f[2]) <= 5
            and int(f[1]) - int(f[3]) <= 5
        ):
            hits[f[0]].append({"target": f[5], "strand": f[4], "start": int(f[7]), "end": int(f[8])})
    # Multiple equivalent graph paths are ambiguous; never take the longest.
    if any(len(v) != 1 for v in hits.values()):
        return None
    values = {name: v[0] for name, v in hits.items()}
    if len({(v["target"], v["strand"]) for v in values.values()}) != 1:
        return None
    return values


def spanning_support(reference, files, work, threads, _min_support, _min_anchor):
    """Count distinct original molecules covering the entire bridge/extension."""
    names = set()
    refseq = next(iter(fasta_dict(reference).values()))
    for technology, reads in files.items():
        if not reads:
            continue
        bamfile = work / (technology + ".support.bam")
        map_reads(reference, [reads], bamfile, threads, "map-hifi" if technology == "hifi" else "map-ont")
        with pysam.AlignmentFile(bamfile) as bam:
            for r in bam.fetch():
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                if r.reference_start > 20 or r.reference_end < len(refseq) - 20:
                    continue
                aligned = r.query_alignment_length or 0
                if not aligned or not r.has_tag("NM"):
                    continue
                if 1 - r.get_tag("NM") / aligned < (0.97 if technology == "hifi" else 0.80):
                    continue
                # A deletion must not be misread as a read spanning new sequence.
                if any(op in (2, 3) and n > max(50, len(refseq) // 20) for op, n in (r.cigartuples or [])):
                    continue
                names.add(technology + ":" + r.query_name)
    return sorted(names)


def propose(chrom, sequence, start, end, kind, bams, work, options, threads):
    work.mkdir()
    a = options["min_anchor"]
    window = options["local_window"]
    if kind == "gap" and (start < a or len(sequence) - end < a or end - start > options["max_span"]):
        return None, "INSUFFICIENT_FLANK_OR_TOO_LONG"
    lo = max(0, start - window)
    hi = min(len(sequence), end + window)
    files = {}
    counts = {}
    for tech, bam in bams.items():
        path = work / (tech + ".fastq")
        n = 0
        with path.open("w", encoding="utf-8") as out:
            # A real missing segment appears as an insertion relative to the
            # N-gapped draft. Use a relaxed seed filter here; the assembled,
            # complete bridge is independently remapped at the strict threshold.
            for r in supported_alignments(
                bam, chrom, lo, hi, options["min_mapq"], 0.80 if tech == "hifi" else 0.75
            ):
                seq = r.get_forward_sequence()
                quals = r.get_forward_qualities()
                if not seq or quals is None:
                    continue
                out.write(
                    "@" + r.query_name + "\n" + seq + "\n+\n" + "".join(chr(q + 33) for q in quals) + "\n"
                )
                n += 1
        counts[tech] = n
        if n:
            files[tech] = path
    if counts.get("hifi", 0) < options["min_support"]:
        return None, "INSUFFICIENT_UNIQUE_HIFI_READS"
    command = ["hifiasm", "-o", work / "asm", "-t", threads, "-f0", "-l0"]
    if "ont" in files:
        command += ["--ul", files["ont"]]
    command += [files["hifi"]]
    with (work / "hifiasm.log").open("w", encoding="utf-8") as out:
        execute(command, stdout=out, stderr=out)
    gfas = list(work.glob("asm*.p_ctg.gfa"))
    if len(gfas) != 1:
        return None, "LOCAL_ASSEMBLY_NO_UNIQUE_PRIMARY_GRAPH"
    gfa_to_fasta(gfas[0], work / "contigs.fa", False)
    candidates = fasta_dict(work / "contigs.fa")
    if not candidates:
        return None, "LOCAL_ASSEMBLY_EMPTY"
    anchors = OrderedDict()
    if kind == "gap":
        anchors.update(left=sequence[start - a : start], right=sequence[end : end + a])
    elif kind == "left":
        anchors["end"] = sequence[:a]
    else:
        anchors["end"] = sequence[-a:]
    if any(set(s) - set("ACGT") for s in anchors.values()):
        return None, "AMBIGUOUS_FLANK"
    hits = anchor_hits(work / "contigs.fa", anchors, work, threads, options["min_identity"])
    if hits is None:
        return None, "AMBIGUOUS_LOCAL_PATH"
    first = next(iter(hits.values()))
    contig = candidates[first["target"]]
    if first["strand"] == "-":
        contig = reverse_complement(contig)
        for hit in hits.values():
            hit["start"], hit["end"] = len(contig) - hit["end"], len(contig) - hit["start"]
    if kind == "gap":
        left, right = hits["left"], hits["right"]
        if left["end"] >= right["start"]:
            return None, "OVERLAP_JOIN_REQUIRES_REVIEW"
        replacement = contig[left["end"] : right["start"]]
        bridge = anchors["left"] + replacement + anchors["right"]
    elif kind == "left":
        replacement = contig[: hits["end"]["start"]]
        bridge = replacement + anchors["end"]
    else:
        replacement = contig[hits["end"]["end"] :]
        bridge = anchors["end"] + replacement
    if not replacement or len(replacement) > options["max_span"] or set(replacement) - set("ACGT"):
        return None, "EMPTY_AMBIGUOUS_OR_TOO_LONG_REPLACEMENT"
    if kind != "gap" and not terminal_repeats(bridge, options["telomere_motif"])[kind]:
        return None, "NO_TERMINAL_ARRAY_IN_EXTENSION"
    write_sequences(work / "bridge.fa", {"bridge": bridge})
    support = spanning_support(work / "bridge.fa", files, work, threads, options["min_support"], a)
    if len(support) < options["min_support"]:
        return None, "INSUFFICIENT_FULL_SPAN_SUPPORT"
    if kind == "gap":
        # quarTeT acts only on a validated single-copy, single-gap local draft.
        # Its permissive default identity and optional join mode are not used.
        write_sequences(
            work / "draft.fa", {"locus": anchors["left"] + "N" * (end - start) + anchors["right"]}
        )
        execute(
            [
                "quartet",
                "GapFiller",
                "-d",
                "draft.fa",
                "-g",
                "bridge.fa",
                "-f",
                a,
                "-l",
                a - 10,
                "-i",
                options["min_identity"] * 100,
                "-m",
                options["max_span"],
                "-p",
                "filled",
                "-t",
                threads,
                "--noplot",
            ],
            cwd=work,
        )
        filled = fasta_dict(work / "filled.genome.filled.fasta")
        if filled != {"locus": bridge}:
            return None, "QUARTET_RECONSTRUCTION_MISMATCH"
    save(
        work / "accepted.json",
        {
            "supporting_molecules": support,
            "selection_counts": counts,
            "replacement_sha256": __import__("hashlib").sha256(replacement.encode()).hexdigest(),
        },
    )
    return {
        "start": start,
        "end": end,
        "replacement": replacement,
        "kind": kind,
        "support": len(support),
        "evidence": work.name,
    }, "ACCEPTED"


def finish(args):
    meta = load(args.metadata)
    options = meta["settings"]["finishing"]
    assessment = meta["settings"]["assessment"]
    delivery, sequences, groups = assembly_bundle(args.delivery)
    if delivery["sample_id"] != meta["sample_id"]:
        raise ValueError("Source assembly belongs to a different sample")
    copies = delivery.get("output_copies", len(set(groups.values())))
    expected = copies * meta["settings"]["chromosomes_per_haplotype"]
    if len(sequences) != expected:
        raise ValueError("Source chromosome count does not match the configured copy target")
    expected_groups = {group: meta["settings"]["chromosomes_per_haplotype"] for group in set(groups.values())}
    if dict(Counter(groups.values())) != expected_groups:
        raise ValueError("Source chromosome counts are incorrect for individual copies")
    destination = Path(args.out_dir).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with assembly_scratch(destination) as scratch:
        stage = Path(scratch) / "result"
        stage.mkdir()
        (stage / "loci").mkdir()
        write_sequences(stage / "input.fa", sequences)
        shutil.copyfile(Path(delivery["base"]) / "unplaced.fa", stage / "unplaced.fa")
        competitive = OrderedDict(sequences)
        for name, seq in fasta_dict(stage / "unplaced.fa").items():
            label = "unplaced_competitor__" + name
            if label in competitive:
                raise ValueError("Unplaced competitor ID collision")
            competitive[label] = seq
        write_sequences(stage / "competitive.fa", competitive)
        before = inspect_sequences(sequences, assessment, options["telomere_motif"])
        needs = (
            options["gap_closing"]
            and any(r["gaps"] for r in before)
            or options["telomere_completion"]
            and any(not r["telomeres"]["left"] or not r["telomeres"]["right"] for r in before)
        )
        targets = []
        for row in before:
            chrom = row["sequence_id"]
            length = row["length"]
            if options["gap_closing"]:
                targets.extend((chrom, g["start"], g["end"], "gap") for g in row["gaps"])
            if options["telomere_completion"]:
                if not row["telomeres"]["left"]:
                    targets.append((chrom, 0, 0, "left"))
                if not row["telomeres"]["right"]:
                    targets.append((chrom, length, length, "right"))
        attempts = []
        edits = {n: [] for n in sequences}
        if needs:
            with ExitStack() as stack:
                bams = {}
                for tech, reads in [("hifi", args.hifi), ("ont", args.ul)]:
                    if not reads:
                        continue
                    path = stage / (tech + ".bam")
                    map_reads(
                        stage / "competitive.fa",
                        reads,
                        path,
                        args.threads,
                        "map-hifi" if tech == "hifi" else "map-ont",
                    )
                    bams[tech] = stack.enter_context(pysam.AlignmentFile(path))
                for i, (chrom, start, end, kind) in enumerate(targets):
                    change, status = propose(
                        chrom,
                        sequences[chrom],
                        start,
                        end,
                        kind,
                        bams,
                        stage / "loci" / f"locus{i:06d}",
                        options,
                        args.threads,
                    )
                    attempts.append(
                        {"sequence_id": chrom, "start": start, "end": end, "kind": kind, "status": status}
                    )
                    if change:
                        edits[chrom].append(change)
        updated = OrderedDict()
        with (stage / "edits.tsv").open("w", encoding="utf-8") as out:
            out.write(
                "sequence_id\tkind\told_start_0based\told_end_exclusive\t"
                "new_start_0based\tnew_end_exclusive\tsupporting_molecules\tevidence\n"
            )
            for chrom, seq in sequences.items():
                changes = sorted(edits[chrom], key=lambda c: (c["start"], c["end"]))
                offset = 0
                last = 0
                pieces = []
                for change in changes:
                    start, end = change["start"], change["end"]
                    if start < last:
                        raise ValueError("Overlapping finishing edits")
                    pieces.extend([seq[last:start], change["replacement"]])
                    last = end
                    out.write(
                        "\t".join(
                            map(
                                str,
                                [
                                    chrom,
                                    change["kind"],
                                    start,
                                    end,
                                    start + offset,
                                    start + offset + len(change["replacement"]),
                                    change["support"],
                                    "loci/" + change["evidence"],
                                ],
                            )
                        )
                        + "\n"
                    )
                    offset += len(change["replacement"]) - (end - start)
                pieces.append(seq[last:])
                updated[chrom] = "".join(pieces)
        stats = write_bundle(stage, updated, groups)
        native_tidk(stage / "genome.fa", stage / "tidk", options["telomere_motif"])
        after = inspect_sequences(updated, assessment, options["telomere_motif"])
        save(
            stage / "manifest.json",
            {
                **assembly_identity(delivery),
                "sample_id": meta["sample_id"],
                "output_copies": copies,
                "groups": groups,
                "sets": stats,
                "source_delivery_sha256": checksum(args.delivery),
                "source_delivery": str(Path(args.delivery).resolve()),
                "status": "COMPLETED",
                "T2T_status": "NOT_ASSESSED",
                "attempts": attempts,
                "before": before,
                "after": after,
                "accepted_changes": sum(len(v) for v in edits.values()),
                "tools": meta["tools"],
                "code": meta["code"],
                "coordinate_policy": (
                    "edits.tsv maps 06 input coordinates to 07 output;"
                    " original AGP is not valid for modified sequences"
                ),
                "phase_policy": (
                    "Preserve labels; unique competitive chromosome mapping"
                    " is a screening criterion, not independent phasing "
                    "proof"
                ),
            },
        )
        # Temporary local coordinates and read BAMs are evidence, not updated AGP.
        write_table(stage / "attempts.tsv", ["sequence_id", "start", "end", "kind", "status"], attempts)
        if not args.keep_intermediates:
            for locus in (stage / "loci").iterdir():
                if (locus / "accepted.json").is_file():
                    retain_files(locus, {"accepted.json", "bridge.fa", "anchors.paf", "hifiasm.log"})
                else:
                    shutil.rmtree(locus)
            keep = {
                "genome.fa",
                "genome.fa.fai",
                "unplaced.fa",
                "manifest.json",
                "edits.tsv",
                "attempts.tsv",
                "loci",
                "tidk",
            }
            keep.update(group + suffix for group in set(groups.values()) for suffix in (".fa", ".fa.fai"))
            retain_files(stage, keep)
        rebase_symlinks(stage, destination)
        publish(stage, destination, scratch)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata")
    parser.add_argument("--assembly-manifest", "--delivery", dest="delivery", required=True)
    parser.add_argument("--reads-manifest")
    parser.add_argument("--hifi", nargs="*", default=[])
    parser.add_argument("--ul", nargs="*", default=[])
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--gap-closing", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--telomere-completion", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--telomere-motif")
    for field in ("min-anchor", "min-mapq", "min-support", "local-window", "max-span"):
        parser.add_argument("--" + field, type=int)
    parser.add_argument("--min-identity", type=float)
    parser.add_argument("--keep-intermediates", action="store_true")
    args = parser.parse_args()
    defaults = yaml.safe_load(
        Path(__file__).resolve().parents[1].joinpath("defaults.yaml").read_text(encoding="utf-8")
    )
    source, _, _ = assembly_bundle(args.delivery)
    meta = (
        load(args.metadata)
        if args.metadata
        else {"sample_id": source["sample_id"], "settings": defaults, "tools": {}, "code": {}}
    )
    options = meta["settings"]["finishing"]
    for field in (
        "gap_closing",
        "telomere_completion",
        "telomere_motif",
        "min_anchor",
        "min_identity",
        "min_mapq",
        "min_support",
        "local_window",
        "max_span",
    ):
        if getattr(args, field) is not None:
            options[field] = getattr(args, field)
    if args.reads_manifest:
        reads = read_bundle(args.reads_manifest, source["sample_id"])
        args.hifi = [
            p
            for record in reads["readsets"].values()
            if record["technology"] == "hifi"
            for p in record["paths"]
        ]
        args.ul = [
            p
            for record in reads["readsets"].values()
            if record["technology"] == "ont"
            for p in record["paths"]
        ]
    meta["settings"]["ploidy"] = source.get("ploidy")
    meta["settings"]["chromosomes_per_haplotype"] = source["expected_chromosomes_per_haplotype"]
    if not options["gap_closing"] and not options["telomere_completion"]:
        destination = Path(args.out_dir).absolute()
        destination.parent.mkdir(parents=True, exist_ok=True)
        with assembly_scratch(destination) as temporary:
            stage = Path(temporary) / "result"
            stage.mkdir()
            passthrough(args.delivery, stage, "DISABLED")
            write_table(
                stage / "edits.tsv",
                [
                    "sequence_id",
                    "kind",
                    "old_start_0based",
                    "old_end_exclusive",
                    "new_start_0based",
                    "new_end_exclusive",
                    "supporting_molecules",
                    "evidence",
                ],
                [],
            )
            publish(stage, destination, temporary)
        return
    if not args.hifi or not options["telomere_motif"]:
        parser.error("Active finishing requires HiFi reads and --telomere-motif")
    if not 0.95 <= options["min_identity"] <= 1 or options["min_support"] < 2 or options["min_anchor"] < 500:
        parser.error("Invalid anchor, identity or support threshold")
    if (
        args.threads < 4
        or not 20 <= options["min_mapq"] <= 60
        or options["local_window"] < max(1000, options["min_anchor"])
        or options["max_span"] < 1000
        or not re.fullmatch("[ACGT]{2,20}", options["telomere_motif"])
    ):
        parser.error("Invalid resources, mapping quality, window, span or telomere motif")
    destination = Path(args.out_dir).absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".finish-config-", dir=destination.parent) as temporary:
        args.metadata = Path(temporary) / "metadata.json"
        save(args.metadata, meta)
        finish(args)


if __name__ == "__main__":
    main()
