"""HapHiC scaffolding and validated chromosome exports for all strategies."""

import csv
import os
import re
import shutil
import tempfile
from collections import Counter
from contextlib import ExitStack
from pathlib import Path

import pysam

from assembly.lib.hic import bam_pairs
from assembly.lib.io import assembly_scratch
from assembly.lib.io import checksum, execute, load, publish, save
from assembly.lib.sequences import (
    complete_delivery,
    fasta_dict,
    fasta_statistics,
    parse_agp,
    reconstruct,
    write_fasta,
)


def separate_scaffolding(args, metadata, stage):
    options = metadata["settings"]["scaffolding"]
    inputs = stage / "inputs"
    inputs.mkdir(parents=True)
    fastas = [inputs / f"hap{i}.fa" for i in (1, 2)]
    with fastas[0].open("w", encoding="utf-8") as one, fastas[1].open("w", encoding="utf-8") as two:
        with pysam.FastxFile(str(args.fasta)) as source:
            for record in source:
                hap = record.name.split("_", 1)[0]
                if hap not in ("hap1", "hap2"):
                    raise ValueError("Per-haplotype scaffolding requires hap1_/hap2_ contig IDs")
                write_fasta(one if hap == "hap1" else two, record.name, record.sequence)
    bams = [inputs / f"hap{i}.bam" for i in (1, 2)]
    with pysam.AlignmentFile(str(args.bam), "rb") as source:
        headers = []
        for i in (1, 2):
            header = source.header.to_dict()
            header["SQ"] = [sq for sq in header["SQ"] if sq["SN"].startswith(f"hap{i}_")]
            headers.append(pysam.AlignmentHeader.from_dict(header))
        with pysam.AlignmentFile(str(bams[0]), "wb", header=headers[0]) as one, pysam.AlignmentFile(
            str(bams[1]), "wb", header=headers[1]
        ) as two:
            for read in source:
                hap = read.reference_name.split("_", 1)[0]
                if hap not in ("hap1", "hap2") or not read.next_reference_name.startswith(hap + "_"):
                    raise ValueError("Per-haplotype linkage contains an inter-haplotype or unmapped mate")
                index = 0 if hap == "hap1" else 1
                (one if index == 0 else two).write(
                    pysam.AlignedSegment.fromstring(read.to_string(), headers[index])
                )
    commands = []
    for i in (1, 2):
        hap = f"hap{i}"
        bam_pairs(bams[i - 1])
        gfa = [p for p in args.gfas if Path(p).name == hap + ".gfa"]
        if len(gfa) != 1:
            raise ValueError(f"Expected one matching {hap}.gfa")
        command = [
            "haphic",
            "pipeline",
            fastas[i - 1].resolve(),
            bams[i - 1].resolve(),
            metadata["settings"]["chromosomes_per_haplotype"],
            "--gfa",
            Path(gfa[0]).resolve(),
            "--RE",
            options["restriction_sites"],
            "--correct_nrounds",
            options["correct_nrounds"],
            "--Ns",
            options["gap_size"],
            "--threads",
            args.threads,
            "--processes",
            min(options["processes"], args.threads),
            "--outdir",
            stage / hap,
        ]
        execute(command)
        commands.append([str(x) for x in command])
    build = stage / "04.build"
    build.mkdir()
    tours = stage / "03.sort/final_tours"
    tours.mkdir(parents=True)
    for filename in ("scaffolds.fa", "scaffolds.raw.agp", "scaffolds.agp"):
        with (build / filename).open("w", encoding="utf-8") as output:
            for i in (1, 2):
                hap = f"hap{i}"
                path = stage / hap / "04.build" / filename
                if filename.endswith(".fa"):
                    for name, seq in fasta_dict(path).items():
                        write_fasta(output, hap + "_" + name, seq)
                else:
                    for line in path.read_text(encoding="utf-8").splitlines():
                        if not line or line.startswith("#"):
                            continue
                        fields = line.split("\t")
                        fields[0] = hap + "_" + fields[0]
                        output.write("\t".join(fields) + "\n")
    for i in (1, 2):
        hap = f"hap{i}"
        for tour in (stage / hap / "03.sort/final_tours").glob("*.tour"):
            shutil.copyfile(tour, tours / (hap + "_" + tour.name))
    (build / "juicebox.sh").write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        "# Each haplotype has its own native Juicebox helper and review coordinates.\n"
        'HERE="$(cd -- "$(dirname -- "$0")" && pwd)"\n'
        'for hap in hap1 hap2; do (cd "$HERE/../$hap/04.build" && bash juicebox.sh); done\n',
        encoding="utf-8",
    )
    save(
        stage / "separate_haplotypes.json",
        {
            "status": "AUTOMATIC_CANDIDATE",
            "commands": commands,
            "phase_validation": "NOT_ASSESSED",
            "shared_physical_molecules": True,
            "limitations": [
                "Structural differences between haplotypes can bias shared Hi-C mappings",
                "No inter-haplotype contacts were constructed; this does not prove phase correctness",
            ],
        },
    )
    return commands


def scaffold(args):
    metadata = load(args.metadata)
    options = metadata["settings"]["scaffolding"]
    destination = Path(args.out_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with assembly_scratch(destination, prefix=".haphic-") as temp:
        stage = Path(temp) / "haphic"
        mode = load(args.linkage).get("hic_reference", "combined")
        command = [
            "haphic",
            "pipeline",
            Path(args.fasta).resolve(),
            Path(args.bam).resolve(),
            metadata["settings"].get("output_copies", 2) * metadata["settings"]["chromosomes_per_haplotype"],
            "--RE",
            options["restriction_sites"],
            "--correct_nrounds",
            options["correct_nrounds"],
            "--Ns",
            options["gap_size"],
            "--threads",
            args.threads,
            "--processes",
            min(options["processes"], args.threads),
            "--outdir",
            stage,
        ]
        if args.gfas:
            command += ["--gfa", ",".join(str(Path(p).resolve()) for p in args.gfas), "--phasing_weight", "1"]
        if mode == "per_haplotype":
            commands = separate_scaffolding(args, metadata, stage)
        else:
            execute(command)
            commands = [[str(x) for x in command]]
        for name in ("scaffolds.fa", "scaffolds.raw.agp", "scaffolds.agp", "juicebox.sh"):
            if not (stage / "04.build" / name).is_file():
                raise ValueError(f"Missing HapHiC build artifact: {name}")
        save(
            stage / "manifest.json",
            {
                "status": "AUTOMATIC_CANDIDATE",
                "command": commands[0] if len(commands) == 1 else None,
                "settings": options,
                "commands_executed": commands,
                "input_fasta_sha256": checksum(args.fasta),
                "input_bam_sha256": load(args.linkage)["bam_sha256"],
                "hic_reference": load(args.linkage).get("hic_reference", "combined"),
                "phase_validation": "NOT_ASSESSED",
                "automatic_correction_rounds_requested": options["correct_nrounds"],
                "human_review": "NOT_REVIEWED",
            },
        )
        # HapHiC embeds temporary paths in Juicebox helper. Replace with published paths.
        for helper in stage.rglob("*.sh"):
            if not helper.is_symlink():
                helper.write_text(
                    helper.read_text(encoding="utf-8").replace(str(stage), str(destination.resolve())),
                    encoding="utf-8",
                )
        for path in stage.rglob("*"):
            if path.is_symlink():
                target = os.readlink(path)
                if target.startswith(str(stage)):
                    path.unlink()
                    path.symlink_to(target.replace(str(stage), str(destination.resolve()), 1))
        publish(stage, destination, temp)


def deliver(args):
    metadata, assembly = load(args.metadata), load(args.assembly)
    options = metadata["settings"]["scaffolding"]
    sequences = fasta_dict(args.fasta)
    if checksum(args.fasta) != assembly["combined_sha256"]:
        raise ValueError("Joint assembly FASTA differs from its source manifest")
    haps = {}
    with Path(args.mapping).open(encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            haps[row["contig_id"]] = row["haplotype"]
    if set(haps) != set(sequences):
        raise ValueError("Haplotype table does not match source FASTA")
    raw_agp = Path(args.haphic) / "04.build/scaffolds.raw.agp"
    automatic = parse_agp(raw_agp, sequences)
    native = fasta_dict(Path(args.haphic) / "04.build/scaffolds.fa")
    if set(native) != set(automatic):
        raise ValueError("Native HapHiC FASTA and AGP objects differ")
    for obj, rows in automatic.items():
        if reconstruct(rows, sequences) != native[obj]:
            raise ValueError(f"Native HapHiC FASTA disagrees with AGP: {obj}")
    # Only objects derived from explicit HapHiC tours are anchored chromosome candidates.
    tours = list((Path(args.haphic) / "03.sort/final_tours").glob("*.tour"))
    anchored = {p.stem.rsplit("_", 1)[0] for p in tours}
    if not anchored or not anchored <= set(automatic):
        raise ValueError("Cannot trace HapHiC chromosome objects to sorted tours")
    reviewed = False
    chosen_agp, objects = raw_agp, automatic
    if options["review_agp"]:
        review = load(options["review_manifest"])
        required = {
            "status",
            "reviewer",
            "reviewed_at",
            "input_fasta_sha256",
            "automatic_agp_sha256",
            "review_agp_sha256",
            "chromosome_objects",
        }
        if (
            not required <= set(review)
            or review["status"] != "APPROVED"
            or not review["reviewer"]
            or not review["reviewed_at"]
        ):
            raise ValueError("Manual review manifest is incomplete or not APPROVED")
        if (
            review["input_fasta_sha256"] != checksum(args.fasta)
            or review["automatic_agp_sha256"] != checksum(raw_agp)
            or review["review_agp_sha256"] != checksum(options["review_agp"])
        ):
            raise ValueError("Manual review is stale or refers to a different assembly/AGP")
        chosen_agp = Path(options["review_agp"])
        objects = parse_agp(chosen_agp, sequences)
        anchored = set(review["chromosome_objects"])
        if not anchored or not anchored <= set(objects):
            raise ValueError("Reviewed chromosome object list is invalid")
        reviewed = True
    strategy = metadata["settings"].get("assembly", {}).get("strategy", "diploid")
    labels = ["reference"] if strategy == "single" else ["hap1", "hap2"]
    assigned, unresolved = {label: [] for label in labels}, []
    for obj, rows in objects.items():
        hapset = {haps[r[5]] for r in rows if r[4] == "W"}
        if obj in anchored and len(hapset) == 1:
            assigned[next(iter(hapset))].append(obj)
        else:
            unresolved.append(obj)
    expected = metadata["settings"]["chromosomes_per_haplotype"]
    counts = {hap: len(v) for hap, v in assigned.items()}
    quality = (
        ("CHROMOSOME_REVIEWED" if reviewed else "CHROMOSOME_CANDIDATE")
        if all(count == expected for count in counts.values())
        else "INCOMPLETE_SCAFFOLDING"
    )
    destination = Path(args.out_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".delivery-", dir=destination.parent) as temp:
        stage = Path(temp) / "delivery"
        stage.mkdir()
        shutil.copyfile(args.fasta, stage / "components.fa")
        names = {}
        for hap, group in assigned.items():
            # Local names reflect length order, not syntenic correspondence to reference chromosomes.
            group.sort(key=lambda o: int(objects[o][-1][2]), reverse=True)
            names.update({o: f"{hap}_Chr{i:02d}" for i, o in enumerate(group, 1)})
        names.update({o: f"unplaced_{i:06d}" for i, o in enumerate(unresolved, 1)})
        with ExitStack() as stack:
            handles = {
                k: stack.enter_context((stage / f"{k}.fa").open("w", encoding="utf-8"))
                for k in ["chromosomes", *[h + ".chromosomes" for h in labels], "unplaced"]
            }
            agp = stack.enter_context((stage / "assembly.agp").open("w", encoding="utf-8"))
            table = stack.enter_context((stage / "chromosome_names.tsv").open("w", encoding="utf-8"))
            table.write("sequence_id\tsource_object\thaplotype\tstatus\tlength\n")
            for obj, name in names.items():
                rows = objects[obj]
                sequence = reconstruct(rows, sequences)
                hap = next((h for h in assigned if obj in assigned[h]), None)
                if hap:
                    write_fasta(handles["chromosomes"], name, sequence)
                    write_fasta(handles[hap + ".chromosomes"], name, sequence)
                else:
                    write_fasta(handles["unplaced"], name, sequence)
                for row in rows:
                    line = "\t".join([name, *row[1:]]) + "\n"
                    agp.write(line)
                table.write(
                    "\t".join(
                        [
                            name,
                            obj,
                            hap or "unresolved",
                            "chromosome_candidate" if hap else "unplaced",
                            str(len(sequence)),
                        ]
                    )
                    + "\n"
                )
        save(
            stage / "review.request.json",
            {
                "status": "APPROVED" if reviewed else "PENDING",
                "reviewer": None,
                "reviewed_at": None,
                "input_fasta_sha256": checksum(args.fasta),
                "automatic_agp_sha256": checksum(raw_agp),
                "review_agp_sha256": None,
                "chromosome_objects": sorted(anchored),
                "instructions": (
                    "Review contact maps and original-coordinate AGP; "
                    "supply external review_agp and a completed manifest. "
                    "Do not approve by changing status alone."
                ),
            },
        )
        if reviewed:
            shutil.copyfile(options["review_manifest"], stage / "review.applied.json")
        result = {
            "sample_id": metadata["sample_id"],
            "assembly_id": assembly["assembly_id"],
            "ploidy": metadata["settings"].get("ploidy", 2),
            "strategy": strategy,
            "output_copies": len(labels),
            "phase_level": assembly.get("phase_level", "HIC_PHASED_UNVALIDATED"),
            "quality_level": quality,
            "human_review": "APPROVED" if reviewed else "NOT_REVIEWED",
            "T2T_status": "NOT_ASSESSED",
            "chromosome_counts": counts,
            "hic_reference": load(Path(args.haphic) / "manifest.json").get("hic_reference", "combined"),
            "phase_validation": "NOT_ASSESSED",
            "expected_chromosomes_per_haplotype": expected,
            "unplaced_objects": len(unresolved),
            "mixed_haplotype_objects": [
                o for o in anchored if len({haps[r[5]] for r in objects[o] if r[4] == "W"}) > 1
            ],
            "input_fasta_sha256": checksum(args.fasta),
            "applied_agp_sha256": checksum(chosen_agp),
            "source_coverage_validation": "PASS_EXACTLY_ONCE_ALL_BASES",
            "native_fasta_agp_concordance": "PASS",
            "sets": {
                name: fasta_statistics(stage / f"{name}.fa")
                for name in ["chromosomes", *[h + ".chromosomes" for h in labels], "unplaced"]
            },
            "limitations": [
                "Chromosome numbers are length ranks, not reference chromosome identities",
                "Hi-C gaps contain Ns; this is not a T2T or gap-free assembly",
                "Hi-C phasing is not independent phase validation; hap labels are not parents",
                "BUSCO, k-mer QV, contamination and telomere completeness are not assessed",
            ],
        }
        if result["hic_reference"] == "per_haplotype":
            result["limitations"].append(
                "Hi-C molecules were reused within each haplotype reference; "
                "absence of inter-haplotype contacts is by construction, "
                "not evidence of correct phasing"
            )
        complete_delivery(stage, result, None if strategy == "single" else "hifiasm", ["HapHiC"])
        if strategy == "single":
            result["phase_validation"] = "NOT_APPLICABLE"
        save(stage / "manifest.json", result)
        with (stage / "SHA256SUMS").open("w", encoding="utf-8") as out:
            for p in sorted(stage.iterdir()):
                if p.is_file() and p.name != "SHA256SUMS":
                    out.write(f"{checksum(p)}  {p.name}\n")
        publish(stage, destination, temp)


def export_polyploid(args):
    meta = load(args.metadata)
    assembly = load(args.assembly)
    base = Path(args.gphase) / "gphase_output/scaffold_hap"
    raw = fasta_dict(args.raw_unitigs)
    if checksum(args.raw_unitigs) != assembly["fasta_sha256"]:
        raise ValueError("Raw unitigs changed")
    # Validate both the inspectable unitig representation and graph-merged contigs.
    components = {}
    objects = {}
    native = {}
    for representation in ("unitig", "contig"):
        components[representation] = fasta_dict(base / f"gphase_final.{representation}.fasta")
        objects[representation] = parse_agp(
            base / f"gphase_final.{representation}.scaffold.agp",
            components[representation],
            require_complete=False,
        )
        native[representation] = fasta_dict(base / f"gphase_final.{representation}.scaffold.fasta")
        if set(objects[representation]) != set(native[representation]):
            raise ValueError("GPhase AGP/FASTA object names disagree")
        for name, rows in objects[representation].items():
            if reconstruct(rows, components[representation]) != native[representation][name]:
                raise ValueError(f"GPhase native FASTA differs from AGP: {name}")
    used_raw = set()
    for name, seq in components["unitig"].items():
        original = re.sub(r"_dup\d+$", "", name)
        if original not in raw or seq != raw[original]:
            raise ValueError(f"GPhase unitig copy differs from its raw sequence: {name}")
        used_raw.add(original)
    chroms = meta["settings"]["chromosomes_per_haplotype"]
    ploidy = meta["settings"]["ploidy"]
    expected = {(c, h) for c in range(1, chroms + 1) for h in range(1, ploidy + 1)}
    groups = Counter()
    names = {}
    unclassified = []
    for name in native["contig"]:
        match = re.match(r"^chr(\d+)g(\d+)(?:_|$)", name)
        if not match:
            unclassified.append(name)
            continue
        pair = tuple(map(int, match.groups()))
        if pair not in expected:
            raise ValueError(f"Unexpected chromosome/haplotype label: {name}")
        groups[pair] += 1
        names[name] = f"Chr{pair[0]:02d}_hap{pair[1]}_part{groups[pair]}"
    complete = set(groups) == expected and all(n == 1 for n in groups.values()) and not unclassified
    # No selection of a largest scaffold to force the expected chromosome count.
    dest = Path(args.out_dir)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with assembly_scratch(dest) as temp:
        stage = Path(temp) / "delivery"
        stage.mkdir()
        labels = ["chromosomes", *[f"hap{h}.chromosomes" for h in range(1, ploidy + 1)], "unplaced"]
        with ExitStack() as stack:
            files = {
                label: stack.enter_context((stage / f"{label}.fa").open("w", encoding="utf-8"))
                for label in labels
            }
            for name, seq in native["contig"].items():
                if name not in names:
                    write_fasta(files["unplaced"], name, seq)
                    continue
                hap = int(re.match(r"^chr\d+g(\d+)", name).group(1))
                output_name = (
                    names[name].removesuffix("_part1")
                    if groups[(int(re.match(r"^chr(\d+)", name).group(1)), hap)] == 1
                    else names[name]
                )
                names[name] = output_name
                write_fasta(files["chromosomes"], output_name, seq)
                write_fasta(files[f"hap{hap}.chromosomes"], output_name, seq)
            for name, seq in raw.items():
                if name not in used_raw:
                    write_fasta(files["unplaced"], "unassigned_unitig_" + name, seq)
        shutil.copyfile(base / "gphase_final.contig.fasta", stage / "components.fa")
        shutil.copyfile(base / "gphase_final.unitig.fasta", stage / "unitig_components.fa")
        shutil.copyfile(base / "gphase_final.unitig.scaffold.agp", stage / "unitig_assembly.agp")
        shutil.copyfile(base / "gphase_final.ctg2utg.txt", stage / "contig_to_unitig.txt")
        with (stage / "assembly.agp").open("w", encoding="utf-8") as out:
            for name, rows in objects["contig"].items():
                for row in rows:
                    out.write("\t".join([names.get(name, name), *row[1:]]) + "\n")
        with (stage / "chromosome_names.tsv").open("w", encoding="utf-8") as out:
            out.write("sequence_id\tsource_object\thaplotype\tstatus\tlength\n")
            for name, seq in native["contig"].items():
                hap = re.match(r"^chr\d+g(\d+)", name) if name in names else None
                out.write(
                    f'{names.get(name,name)}\t{name}\t{("hap"+hap.group(1)) if hap else "unresolved"}\t'
                    f'{"chromosome_candidate" if hap else "unplaced"}\t{len(seq)}\n'
                )
            for name, seq in raw.items():
                if name not in used_raw:
                    out.write(f"unassigned_unitig_{name}\t{name}\tunresolved\tunplaced\t{len(seq)}\n")
        execute(["samtools", "faidx", stage / "unitig_components.fa"])
        quality = "POLYPLOID_CHROMOSOME_CANDIDATE" if complete else "INCOMPLETE_SCAFFOLDING"
        result = {
            "sample_id": meta["sample_id"],
            "assembly_id": assembly["assembly_id"],
            "ploidy": ploidy,
            "strategy": "polyploid",
            "output_copies": ploidy,
            "expected_chromosomes_per_haplotype": chroms,
            "quality_level": quality,
            "human_review": "NOT_REVIEWED",
            "phase_level": "GPHASE_CHROMOSOME_PHASED_UNVALIDATED",
            "T2T_status": "NOT_ASSESSED",
            "chromosome_group_counts": {f"Chr{c:02d}_hap{h}": groups[(c, h)] for c, h in sorted(expected)},
            "chromosome_counts": {
                f"hap{h}": sum(n for (c, k), n in groups.items() if k == h) for h in range(1, ploidy + 1)
            },
            "unassigned_raw_unitigs": len(set(raw) - used_raw),
            "unclassified_scaffolds": unclassified,
            "native_fasta_agp_concordance": "PASS_BOTH_UNITIG_AND_CONTIG_REPRESENTATIONS",
            "unitig_copy_sequence_validation": "PASS_IDENTICAL_TO_RAW_INPUT",
            "source_coverage_validation": (
                "GRAPH_COORDINATES_BOUNDS_NONOVERLAP; " "RAW_GRAPH_OVERLAP_TRIMMING_ALLOWED"
            ),
            "sets": {label: fasta_statistics(stage / f"{label}.fa") for label in labels},
            "limitations": [
                "Copy labels are chromosome-local; no parental or C/D/K/L subgenome identity assigned",
                "Copy number and phase accuracy require independent k-mer and sequence validation",
                (
                    "Graph overlap trimming and duplicated collapsed unitigs follow GPhase;"
                    " no claim of exactly-once raw graph bases"
                ),
                "Unassigned raw unitigs are retained, including possible graph redundancy",
                "Scaffolds can contain Ns; telomere and T2T completeness are not certified",
                "BUSCO, k-mer QV, contamination and manual Hi-C curation remain required",
            ],
        }
        complete_delivery(stage, result, "GPhase", ["YaHS", "HapHiC.sort"])
        save(stage / "manifest.json", result)
        save(
            stage / "review.request.json",
            {
                "status": "NOT_REVIEWED",
                "assembly_id": assembly["assembly_id"],
                "components_sha256": checksum(stage / "components.fa"),
                "automatic_agp_sha256": checksum(stage / "assembly.agp"),
                "instructions": (
                    "Review remapped Hi-C contacts and graph/AGP; reviewed "
                    "polyploid AGP import is not implemented"
                ),
            },
        )
        with (stage / "SHA256SUMS").open("w", encoding="utf-8") as out:
            for path in sorted(stage.iterdir()):
                if path.is_file() and path.name != "SHA256SUMS":
                    out.write(f"{checksum(path)}  {path.name}\n")
        publish(stage, dest, temp)
