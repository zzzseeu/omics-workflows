#!/usr/bin/env python3
"""Apply the agreed technology-specific read QC and publish one clean-read manifest."""
import argparse
import shutil

# Direct script execution also works outside the repository.
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
# Imports below follow the standalone CLI path bootstrap.
# pylint: disable=wrong-import-position
from assembly.lib.io import assembly_scratch, checksum, execute, load, pipeline, publish, save
from assembly.lib.postassembly import artifact, fastq_statistics, read_bundle, write_table


def filter_paired(record, _work, output, reports, options, threads):
    files = [output / f'{record["readset_id"]}.R{mate}.fastq.gz' for mate in (1, 2)]
    execute(
        [
            "fastp",
            "-i",
            record["paths"][0],
            "-I",
            record["paths"][1],
            "-o",
            files[0],
            "-O",
            files[1],
            "-w",
            min(threads, 16),
            "--length_required",
            options["fastp_min_length"],
            "--qualified_quality_phred",
            options["fastp_quality"],
            "--json",
            reports / "fastp.json",
            "--html",
            reports / "fastp.html",
        ]
    )
    return files


def filter_long(record, work, output, reports, options, threads):
    source = Path(record["paths"][0])
    target = output / f'{record["readset_id"]}.fastq.gz'
    execute(
        [
            "NanoPlot",
            "--fastq",
            source,
            "--outdir",
            reports,
            "--threads",
            threads,
            "--no_static",
            "--plots",
            "dot",
            "--tsv_stats",
        ]
    )
    if record["technology"] == "hifi":
        native = work / "adapter"
        native.mkdir()
        # The upstream shell script decompresses its inputs. Never give it a
        # link to an immutable producer or a user-owned raw input.
        shutil.copyfile(source, native / "reads.fastq.gz")
        execute(
            [
                "hifiadapterfilt.sh",
                "-p",
                "reads",
                "-t",
                threads,
                "-l",
                options["hifi_adapter_length"],
                "-m",
                options["hifi_adapter_identity"],
                "-o",
                "filtered",
            ],
            cwd=native,
        )
        candidates = list(native.rglob("*.filt.fastq.gz"))
        if len(candidates) != 1:
            raise ValueError("HiFiAdapterFilt did not produce exactly one filtered FASTQ")
        shutil.move(candidates[0], target)
        for stat in native.rglob("*.stats"):
            shutil.copyfile(stat, reports / "adapter.stats")
    else:
        pipeline(
            [
                ["gzip", "-dc", source],
                ["NanoFilt", "-q", options["ont_min_quality"], "-l", options["ont_min_length"]],
                ["pigz", "-p", max(1, threads - 2), "-c"],
            ],
            target,
        )
    # HTML reports embed their figures; retain HTML assets and summary only.
    for path in reports.iterdir():
        if path.is_file() and path.suffix not in (".html", ".txt", ".tsv", ".stats"):
            path.unlink()
    return [target]


def quality_sample(args):
    metadata = load(args.metadata)
    prepared = read_bundle(args.inputs_manifest, metadata["sample_id"])
    options = metadata["settings"]["read_qc"]
    destination = Path(args.out_dir).absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with assembly_scratch(destination) as temporary:
        stage = Path(temporary) / "result"
        output = stage / "reads"
        output.mkdir(parents=True)
        records = {}
        for rid, record in prepared["readsets"].items():
            work = Path(temporary) / rid
            work.mkdir()
            reports = stage / "reports" / rid
            reports.mkdir(parents=True)
            operation = filter_paired if record["technology"] in ("hic", "wgs") else filter_long
            files = operation(record, work, output, reports, options, args.threads)
            statistics = fastq_statistics(files)
            records[rid] = {key: record[key] for key in ("readset_id", "library_id", "technology")}
            records[rid].update(
                files=[artifact(path, stage) for path in files],
                statistics=statistics,
                input_statistics=record["statistics"],
                filtering=(
                    "fastp"
                    if len(files) == 2
                    else "HiFiAdapterFilt" if record["technology"] == "hifi" else "NanoFilt"
                ),
            )
        write_table(
            stage / "reads.tsv",
            [
                "readset_id",
                "technology",
                "filtering",
                "num_reads",
                "num_bases",
                "mean_length",
                "read_n50",
                "num_pairs",
            ],
            [{**record, **record["statistics"]} for record in records.values()],
        )
        save(
            stage / "manifest.json",
            {
                "sample_id": metadata["sample_id"],
                "status": "PASS",
                "readsets": records,
                "parameters": options,
                "input_manifest_sha256": checksum(args.inputs_manifest),
                "tools": metadata.get("tools", {}),
            },
        )
        publish(stage, destination, temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs-manifest", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--threads", type=int, default=8)
    quality_sample(parser.parse_args())


if __name__ == "__main__":
    main()
