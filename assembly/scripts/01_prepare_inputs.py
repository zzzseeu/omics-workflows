#!/usr/bin/env python3
"""Check and prepare all read technologies through one streaming entry point."""
import argparse
import gzip
import shutil
import sqlite3

# Direct script execution also works outside the repository.
import sys
from pathlib import Path
from types import SimpleNamespace

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
# Imports below follow the standalone CLI path bootstrap.
# pylint: disable=wrong-import-position
from assembly.lib.io import assembly_scratch, checksum, load, publish, save, source_check
from assembly.lib.postassembly import artifact, fastq_statistics, write_table
from assembly.lib.reads import prepare_long_reads as prepare
from assembly.lib.reads import reuse


def paired_files(spec, folder):
    result = []
    sources = []
    for number in (1, 2):
        source = Path(spec[f"file{number}"])
        observed = checksum(source, "md5")
        expected = spec["expected_md5" if number == 1 else "expected_md5_2"]
        if expected and expected != observed:
            raise ValueError(f"MD5 mismatch: {source}")
        target = folder / f'{spec["readset_id"]}.R{number}.fastq.gz'
        with source.open("rb") as handle:
            compressed = handle.read(2) == b"\x1f\x8b"
        if compressed:
            target.symlink_to(source.resolve())
        else:
            with source.open("rb") as handle, gzip.open(target, "wb", compresslevel=3) as out:
                shutil.copyfileobj(handle, out)
        result.append(target)
        sources.append(
            {
                "file": str(source),
                "md5": observed,
                "sha256": checksum(source),
                "provided_md5_status": "PASS" if expected else "NOT_PROVIDED",
            }
        )
    return result, sources


def audit_ids(readsets, stage):
    """Use disk-backed keys rather than keeping every long-read ID in memory."""
    database = stage / "read_ids.sqlite"
    with sqlite3.connect(database) as connection:
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute(
            "CREATE TABLE reads (technology TEXT, id TEXT, PRIMARY KEY(technology,id)) WITHOUT ROWID"
        )
        for record in readsets.values():
            if record["technology"] not in ("hifi", "ont"):
                continue
            with pysam.FastxFile(str(stage / record["files"][0]["file"])) as reads:
                for read in reads:
                    try:
                        connection.execute(
                            "INSERT INTO reads VALUES (?, ?)", (record["technology"], read.name)
                        )
                    except sqlite3.IntegrityError as error:
                        raise ValueError(f'Duplicate {record["technology"]} read ID: {read.name}') from error
            connection.commit()
    database.unlink()


def prepare_sample(args):
    metadata = load(args.metadata)
    destination = Path(args.out_dir).absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with assembly_scratch(destination) as temporary:
        stage = Path(temporary) / "result"
        reads_dir = stage / "reads"
        reads_dir.mkdir(parents=True)
        readsets = {}
        for rid, source in metadata["readsets"].items():
            spec = dict(source)
            source_check(spec)
            if spec["technology"] in ("hic", "wgs"):
                files, provenance = paired_files(spec, reads_dir)
            else:
                spec_path = Path(temporary) / f"{rid}.json"
                save(spec_path, spec)
                fastq = reads_dir / f"{rid}.fastq.gz"
                manifest = Path(temporary) / f"{rid}.prepared.json"
                call = SimpleNamespace(spec=spec_path, fastq=fastq, manifest=manifest, threads=args.threads)
                (reuse if spec.get("prepared_reuse") else prepare)(call)
                files = [fastq]
                provenance = load(manifest)
                # Paths in the embedded preparation record must survive publication.
                provenance["prepared_fastq"] = str(destination / fastq.relative_to(stage))
            stats = fastq_statistics(files)
            source_check(spec)
            readsets[rid] = {
                "readset_id": rid,
                "library_id": spec["library_id"],
                "technology": spec["technology"],
                "statistics": stats,
                "files": [artifact(path, stage) for path in files],
                "sources": provenance,
            }
        if metadata["settings"]["assembly"]["strategy"] == "polyploid":
            audit_ids(readsets, stage)
        write_table(
            stage / "reads.tsv",
            ["readset_id", "technology", "num_reads", "num_bases", "num_pairs"],
            [{**record, **record["statistics"]} for record in readsets.values()],
        )
        save(
            stage / "manifest.json",
            {
                "sample_id": metadata["sample_id"],
                "status": "PASS",
                "readsets": readsets,
                "long_read_ids": (
                    "CHECKED"
                    if metadata["settings"]["assembly"]["strategy"] == "polyploid"
                    else "NOT_REQUESTED"
                ),
            },
        )
        publish(stage, destination, temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--threads", type=int, default=2)
    prepare_sample(parser.parse_args())


if __name__ == "__main__":
    main()
