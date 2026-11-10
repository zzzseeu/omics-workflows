#!/usr/bin/env python3
"""Native component smoke checks on synthetic data; not biological validation."""
import importlib
import argparse
from pathlib import Path
import random
import subprocess
import sys
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

# Local scripts must be imported after adding their repository path.
# pylint: disable=wrong-import-position,wrong-import-order
from assembly.lib.assessment import mapping_metrics, quast_metrics, telomere_metrics
from assembly.lib.sequences import fasta_dict, fasta_statistics
from assembly.lib.io import checksum, load, save, signature

prepare_sample = importlib.import_module("assembly.scripts.01_prepare_inputs").prepare_sample
from assembly.lib.postassembly import artifact, assembly_bundle


def invoke(script, *arguments):
    subprocess.run(
        [sys.executable, str(ROOT / "assembly/scripts" / script), *map(str, arguments)], check=True
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    rng = random.Random(41391)
    genome = "".join(rng.choices("ACGT", k=100000))
    fastq = output / "raw.fastq"
    with fastq.open("w", encoding="utf-8") as handle:
        for index in range(400):
            start = rng.randrange(88001)
            read = genome[start : start + 12000]
            handle.write(f"@molecule{index}\n{read}\n+\n" + "I" * len(read) + "\n")
    source = {
        "sample_id": "synthetic",
        "readset_id": "hifi",
        "library_id": "hifi",
        "technology": "hifi",
        "format": "fastq",
        "file1": str(fastq),
        "source_signature": signature(fastq),
        "expected_md5": None,
    }
    spec = output / "prepare.json"
    save(
        spec,
        {
            "sample_id": "synthetic",
            "readsets": {"hifi": source},
            "settings": {"assembly": {"strategy": "single"}},
        },
    )
    prepare_sample(SimpleNamespace(metadata=spec, out_dir=output / "01_inputs", threads=4))
    prepared = load(output / "01_inputs/manifest.json")
    # Error-free synthetic fixture is declared explicitly; NanoPlot/adapter
    # filtering is tested separately and is not claimed by this smoke test.
    reads = output / "clean_fixture"
    reads.mkdir()
    (reads / "reads.fastq.gz").symlink_to(output / "01_inputs/reads/hifi.fastq.gz")
    save(
        reads / "manifest.json",
        {
            "sample_id": "synthetic",
            "status": "PASS",
            "readsets": {
                "hifi": {
                    "readset_id": "hifi",
                    "library_id": "hifi",
                    "technology": "hifi",
                    "statistics": prepared["readsets"]["hifi"]["statistics"],
                    "files": [artifact(reads / "reads.fastq.gz", reads)],
                }
            },
            "validation": "SYNTHETIC_CLEAN_FIXTURE_NOT_TECHNOLOGY_QC",
        },
    )
    meta = output / "assembly.json"
    save(
        meta,
        {
            "sample_id": "synthetic",
            "settings": {
                "ploidy": 2,
                "output_copies": 1,
                "assembly": {
                    "strategy": "single",
                    "mode": "hifi_only",
                    "bloom_filter_bits": 0,
                    "purge_level": 0,
                },
                "scaffolding": {"enabled": False},
            },
            "tools": {"hifiasm": {"path": str(Path(sys.prefix) / "bin/hifiasm")}},
            "code": {},
        },
    )
    contigs = output / "03_assembly/contigs"
    invoke(
        "03_assemble_genome.py",
        "--stage",
        "contigs",
        "--metadata",
        meta,
        "--reads-manifest",
        reads / "manifest.json",
        "--out-dir",
        contigs,
        "--threads",
        4,
    )
    sequences = fasta_dict(contigs / "primary.fa")
    if sum(map(len, sequences.values())) < 80000:
        raise AssertionError("Unexpectedly incomplete synthetic primary assembly")
    if list(contigs.glob("*.bin")):
        raise AssertionError("Assembler caches were not pruned")
    delivery = output / "delivery_fixture"
    delivery.mkdir()
    for label in ("chromosomes", "reference.chromosomes"):
        (delivery / (label + ".fa")).symlink_to(contigs / "primary.fa")
        subprocess.run(["samtools", "faidx", str(delivery / (label + ".fa"))], check=True)
    (delivery / "unplaced.fa").write_text("", encoding="utf-8")
    save(
        delivery / "manifest.json",
        {
            "sample_id": "synthetic",
            "ploidy": 2,
            "strategy": "single",
            "output_copies": 1,
            "expected_chromosomes_per_haplotype": len(sequences),
            "quality_level": "NATIVE_CONTIG_FIXTURE_NOT_HIC_SCAFFOLDED",
            "sets": {
                label: fasta_statistics(delivery / (label + ".fa"))
                for label in ("chromosomes", "reference.chromosomes", "unplaced")
            },
        },
    )
    finished = output / "04_finishing"
    invoke(
        "04_finish_chromosomes.py", "--assembly-manifest", delivery / "manifest.json", "--out-dir", finished
    )
    polished = output / "05_polishing"
    invoke(
        "05_polish_chromosomes.py",
        "--assembly-manifest",
        finished / "manifest.json",
        "--out-dir",
        polished,
        "--skip",
    )
    _, unchanged, groups = assembly_bundle(polished / "manifest.json")
    if unchanged != sequences or set(groups.values()) != {"reference"}:
        raise AssertionError("Disabled post-assembly stages changed sequence/copy identity")
    quast = quast_metrics(polished / "genome.fa", "reference", output / "quast", 4, 1)
    mapping = mapping_metrics(
        polished / "genome.fa",
        [output / "01_inputs/reads/hifi.fastq.gz"],
        sequences,
        output / "mapping",
        4,
        1,
        10000,
    )
    rows, telomeres, _ = telomere_metrics(
        polished / "genome.fa", output / "tidk", "TTTAGGG", 10000, 200, 1000
    )
    if len(rows) != 2 * len(sequences) or mapping["mapping_percent"] < 95 or "N50" not in quast:
        raise AssertionError("Incomplete native assessment components")
    save(
        output / "validation.json",
        {
            "status": "PASS",
            "assembly_bases": sum(map(len, sequences.values())),
            "mapping": mapping,
            "tidk": telomeres,
            "final_sha256": checksum(polished / "genome.fa"),
            "checks": [
                "prepare_inputs",
                "single_primary_assembly",
                "cache_cleanup",
                "reference_identity",
                "finishing_passthrough",
                "polishing_passthrough",
                "QUAST",
                "HiFi_mosdepth",
                "tidk_ends",
            ],
            "limitations": [
                "Synthetic component validation only",
                "No Hi-C scaffolding tested by this job",
                "NanoPlot/HiFiAdapterFilt/LAI/BUSCO/Merqury full pipeline not executed",
            ],
        },
    )


if __name__ == "__main__":
    main()
