"""Seven-stage contracts: copy identity, input integrity and final assessment state."""

import importlib

from contextlib import ExitStack
import csv
import gzip
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

# Local scripts must be imported after adding their repository path.
# pylint: disable=wrong-import-position,wrong-import-order
from common.workflow import script_sources

assess_assembly = importlib.import_module("assembly.scripts.06_assess_assembly")
from assembly.lib import assessment as assessment_tools
from assembly.runner import resolve_config, verify_cached_outputs, verify_source_locks
from assembly.lib.assembly import run_gphase
from assembly.lib.assessment import interval_length, ltr_fractions, number
from assembly.lib.sequences import fasta_statistics
from assembly.lib.io import checksum, load, save, signature

prepare_sample = importlib.import_module("assembly.scripts.01_prepare_inputs").prepare_sample
merge_wgs = importlib.import_module("assembly.scripts.05_polish_chromosomes").merge_wgs
quality_reads = importlib.import_module("assembly.scripts.02_quality_reads")
filter_long = quality_reads.filter_long
from assembly.lib.postassembly import artifact, assembly_bundle, fastq_statistics, passthrough, read_bundle


# Each method tests a distinct workflow contract.
# pylint: disable-next=too-many-public-methods
class SevenStageTests(unittest.TestCase):
    def setUp(self):
        # unittest cleanup owns this directory across the complete test method.
        # pylint: disable-next=consider-using-with
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.path = Path(temporary.name)

    def reads(self, name="reads.fastq", identifier="read", sequence="ACGT"):
        path = self.path / name
        text = f"@{identifier}\n{sequence}\n+\n" + "I" * len(sequence) + "\n"
        if name.endswith(".gz"):
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(text)
        else:
            path.write_text(text, encoding="utf-8")
        return path

    def test_plain_fastq_preparation_is_binary_safe_and_preserves_ids(self):
        source = self.reads()
        metadata = self.path / "metadata.json"
        save(
            metadata,
            {
                "sample_id": "S",
                "settings": {"assembly": {"strategy": "single"}},
                "readsets": {
                    "h": {
                        "sample_id": "S",
                        "readset_id": "h",
                        "library_id": "H",
                        "technology": "hifi",
                        "format": "fastq",
                        "file1": str(source),
                        "expected_md5": None,
                        "source_signature": signature(source),
                    }
                },
            },
        )
        output = self.path / "01_inputs"
        prepare_sample(SimpleNamespace(metadata=metadata, out_dir=output, threads=1))
        result = read_bundle(output / "manifest.json", "S")
        with gzip.open(result["readsets"]["h"]["paths"][0], "rt", encoding="utf-8") as handle:
            self.assertEqual(handle.read(), source.read_text())
        self.assertEqual(result["readsets"]["h"]["statistics"]["num_reads"], 1)

    def test_pair_order_and_unequal_counts_cannot_pass(self):
        first = self.reads("r1.fq", "a/1")
        second = self.reads("r2.fq", "b/2")
        with self.assertRaisesRegex(ValueError, "IDs/order"):
            fastq_statistics([first, second])
        second.write_text("", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "counts differ"):
            fastq_statistics([first, second])

    def test_mixed_plain_and_gzip_wgs_remains_valid_gzip(self):
        first, second = self.reads("first.fq"), self.reads("second.fq.gz", "second")
        merged = self.path / "merged.fq.gz"
        merge_wgs([first, second, first], merged)
        with gzip.open(second, "rt") as handle:
            expected = first.read_text() + handle.read() + first.read_text()
        with gzip.open(merged, "rt") as handle:
            self.assertEqual(handle.read(), expected)

    def test_adapter_filter_never_mutates_the_prepared_input(self):
        source = self.reads("prepared.fq.gz")
        before = source.read_bytes()
        work, output, reports = [self.path / name for name in ("work", "output", "reports")]
        for folder in (work, output, reports):
            folder.mkdir()

        def native(command, **options):
            if command[0] == "hifiadapterfilt.sh":
                native_input = options["cwd"] / "reads.fastq.gz"
                self.assertFalse(native_input.is_symlink())
                native_input.unlink()
                (options["cwd"] / "reads.filt.fastq.gz").write_bytes(before)

        with patch.object(quality_reads, "execute", native):
            files = filter_long(
                {"paths": [source], "readset_id": "h", "technology": "hifi"},
                work,
                output,
                reports,
                {"hifi_adapter_length": 44, "hifi_adapter_identity": 97},
                4,
            )
        self.assertEqual(source.read_bytes(), before)
        self.assertEqual(files[0].read_bytes(), before)

    def test_mosdepth_coverage_uses_named_threshold_not_region_label(self):
        reference = self.path / "genome.fa"
        reference.write_text(">Chr01\nACGT\n")

        def native(_command):
            folder = self.path / "mapping"
            with gzip.open(folder / "coverage.thresholds.bed.gz", "wt") as handle:
                handle.write("#chrom\tstart\tend\tregion\t1X\nChr01\t0\t4\tunknown\t3\n")
            (folder / "coverage.mosdepth.summary.txt").write_text(
                "chrom\tlength\tbases\tmean\tmin\tmax\nChr01\t4\t8\t2\t0\t4\n"
            )

        alignment = SimpleNamespace(
            is_secondary=False, is_supplementary=False, is_unmapped=False, has_tag=lambda _tag: False
        )
        with ExitStack() as stack:
            stack.enter_context(patch.object(assessment_tools, "map_reads"))
            stack.enter_context(patch.object(assessment_tools, "execute", native))
            mocked = stack.enter_context(patch.object(assessment_tools.pysam, "AlignmentFile"))
            mocked.return_value.__enter__.return_value.fetch.return_value = [alignment]
            result = assessment_tools.mapping_metrics(
                reference, [], {"Chr01": "ACGT"}, self.path / "mapping", 4, 1, 10000
            )
        self.assertEqual(result["covered_fraction"], 0.75)
        self.assertEqual(result["mapping_percent"], 100)

    def test_mutated_clean_reads_and_wrong_sample_are_rejected(self):
        source = self.reads("clean.fq.gz")
        manifest = self.path / "manifest.json"
        save(
            manifest,
            {"sample_id": "S", "status": "PASS", "readsets": {"h": {"files": [artifact(source, self.path)]}}},
        )
        with self.assertRaisesRegex(ValueError, "different sample"):
            read_bundle(manifest, "OTHER")
        source.write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            read_bundle(manifest, "S")

    def test_directory_timestamp_cannot_hide_deleted_reads(self):
        folder = self.path / "samples/S/02_read_qc"
        folder.mkdir(parents=True)
        (folder / ".snakemake_timestamp").touch()
        source = self.reads("clean.fq.gz")
        save(
            folder / "manifest.json",
            {
                "sample_id": "S",
                "status": "PASS",
                "readsets": {"h": {"files": [{"file": "missing.fq.gz", "sha256": checksum(source)}]}},
            },
        )
        with self.assertRaisesRegex(ValueError, "Invalid cached output"):
            verify_cached_outputs({"sample_order": ["S"], "run_dir": str(self.path)})

    def test_changed_gphase_source_is_rejected_before_launch(self):
        root = self.path / "opt/GPhase"
        (root / "source").mkdir(parents=True)
        (root / "bin").mkdir()
        source = root / "source/gphase.py"
        source.write_text("original")
        save(root / "SOURCE_LOCK.json", {"source_files": {"gphase.py": checksum(source)}})
        tools = {"gphase": {"path": str(root / "bin/gphase")}}
        verify_source_locks(self.path, tools)
        self.assertIn("source_lock_sha256", tools["gphase"])
        source.write_text("changed")
        with self.assertRaisesRegex(ValueError, "Locked GPhase source changed"):
            verify_source_locks(self.path, tools)

    def configure(self, strategy, ploidy):
        (self.path / "samples.tsv").write_text("sample_id\tspecies\tprofile\nS\tOryza\t\n")
        lines = ["sample_id\treadset_id\tlibrary_id\ttechnology\tformat\tfile1\tfile2"]
        technologies = ["hifi", "hic", "ont"] if strategy == "polyploid" else ["hifi", "hic"]
        for technology in technologies:
            first = self.reads(technology + ".fq", technology + "/1")
            mate = self.reads("hic2.fq", "hic/2") if technology == "hic" else None
            lines.append(f'S\t{technology}\t{technology}\t{technology}\tfastq\t{first}\t{mate or ""}')
        (self.path / "readsets.tsv").write_text("\n".join(lines) + "\n")
        config = {
            "schema_version": 1,
            "samples_table": "samples.tsv",
            "readsets_table": "readsets.tsv",
            "defaults": {
                "ploidy": ploidy,
                "chromosomes_per_haplotype": 12,
                "assembly": {
                    "strategy": strategy,
                    "mode": "hifi_ul_hic" if strategy == "polyploid" else "hifi_hic",
                    "bloom_filter_bits": 0,
                },
                "scaffolding": {"enabled": True},
            },
        }
        path = self.path / "config.yaml"
        path.write_text(yaml.safe_dump(config))
        return resolve_config(path, ROOT, self.path / "run", None, 8, 64000)

    def test_single_output_does_not_change_biological_ploidy(self):
        settings = self.configure("single", 2)["samples"]["S"]["settings"]
        self.assertEqual(settings["ploidy"], 2)
        self.assertEqual(settings["output_copies"], 1)
        self.assertEqual(settings["assembly"]["representation"], "primary_alternate")

    def test_polyploid_copy_counts_are_parameterized(self):
        for count in (3, 4, 6):
            with self.subTest(ploidy=count):
                settings = self.configure("polyploid", count)["samples"]["S"]["settings"]
                self.assertEqual(settings["output_copies"], count)
                self.assertEqual(settings["scaffolding"]["backend"], "gphase")

    def test_final_assessment_is_independent_of_finishing(self):
        settings = self.configure("diploid", 2)["samples"]["S"]["settings"]
        self.assertFalse(settings["finishing"]["enabled"])
        self.assertTrue(settings["assessment"]["enabled"])

    def delivery(self):
        folder = self.path / "delivery"
        folder.mkdir()
        sequence = "CCCTAAA" * 40 + "ACGT" * 20 + "TTTAGGG" * 40
        for label in ("chromosomes", "reference.chromosomes"):
            (folder / f"{label}.fa").write_text(">Chr01\n" + sequence + "\n")
            (folder / f"{label}.fa.fai").write_text(
                f"Chr01\t{len(sequence)}\t7\t{len(sequence)}\t{len(sequence)+1}\n"
            )
        (folder / "unplaced.fa").write_text("")
        save(
            folder / "manifest.json",
            {
                "sample_id": "S",
                "strategy": "single",
                "ploidy": 2,
                "output_copies": 1,
                "expected_chromosomes_per_haplotype": 1,
                "sets": {
                    label: fasta_statistics(folder / f"{label}.fa")
                    for label in ("chromosomes", "reference.chromosomes", "unplaced")
                },
            },
        )
        return folder / "manifest.json", sequence

    def test_disabled_stages_link_without_duplicate_sequence_bytes(self):
        manifest, sequence = self.delivery()
        first, second = self.path / "04_finishing", self.path / "05_polishing"
        first.mkdir()
        passthrough(manifest, first, "DISABLED")
        second.mkdir()
        passthrough(first / "manifest.json", second, "DISABLED")
        self.assertTrue((second / "genome.fa").is_symlink())
        record, genome, groups = assembly_bundle(second / "manifest.json")
        self.assertEqual(genome, {"Chr01": sequence})
        self.assertEqual(groups, {"Chr01": "reference"})
        self.assertEqual(record["ploidy"], 2)
        self.assertEqual(record["output_copies"], 1)

    def test_unplaced_changes_are_detected(self):
        manifest, _ = self.delivery()
        (manifest.parent / "unplaced.fa").write_text(">extra\nACGT\n")
        with self.assertRaisesRegex(ValueError, "Unplaced FASTA checksum"):
            assembly_bundle(manifest)

    def assessment_arguments(self, manifest):
        source = self.reads("clean.fastq.gz")
        reads = self.path / "reads.json"
        save(
            reads,
            {
                "sample_id": "S",
                "status": "PASS",
                "readsets": {
                    "h": {
                        "readset_id": "h",
                        "library_id": "H",
                        "technology": "hifi",
                        "files": [artifact(source, self.path)],
                        "statistics": {"num_reads": 1},
                    }
                },
            },
        )
        lineage = self.path / "lineage"
        lineage.mkdir()
        (lineage / "dataset.cfg").write_text("name=fixture")
        return SimpleNamespace(
            assembly_manifest=str(manifest),
            reads_manifest=str(reads),
            busco_lineage=str(lineage),
            metadata=None,
            out_dir=str(self.path / "06_assessment"),
            threads=4,
            mem_mb=64000,
            kmer_size=21,
            min_mapq=1,
            coverage_window=10000,
            quast_min_contig=1,
            terminal_window=10000,
            terminal_max_distance=1000,
            telomere_min_length=200,
            telomere_motif="TTTAGGG",
            keep_intermediates=False,
        )

    def assessment_mocks(self, stack, sequence):
        def quast(_fasta, label, folder, _threads, _minimum):
            folder.mkdir(parents=True)
            (folder / "report.tsv").write_text(f"Assembly\t{label}\nN50\t{len(sequence)}\n")
            return {"N50": {label: str(len(sequence))}, "L50": {label: "1"}}

        def busco(_fasta, label, folder, _lineage, _threads):
            folder.mkdir(parents=True)
            summary, table = folder / f"{label}.json", folder / f"{label}.tsv"
            summary.write_text("{}")
            table.write_text("fixture\tComplete\n")
            return dict(zip("CSDFM", [100, 100, 0, 0, 0])), summary, [table]

        def mapping(_reference, _reads, _genome, folder, *_options):
            folder.mkdir(parents=True)
            for name in (
                "chromosome_coverage.tsv",
                "coverage.regions.bed.gz",
                "coverage.mosdepth.global.dist.txt",
                "hifi.bam",
            ):
                (folder / name).write_bytes(b"fixture")
            return {
                "status": "COMPLETED",
                "mapping_percent": 100,
                "split_read_percent": 0,
                "covered_fraction": 1,
            }

        def structure(_source, _genome, _reads, _metadata, _work, output, _threads):
            output.mkdir()
            for name in ("contact_map.pdf", "separate_plots.pdf"):
                (output / name).write_bytes(b"%PDF fixture")
            return {
                "status": "COMPLETED",
                "anchored_non_N_fraction": 1,
                "unplaced_length": 0,
                "human_review": "NOT_REVIEWED",
            }

        def telomeres(_genome, folder, *_parameters):
            folder.mkdir()
            table = folder / "telomeres.ends.tsv"
            table.write_text("id\tside\tstatus\nChr01\tstart\tpresent\nChr01\tend\tpresent\n")
            with table.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            return rows, {"t2t": 1, "unit": "TTTAGGG"}, table

        for name, function in [
            ("quast_metrics", quast),
            ("busco_metrics", busco),
            ("mapping_metrics", mapping),
            ("structural_metrics", structure),
            ("telomere_metrics", telomeres),
        ]:
            stack.enter_context(patch.object(assess_assembly, name, function))
        stack.enter_context(
            patch.object(
                assess_assembly,
                "lai_metrics",
                return_value=({"status": "NOT_APPLICABLE", "reason": "fixture lacks LTR content"}, None),
            )
        )

    def test_final_contract_missing_wgs_is_not_pass_and_does_not_certify_t2t(self):
        manifest, sequence = self.delivery()
        args = self.assessment_arguments(manifest)
        with ExitStack() as stack:
            self.assessment_mocks(stack, sequence)
            assess_assembly.assess(args)
        folder = Path(args.out_dir)
        result = load(folder / "assessment.json")
        self.assertEqual(set(result["evaluations"]), set(assess_assembly.CATEGORIES))
        self.assertEqual(result["status"], "PARTIAL")
        self.assertEqual(result["evaluations"]["kmers"]["status"], "MISSING_DATA")
        self.assertEqual(result["T2T_status"], "NOT_CERTIFIED")
        self.assertEqual(result["phasing"], "NOT_APPLICABLE")
        self.assertFalse(list(folder.rglob("*.bam")))
        self.assertFalse((folder / "intermediates").exists())
        self.assertTrue((folder / "telomeres.tsv").is_file())

    def test_native_failure_does_not_publish_completed_assessment(self):
        manifest, sequence = self.delivery()
        args = self.assessment_arguments(manifest)
        with ExitStack() as stack:
            self.assessment_mocks(stack, sequence)
            stack.enter_context(
                patch.object(assess_assembly, "busco_metrics", side_effect=RuntimeError("failed native tool"))
            )
            with self.assertRaisesRegex(RuntimeError, "failed native tool"):
                assess_assembly.assess(args)
        self.assertFalse((Path(args.out_dir) / "assessment.json").exists())

    def test_native_telomere_flags_and_missing_end_are_preserved(self):
        rows = [
            {"id": "Chr01", "side": "start", "status": "wrong_strand"},
            {"id": "Chr01", "side": "end", "status": "not_terminal"},
        ]
        ends = assess_assembly.validate_telomeres(rows, {"Chr01": "ACGT"})
        self.assertEqual(ends[("Chr01", "start")]["status"], "wrong_strand")
        with self.assertRaisesRegex(ValueError, "both ends"):
            assess_assembly.validate_telomeres(rows[:1], {"Chr01": "ACGT"})
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            assess_assembly.validate_telomeres(rows + rows[:1], {"Chr01": "ACGT"})

    def test_repeat_intervals_are_not_double_counted(self):
        self.assertEqual(interval_length([(0, 20), (10, 30), (30, 40)]), 40)
        intact = self.path / "pass.list"
        annotation = self.path / "repeat.out"
        intact.write_text("seq000001:1..20\tpass\nseq000001:11..30\tpass\n")
        annotation.write_text(
            "100 0 0 0 seq000001 1 80 0 + family LTR/Gypsy\n"
            "100 0 0 0 seq000001 81 90 0 + family LTR/Gypsy\n"
        )
        values = ltr_fractions(intact, annotation, {"seq000001": 100})
        self.assertEqual(values, {"intact_fraction": 0.3, "total_ltr_fraction": 0.8})
        with self.assertRaisesRegex(ValueError, "NaN"):
            number("nan")

    def test_cli_help_works_for_all_unified_entries(self):
        for name in (
            "01_prepare_inputs.py",
            "02_quality_reads.py",
            "03_assemble_genome.py",
            "04_finish_chromosomes.py",
            "05_polish_chromosomes.py",
            "06_assess_assembly.py",
            "07_write_assembly_report.py",
        ):
            result = subprocess.run(
                [sys.executable, str(ROOT / "assembly/scripts" / name), "--help"],
                capture_output=True,
                text=True,
                check=False,
                cwd=self.path,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("--out-dir", result.stdout)

    def test_exactly_seven_numbered_entries_and_shared_dependency_boundaries(self):
        entries = sorted((ROOT / "assembly/scripts").glob("[0-9][0-9]_*.py"))
        self.assertEqual([path.name[:3] for path in entries], [f"{i:02d}_" for i in range(1, 8)])
        assembler = set(script_sources(ROOT, "assembly/scripts/03_assemble_genome.py"))
        assessment = set(script_sources(ROOT, "assembly/scripts/06_assess_assembly.py"))
        self.assertIn(ROOT / "assembly/lib/assembly.py", assembler)
        self.assertIn(ROOT / "assembly/lib/scaffolding.py", assembler)
        self.assertIn(ROOT / "assembly/lib/__init__.py", assembler)
        self.assertNotIn(ROOT / "assembly/lib/assessment.py", assembler)
        self.assertIn(ROOT / "assembly/lib/assessment.py", assessment)
        self.assertNotIn(ROOT / "assembly/lib/assembly.py", assessment)

    def test_gphase_shared_function_preserves_ploidy_and_published_links(self):
        calls = []
        for copies in (3, 4, 6):
            with self.subTest(ploidy=copies):
                folder = self.path / str(copies)
                folder.mkdir()
                fasta, graph = folder / "unitigs.fa", folder / "unitigs.gfa"
                fasta.write_text(">u\nACGT\n")
                graph.write_text("S\tu\tACGT\n")
                (folder / "copies.tsv").write_text("u\t1\n")
                metadata, assembly = folder / "metadata.json", folder / "assembly.json"
                save(
                    metadata,
                    {
                        "settings": {
                            "ploidy": copies,
                            "chromosomes_per_haplotype": 12,
                            "scaffolding": {"processes": 4, "restriction_sites": "GATC"},
                        }
                    },
                )
                save(
                    assembly,
                    {
                        "assembly_id": "fixture",
                        "fasta_sha256": checksum(fasta),
                        "gfa_sha256": checksum(graph),
                    },
                )
                calls.clear()

                def native(command, cwd):
                    calls.append(command)
                    scaffold = cwd / "gphase_output/scaffold_hap"
                    scaffold.mkdir(parents=True)
                    for name in (
                        "unitig.fasta",
                        "unitig.scaffold.agp",
                        "unitig.scaffold.fasta",
                        "contig.fasta",
                        "contig.scaffold.agp",
                        "contig.scaffold.fasta",
                        "ctg2utg.txt",
                    ):
                        (scaffold / ("gphase_final." + name)).write_text("fixture\n")
                    (cwd / "internal_link").symlink_to(scaffold / "gphase_final.unitig.fasta")

                output = folder / "gphase"
                with patch("assembly.lib.assembly.execute", native):
                    run_gphase(
                        SimpleNamespace(
                            metadata=metadata,
                            assembly=assembly,
                            fasta=fasta,
                            gfa=graph,
                            copies=folder / "copies.tsv",
                            bam=folder / "hic.bam",
                            out_dir=output,
                            threads=8,
                        )
                    )
                command = calls[0]
                self.assertEqual(command[command.index("--n_hap") + 1], copies)
                self.assertEqual(command[command.index("--processes") + 1], 4)
                self.assertEqual((output / "internal_link").read_text(), "fixture\n")
                self.assertEqual(load(output / "manifest.json")["human_review"], "NOT_REVIEWED")


if __name__ == "__main__":
    unittest.main()
