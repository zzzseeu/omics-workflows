"""Configuration checks must reject invalid overrides before tools are started."""

import copy
import unittest

from omics_workflows.common.resources import workflow_defaults
from omics_workflows.common.workflow import ConfigError
from omics_workflows.assembly.runner import validate_config, validate_settings


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.config = {"samples_table": "samples.tsv", "readsets_table": "reads.tsv"}

    def test_required_fields_and_types_are_checked(self):
        for field, value in (
            ("samples_table", ""),
            ("samples_table", []),
            ("readsets_table", None),
        ):
            with self.subTest(field=field, value=value):
                config = {**self.config, field: value}
                with self.assertRaisesRegex(ConfigError, field):
                    validate_config(config)
        for field in self.config:
            config = dict(self.config)
            del config[field]
            with self.subTest(missing=field), self.assertRaisesRegex(ConfigError, field):
                validate_config(config)

    def test_all_override_scopes_are_checked_including_unused_profiles(self):
        for scope in ("defaults", "analysis_profiles", "samples"):
            with self.subTest(scope=scope):
                config = copy.deepcopy(self.config)
                settings = {"assembly": {"thread": 4}}
                config[scope] = settings if scope == "defaults" else {"unused": settings}
                with self.assertRaisesRegex(ConfigError, rf"{scope}.*assembly.*thread"):
                    validate_config(config)

    def test_invalid_types_and_bounds_are_rejected_with_field_paths(self):
        cases = (
            ("assembly", "threads", True),
            ("assembly", "threads", "8"),
            ("assembly", "threads", 8.0),
            ("assembly", "mem_mb", 511),
            ("assembly", "bloom_filter_bits", 40),
            ("assembly", "purge_level", -1),
            ("read_qc", "ont_min_quality", float("nan")),
            ("read_qc", "hifi_adapter_identity", 101),
            ("read_qc", "fastp_quality", -1),
            ("scaffolding", "enabled", "false"),
            ("scaffolding", "threads", 3),
            ("scaffolding", "restriction_sites", "gatc"),
            ("scaffolding", "enzyme_provenance", ""),
            ("scaffolding", "mapq", 61),
            ("finishing", "telomere_motif", "TTTAGGG\n"),
            ("finishing", "min_identity", 0.94),
            ("finishing", "min_support", 1),
            ("polishing", "on_missing_wgs", "ignore"),
            ("polishing", "min_contig_length", 999),
            ("assessment", "min_busco_complete", 101),
            ("assessment", "kmer_size", 32),
            ("assessment", "busco", 0),
            ("assessment", "busco_lineage", False),
            ("assessment", "terminal_window", 99),
            ("assessment", "terminal_min_repeats", 2),
        )
        for group, field, value in cases:
            with self.subTest(group=group, field=field, value=value):
                config = {**self.config, "defaults": {group: {field: value}}}
                with self.assertRaisesRegex(ConfigError, rf"defaults\.{group}\.{field}"):
                    validate_config(config)

    def test_only_implemented_data_modes_are_allowed(self):
        for mode in ("ont_only", "ont_hic", "hifi_ul", "unknown"):
            with self.subTest(mode=mode), self.assertRaisesRegex(ConfigError, "assembly.mode"):
                validate_config({**self.config, "defaults": {"assembly": {"mode": mode}}})

    def test_nullable_metadata_and_supported_parameter_boundaries_remain_valid(self):
        settings = {
            "ploidy": None,
            "genome_size_bp": None,
            "chromosomes_per_haplotype": None,
            "assembly": {"strategy": None, "bloom_filter_bits": 0, "purge_level": None},
            "read_qc": {"ont_min_quality": 7.5, "hifi_adapter_identity": 97.5},
            "scaffolding": {
                "restriction_sites": "GATC,GANTC",
                "review_agp": None,
                "review_manifest": None,
                "min_retained_fraction": 0,
            },
            "finishing": {"telomere_motif": None, "source_delivery": None, "min_identity": 1},
            "assessment": {"busco_lineage": None, "min_busco_complete": 100},
        }
        before = copy.deepcopy(settings)
        validate_config({**self.config, "defaults": settings})
        self.assertEqual(settings, before)

    def test_module_defaults_and_every_declared_section_pass_python_validation(self):
        import yaml

        defaults = yaml.safe_load(workflow_defaults("assembly"))
        validate_settings(defaults, "module_defaults")
        for group in ("assembly", "read_qc", "scaffolding", "finishing", "polishing", "assessment"):
            with self.subTest(group=group):
                with self.assertRaisesRegex(ConfigError, group):
                    validate_settings({group: None}, "defaults")
