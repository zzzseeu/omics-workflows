"""Shared workflow utilities must work without genome-assembly definitions."""

import ast
import fcntl
import json
import signal
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from omics_workflows.common import workflow


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        # unittest cleanup spans the complete test method.
        # pylint: disable-next=consider-using-with
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def test_shared_module_does_not_import_an_analysis_module(self):
        tree = ast.parse(Path(workflow.__file__).read_text(encoding="utf-8"))
        imports = [node.module for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)]
        self.assertFalse(any(name and name.startswith("omics_workflows.workflows.") for name in imports))

    def test_configuration_uses_the_supplied_python_validator_and_defaults(self):
        config = self.root / "config.yaml"
        defaults = self.root / "defaults.yaml"

        def validate(values):
            """由示例模块直接声明允许的字段和值，不加载外部校验文件。"""
            workflow.check_fields(values, {"assay", "options"}, required={"assay"})
            workflow.check_value(values["assay"], str, "assay", choices=("expression",))
            if "options" in values:
                workflow.check_fields(values["options"], {"threads", "paired"}, path="options")

        defaults.write_text("options: {threads: 2, paired: true}\n", encoding="utf-8")
        config.write_text("assay: expression\noptions: {threads: 4}\n", encoding="utf-8")
        merged = workflow.merge(
            workflow.load_defaults(defaults), workflow.load_configuration(config, validate)
        )
        self.assertEqual(merged["options"], {"threads": 4, "paired": True})
        self.assertEqual(workflow.load_defaults(defaults)["options"]["threads"], 2)
        config.write_text("assay: unsupported\n", encoding="utf-8")
        with self.assertRaisesRegex(workflow.ConfigError, "Configuration validation failed"):
            workflow.load_configuration(config, validate)

    def test_yaml_root_and_parse_errors_are_reported_before_the_module_validator(self):
        config = self.root / "config.yaml"
        validate = Mock()
        for value in ("", "[]\n", "null\n", "input: [\n"):
            with self.subTest(value=value):
                config.write_text(value, encoding="utf-8")
                with self.assertRaisesRegex(workflow.ConfigError, "Configuration validation failed"):
                    workflow.load_configuration(config, validate)
        validate.assert_not_called()

    def test_legacy_normalization_runs_before_the_python_validator(self):
        config = self.root / "config.yaml"
        config.write_text("legacy: value\n", encoding="utf-8")
        normalize = Mock(return_value={"current": "value"})
        validate = Mock()
        result = workflow.load_configuration(config, validate, normalize)
        normalize.assert_called_once_with({"legacy": "value"})
        validate.assert_called_once_with({"current": "value"})
        self.assertEqual(result, {"current": "value"})

    def test_field_checks_reject_wrong_structure_keys_and_missing_fields(self):
        for value in (None, [], {1: "value"}, {"unknown": 1}):
            with self.subTest(value=value), self.assertRaises(workflow.ConfigError):
                workflow.check_fields(value, {"required"}, required={"required"})

    def test_value_checks_do_not_coerce_types_and_reject_nonfinite_numbers(self):
        for value in (True, False, "8", 8.0):
            with self.subTest(value=value), self.assertRaises(workflow.ConfigError):
                workflow.check_value(value, int, "threads", minimum=1)
        for value in (float("nan"), float("inf"), -float("inf")):
            with self.subTest(value=value), self.assertRaises(workflow.ConfigError):
                workflow.check_value(value, (int, float), "quality", minimum=0)
        workflow.check_value(None, (int, type(None)), "optional_count", minimum=1)
        workflow.check_value(7.5, (int, float), "quality", minimum=0)

    def test_sample_selection_and_technology_values_are_supplied_by_caller(self):
        samples, selected = workflow.select_samples([{"sample_id": "A"}, {"sample_id": "B"}], ["B"])
        self.assertEqual(selected, ["B"])
        row = {"sample_id": "B", "readset_id": "r", "library_id": "L", "technology": "rna"}
        self.assertEqual(workflow.group_readsets([row], samples, {"rna"})["B"], [row])
        with self.assertRaisesRegex(workflow.ConfigError, "unknown technology"):
            workflow.group_readsets([row], samples, {"protein"})

    def test_tools_can_be_checked_in_an_independent_environment(self):
        prefix = self.root / "environment"
        binary = prefix / "bin/tool"
        binary.parent.mkdir(parents=True)
        binary.write_text("#!/bin/sh\necho example-tool-1\n", encoding="utf-8")
        binary.chmod(0o755)
        with patch.dict("os.environ", {"PATH": str(binary.parent)}):
            tools, missing = workflow.inspect_tools({"tool", "missing"}, prefix, {"tool": ["--version"]})
        self.assertEqual(missing, ["missing"])
        self.assertEqual(tools["tool"]["version"], "example-tool-1")
        self.assertEqual(tools["tool"]["sha256"], workflow.file_checksum(binary))
        with self.assertRaisesRegex(workflow.ConfigError, "missing"):
            workflow.require_available(missing)

    def test_extension_lock_layout_is_supplied_by_caller(self):
        root = self.root / "extension"
        root.mkdir()
        source = root / "helper.py"
        source.write_text("original", encoding="utf-8")
        workflow.write_stable(
            root / "SOURCE_LOCK.json", {"content": {"helper.py": workflow.file_checksum(source)}}
        )
        tools = {"custom": {"path": str(root / "entry")}}
        specification = [{"root": root, "files_key": "content", "label": "example"}]
        workflow.verify_source_locks(tools, specification)
        self.assertIn("source_lock_sha256", tools["custom"])
        source.write_text("changed", encoding="utf-8")
        with self.assertRaisesRegex(workflow.ConfigError, "Locked example source changed"):
            workflow.verify_source_locks(tools, specification)

    def test_stage_records_track_only_their_tools_and_local_source_dependencies(self):
        repo = self.root / "source"
        package = repo / "otherflow"
        package.mkdir(parents=True)
        (package / "__init__.py").write_text("", encoding="utf-8")
        (package / "first.py").write_text("from otherflow.helper import VALUE\n", encoding="utf-8")
        (package / "helper.py").write_text("VALUE = 1\n", encoding="utf-8")
        (package / "second.py").write_text("VALUE = 2\n", encoding="utf-8")
        (package / "first.smk").write_text("# first\n", encoding="utf-8")
        (package / "second.smk").write_text("# second\n", encoding="utf-8")
        plan = {
            "tools": {"reader": {"version": "1"}, "writer": {"version": "2"}},
            "samples": {"A": {"sample_id": "A"}},
        }
        stages = {
            "first": {"rule": "otherflow/first.smk", "entries": ["otherflow/first.py"], "tools": {"reader"}},
            "second": {"rule": "otherflow/second.smk", "entries": ["otherflow/second.py"], "tools": set()},
        }

        def metadata(sample, _stage, tools, code):
            return {"sample_id": sample["sample_id"], "tools": tools, "code": code}

        output = self.root / "run"
        workflow.write_plan_metadata(plan, repo, output, stages, metadata)
        first = output / "config/samples/A.first.json"
        second = output / "config/samples/A.second.json"
        before = json.loads(first.read_text(encoding="utf-8"))
        second_before = workflow.source_signature(second)
        self.assertEqual(set(before["tools"]), {"reader"})
        self.assertIn("otherflow/helper.py", before["code"])
        self.assertEqual(json.loads(second.read_text(encoding="utf-8"))["tools"], {})
        (package / "helper.py").write_text("VALUE = 3\n", encoding="utf-8")
        workflow.write_plan_metadata(plan, repo, output, stages, metadata)
        self.assertNotEqual(json.loads(first.read_text(encoding="utf-8"))["code"], before["code"])
        self.assertEqual(workflow.source_signature(second), second_before)

    def test_cached_manifest_uses_the_supplied_result_validator(self):
        folder = self.root / "run/samples/A/01_matrix"
        folder.mkdir(parents=True)
        (folder / ".snakemake_timestamp").touch()
        workflow.write_stable(folder / "matrix.json", {"sample_id": "A", "format": "counts"})
        validator = Mock()
        plan = {"run_dir": str(self.root / "run"), "sample_order": ["A"]}
        workflow.verify_cached_outputs(plan, {"01_matrix": ("matrix.json", validator)})
        validator.assert_called_once_with(folder / "matrix.json", "A")
        workflow.write_stable(folder / "matrix.json", {"sample_id": "WRONG"})
        with self.assertRaisesRegex(workflow.ConfigError, "another sample"):
            workflow.verify_cached_outputs(plan, {"01_matrix": ("matrix.json", validator)})

    def plan(self):
        return {
            "run_dir": str(self.root / "run"),
            "sample_order": ["A"],
            "plan_id": "example",
            "target": "counts",
            "budget": {"cores": 2, "mem_mb": 4096},
        }

    def test_executor_accepts_another_workflow_target_and_records_planning(self):
        process = Mock(pid=12345)
        process.wait.return_value = 0
        prepare = Mock()
        with patch.object(workflow.subprocess, "Popen", return_value=process) as launch:
            code = workflow.execute_workflow(
                self.plan(), self.root / "otherflow/Snakefile", "counts", prepare, True
            )
        self.assertEqual(code, 0)
        prepare.assert_called_once()
        command = launch.call_args.args[0]
        self.assertIn("counts", command)
        self.assertIn("--dry-run", command)
        self.assertIn("mem_mb=4096", command)
        status = json.loads((self.root / "run/dryrun_status.json").read_text(encoding="utf-8"))
        self.assertEqual(status["status"], "PLANNED")

    def test_nonzero_workflow_exit_is_not_reported_as_completed(self):
        process = Mock(pid=12345)
        process.wait.return_value = 7
        with patch.object(workflow.subprocess, "Popen", return_value=process):
            self.assertEqual(workflow.execute_workflow(self.plan(), self.root / "Snakefile", "counts"), 7)
        status = json.loads((self.root / "run/run_status.json").read_text(encoding="utf-8"))
        self.assertEqual((status["status"], status["exit_code"]), ("FAILED", 7))

    def test_cancellation_reaches_the_child_process_group(self):
        process = Mock(pid=12345)
        handlers = {}

        def install(signum, handler):
            previous = handlers.get(signum)
            handlers[signum] = handler
            return previous

        def wait(timeout):
            self.assertEqual(timeout, 1)
            handlers[signal.SIGTERM](signal.SIGTERM, None)
            return 0

        process.wait.side_effect = wait
        with (
            patch.object(workflow.subprocess, "Popen", return_value=process),
            patch.object(workflow.signal, "signal", side_effect=install),
            patch.object(workflow.os, "killpg") as kill,
        ):
            code = workflow.execute_workflow(self.plan(), self.root / "Snakefile", "counts")
        self.assertEqual(code, 128 + signal.SIGTERM)
        kill.assert_called_once_with(12345, signal.SIGTERM)

    def test_existing_launcher_lock_prevents_a_second_run(self):
        folder = self.root / "run"
        folder.mkdir()
        with (folder / ".launcher.lock").open("a", encoding="utf-8") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with self.assertRaisesRegex(workflow.ConfigError, "active launcher"):
                workflow.execute_workflow(self.plan(), self.root / "Snakefile", "counts")
