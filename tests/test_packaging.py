"""Installed command, module registration, resources and environment handoff contracts."""

import importlib
import io
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

from omics_workflows import __version__, cli
from omics_workflows.common.arguments import argument_parser
from omics_workflows.common.resources import package_directory, source_boundary
from omics_workflows.common.workflow import script_sources
from omics_workflows.assembly import runner


class PackagingTests(unittest.TestCase):
    @staticmethod
    def inspect_command_help(module):
        parsers = []

        def capture(*args, **kwargs):
            parser = argument_parser(*args, **kwargs)
            parsers.append(parser)
            return parser

        output = io.StringIO()
        with (
            patch.object(module, "argument_parser", side_effect=capture),
            patch.object(sys, "argv", [module.__name__, "--help"]),
            redirect_stdout(output),
        ):
            try:
                module.main()
            except SystemExit as status:
                if status.code != 0:
                    raise
        return parsers, output.getvalue()

    def test_every_installed_command_parameter_has_a_description_and_type(self):
        modules = [cli, runner]
        modules.extend(
            importlib.import_module(f"omics_workflows.assembly.scripts.{name}")
            for name in (
                "01_prepare_inputs",
                "02_quality_reads",
                "03_assemble_genome",
                "04_finish_chromosomes",
                "05_polish_chromosomes",
                "06_assess_assembly",
                "07_write_assembly_report",
            )
        )
        modules.append(importlib.import_module("omics_workflows.assembly.examples.prepare_demo"))
        for module in modules:
            with self.subTest(command=module.__name__):
                parsers, output = self.inspect_command_help(module)
                self.assertEqual(len(parsers), 1)
                # argparse exposes no public iterator for registered actions.
                # pylint: disable-next=protected-access
                for action in parsers[0]._actions:
                    self.assertIn("类型：", action.help, action.dest)
                    if action.required:
                        self.assertIn("必填", action.help, action.dest)
                    for option in action.option_strings:
                        self.assertIn(option, output)

    def test_assembly_help_changes_preserve_argument_types_and_defaults(self):
        with patch.object(runner, "execute", return_value=0) as execute:
            status = runner.main(
                [
                    "--config",
                    "config/assembly.yaml",
                    "--run-dir",
                    "runs/demo",
                    "--cores",
                    "8",
                    "--mem-mb",
                    "64000",
                    "--samples",
                    "rice01",
                    "rice02",
                    "--check",
                ]
            )
        self.assertEqual(status, 0)
        args = execute.call_args.args[0]
        self.assertEqual(args.config, Path("config/assembly.yaml"))
        self.assertEqual(args.run_dir, Path("runs/demo"))
        self.assertEqual((args.cores, args.mem_mb), (8, 64000))
        self.assertEqual(args.samples, ["rice01", "rice02"])
        self.assertEqual(args.target, "report")
        self.assertTrue(args.check)
        self.assertFalse(args.dry_run)

    def test_version_and_list_do_not_import_analysis_modules(self):
        output = io.StringIO()
        with patch.object(cli.importlib, "import_module") as load, redirect_stdout(output):
            self.assertEqual(cli.main(["list"]), 0)
            with self.assertRaises(SystemExit) as status:
                cli.main(["--version"])
        self.assertEqual(status.exception.code, 0)
        self.assertIn("assembly", output.getvalue())
        self.assertIn(__version__, output.getvalue())
        load.assert_not_called()

    def test_another_registered_workflow_receives_its_own_arguments(self):
        entry = Mock(return_value=7)
        with (
            patch.dict(cli.WORKFLOWS, {"expression": ("example.expression", "Expression analysis")}),
            patch.object(cli.importlib, "import_module", return_value=SimpleNamespace(main=entry)) as load,
        ):
            result = cli.main(["expression", "--design", "design.tsv"])
        self.assertEqual(result, 7)
        load.assert_called_once_with("example.expression")
        entry.assert_called_once_with(["--design", "design.tsv"], prog="omics-workflows expression")

    def test_all_runtime_rules_and_templates_are_package_resources(self):
        with package_directory() as package:
            assembly = package / "assembly"
            for name in ("Snakefile", "defaults.yaml"):
                self.assertTrue((assembly / name).is_file(), name)
            self.assertFalse((assembly / "schema.yaml").exists())
            self.assertEqual(len(list((assembly / "rules").glob("*.smk"))), 7)
            for preset in ("single_hic", "hifi_only", "hifi_ul_hic", "polyploid"):
                for name in ("config.yaml", "samples.tsv", "readsets.tsv"):
                    self.assertTrue((assembly / "examples" / preset / name).is_file(), (preset, name))

    def test_editable_source_boundary_protects_the_whole_checkout(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            package = root / "example"
            package.mkdir(parents=True)
            (root / "pyproject.toml").touch()
            self.assertEqual(source_boundary(package), root)
            installed = root / "site-packages/example"
            self.assertEqual(source_boundary(installed), installed)

    def test_namespaced_source_tracking_stays_inside_the_package(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            package = root / "example"
            package.mkdir()
            (package / "__init__.py").touch()
            entry = package / "entry.py"
            entry.write_text("from example.helper import VALUE\nimport external\n", encoding="utf-8")
            helper = package / "helper.py"
            helper.write_text("VALUE = 1\n", encoding="utf-8")
            external = root / "external.py"
            external.write_text("raise RuntimeError('external source')\n", encoding="utf-8")
            sources = script_sources(package, entry)
            self.assertIn(helper, sources)
            self.assertNotIn(external, sources)

    def test_environment_handoff_uses_the_active_installation_without_a_checkout_entry(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory)
            python = prefix / "bin/python"
            python.parent.mkdir()
            python.touch()
            arguments = [
                "--config",
                "config.yaml",
                "--run-dir",
                "run001",
                "--cores",
                "2",
                "--mem-mb",
                "4096",
                "--env-prefix",
                str(prefix),
                "--check",
            ]
            with patch.object(
                runner.os, "execve", side_effect=RuntimeError("environment handoff")
            ) as execute:
                with self.assertRaisesRegex(RuntimeError, "environment handoff"):
                    runner.main(arguments)
            executable, command, environment = execute.call_args.args
            self.assertEqual(executable, str(python))
            self.assertEqual(command[1:4], ["-m", "omics_workflows.cli", "assembly"])
            self.assertEqual(command[4:], arguments)
            with package_directory() as package:
                self.assertEqual(environment["PYTHONPATH"].split(os.pathsep)[0], str(package.parent))
