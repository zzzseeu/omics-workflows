"""Public command and project-boundary contracts for the simplified source layout."""

import os
import subprocess
import tempfile
import unittest
from pathlib import Path

import omics_workflows

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = Path(omics_workflows.__file__).resolve().parent


class InterfaceTests(unittest.TestCase):
    def invoke(self, *args):
        with tempfile.TemporaryDirectory() as directory:
            return subprocess.run(
                ["/usr/bin/python3", "-m", "omics_workflows.cli", *args],
                cwd=directory,
                env=dict(os.environ, PYTHONPATH=str(PACKAGE.parent), PYTHONDONTWRITEBYTECODE="1"),
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )

    def test_root_help_without_analysis_environment(self):
        result = self.invoke("--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("assembly", result.stdout)

    def test_module_help_without_analysis_environment(self):
        result = self.invoke("assembly", "--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("--env-prefix", result.stdout)
        self.assertIn("--target", result.stdout)

    def test_unknown_workflow_is_rejected(self):
        result = self.invoke("unknown_workflow")
        self.assertEqual(result.returncode, 2)
        self.assertIn("invalid choice", result.stderr)

    def test_direct_assembly_help_outside_repository_without_analysis_environment(self):
        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run(
                [
                    "/usr/bin/python3",
                    str(PACKAGE / "assembly/runner.py"),
                    "--help",
                ],
                cwd=directory,
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--config", result.stdout)

    def test_site_submission_scripts_are_project_owned(self):
        self.assertFalse(list(ROOT.rglob("*.csub.sh")))
        for path in ROOT.rglob("*.sh"):
            self.assertNotIn("#CSUB", path.read_text())
        self.assertFalse((ROOT / "launchers").exists())
        self.assertEqual(
            {
                p.name
                for p in ROOT.iterdir()
                if p.is_dir()
                and not p.name.startswith(".")
                and not p.name.endswith(".egg-info")
                and p.name not in {"build", "dist"}
            },
            {"omics_workflows", "env", "tests", "docs"},
        )


if __name__ == "__main__":
    unittest.main()
