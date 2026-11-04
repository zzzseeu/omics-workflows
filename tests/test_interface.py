"""Public command and project-boundary contracts for the simplified source layout."""
from pathlib import Path
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]


class InterfaceTests(unittest.TestCase):
    def invoke(self, *args):
        return subprocess.run(['/usr/bin/python3', str(ROOT/'run'), *args],
                              capture_output=True, text=True, timeout=30)

    def test_root_help_without_analysis_environment(self):
        result = self.invoke('--help')
        self.assertEqual(result.returncode, 0)
        self.assertIn('assembly', result.stdout)

    def test_module_help_without_analysis_environment(self):
        result = self.invoke('assembly', '--help')
        self.assertEqual(result.returncode, 0)
        self.assertIn('--env-prefix', result.stdout)
        self.assertIn('--target', result.stdout)

    def test_unknown_workflow_is_rejected(self):
        result = self.invoke('unknown_workflow')
        self.assertEqual(result.returncode, 2)
        self.assertIn('invalid choice', result.stderr)

    def test_site_submission_scripts_are_project_owned(self):
        self.assertFalse(list(ROOT.rglob('*.csub.sh')))
        for path in ROOT.rglob('*.sh'):
            self.assertNotIn('#CSUB', path.read_text())
        self.assertFalse((ROOT/'launchers').exists())
        self.assertEqual({p.name for p in ROOT.iterdir() if p.is_dir() and not p.name.startswith('.')},
                         {'assembly', 'common', 'env', 'tests', 'docs'})


if __name__ == '__main__':
    unittest.main()
