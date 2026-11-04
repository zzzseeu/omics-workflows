"""Renaming must preserve sequence bytes, historical evidence and recoverability."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('layout_migration',
                                            ROOT / 'assembly/migrate_layout.py')
layout = importlib.util.module_from_spec(spec)
spec.loader.exec_module(layout)


class MigrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.run = Path(self.tmp.name) / 'run'
        for directory in ['config', 'samples/A/inputs', 'samples/A/qc/reads',
                          'samples/A/qc/assembly', 'samples/A/assembly/hifiasm',
                          'samples/A/report', 'reports']:
            (self.run / directory).mkdir(parents=True, exist_ok=True)
        layout.save(self.run / 'config/plan.json', {'sample_order': ['A']})
        layout.save(self.run / 'run_status.json', {'status': 'COMPLETED', 'exit_code': 0,
                                                 'command': ['historical-command']})
        self.fasta = self.run / 'samples/A/assembly/hifiasm/primary.fa'
        self.fasta.write_text('>a\nACGT\n')
        self.fasta.with_name('primary.fa.fai').write_text('a\t4\t3\t4\t5\n')
        self.fasta.with_name('primary.gfa').write_text('S\ta\tACGT\n')
        self.fasta.with_name('cache.bin').write_bytes(b'old-cache-do-not-change')
        layout.save(self.fasta.with_name('manifest.json'), {'assembly_id': 'original-id',
                    'command': [str(self.fasta)], 'input_manifests': [
                        {'prepared_fastq': 'samples/A/inputs/reads.fastq.gz'}]})
        report = {'assembly_id': 'original-id', 'assembly_manifest': str(self.fasta.with_name('manifest.json')),
                  'sets': [{'fasta': str(self.fasta), 'fai': str(self.fasta) + '.fai',
                            'sha256': layout.sha256(self.fasta),
                            'fai_sha256': layout.sha256(Path(str(self.fasta) + '.fai'))}]}
        layout.save(self.run / 'samples/A/report/manifest.json', report)
        (self.run / 'samples/A/report/report.html').write_text('historical report')
        layout.save(self.run / 'reports/batch.json', {'samples': [
            {'manifest': str(self.run / 'samples/A/report/manifest.json')}]})

    def test_plan_does_not_move_files(self):
        result = layout.migrate(self.run)
        self.assertEqual(result['status'], 'PLAN_ONLY')
        self.assertTrue(self.fasta.is_file())
        self.assertFalse((self.run / 'layout_migration').exists())

    def test_migration_preserves_and_updates_references(self):
        before = layout.stat_record(self.fasta)
        result = layout.migrate(self.run, apply=True)
        path = self.run / 'samples/A/04_assembly/hifiasm/primary.fa'
        self.assertEqual(layout.stat_record(path), before)
        self.assertEqual(result['status'], 'MIGRATED')
        manifest = layout.load(self.run / 'samples/A/10_report/manifest.json')
        self.assertEqual(manifest['sets'][0]['fasta'], str(path))
        native = layout.load(path.with_name('manifest.json'))
        self.assertEqual(native['command'], [str(self.fasta)])
        self.assertEqual(native['input_manifests'][0]['prepared_fastq'],
                         'samples/A/01_inputs/reads.fastq.gz')
        self.assertTrue((self.run / 'layout_migration/original_metadata/samples/A/report/manifest.json').exists())
        self.assertEqual(layout.migrate(self.run, apply=True)['status'], 'MIGRATED')

    def test_running_run_is_rejected(self):
        layout.save(self.run / 'run_status.json', {'status': 'RUNNING'})
        with self.assertRaisesRegex(ValueError, 'completed'):
            layout.migrate(self.run, apply=True)

    def test_resume_after_partial_move(self):
        original = Path.replace

        def interrupted(path, target):
            if path == self.run / 'samples/A/assembly':
                raise OSError('simulated interruption')
            return original(path, target)

        with patch.object(Path, 'replace', interrupted):
            with self.assertRaisesRegex(OSError, 'interruption'):
                layout.migrate(self.run, apply=True)
        self.assertEqual(layout.load(self.run / 'layout_migration/status.json')['status'], 'MIGRATING')
        self.assertEqual(layout.migrate(self.run, apply=True)['status'], 'MIGRATED')

    def test_target_collision_is_rejected(self):
        (self.run / 'samples/A/04_assembly').mkdir()
        with self.assertRaisesRegex(ValueError, 'already exists'):
            layout.migrate(self.run, apply=True)
        self.assertTrue(self.fasta.is_file())


if __name__ == '__main__':
    unittest.main()
