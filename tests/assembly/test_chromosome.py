"""Biologically material delivery gates: exact coverage, haplotypes and stale reviews."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'assembly/scripts'))
import agp_utils
import io_utils
from validate_hic_pairs import paired_input


class AgpTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)/'test.agp'
        self.sequences = {'a': 'ACGTA', 'b': 'GGTT'}

    def parse(self, rows):
        self.path.write_text(rows)
        return agp_utils.parse_agp(self.path, self.sequences)

    def test_corrected_split_reverse_complement_and_gap(self):
        rows = self.parse('chr1\t1\t2\t1\tW\ta\t1\t2\t+\n'
                          'chr1\t3\t5\t2\tU\t3\tscaffold\tyes\tproximity_ligation\n'
                          'chr1\t6\t8\t3\tW\ta\t3\t5\t-\n'
                          'unplaced\t1\t4\t1\tW\tb\t1\t4\t+\n')
        self.assertEqual(agp_utils.reconstruct(rows['chr1'], self.sequences), 'ACNNNTAC')

    def test_sequence_loss_fails(self):
        with self.assertRaisesRegex(ValueError, 'loses source'):
            self.parse('chr1\t1\t5\t1\tW\ta\t1\t5\t+\n')

    def test_duplicate_interval_fails(self):
        with self.assertRaisesRegex(ValueError, 'duplicates or loses'):
            self.parse('x\t1\t5\t1\tW\ta\t1\t5\t+\n'
                       'y\t1\t5\t1\tW\ta\t1\t5\t+\n'
                       'z\t1\t4\t1\tW\tb\t1\t4\t+\n')

    def test_component_bounds_fail(self):
        with self.assertRaisesRegex(ValueError, 'out of bounds'):
            self.parse('x\t1\t6\t1\tW\ta\t1\t6\t+\n')

    def test_object_discontinuity_fails(self):
        with self.assertRaisesRegex(ValueError, 'do not tile'):
            self.parse('x\t2\t6\t1\tW\ta\t1\t5\t+\n')

    def test_unknown_orientation_fails(self):
        with self.assertRaisesRegex(ValueError, 'orientation'):
            self.parse('x\t1\t5\t1\tW\ta\t1\t5\t?\n')


class PairedTests(unittest.TestCase):
    def check_pair(self, first, second):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            a,b = root/'R1.fq',root/'R2.fq'
            a.write_text(first)
            b.write_text(second)
            spec = {'file1':str(a),'file2':str(b),'source_signature':io_utils.signature(a),
                    'source_signature_2':io_utils.signature(b), 'expected_md5':None,
                    'expected_md5_2':None,'sample_id':'S','readset_id':'R','library_id':'L'}
            io_utils.save(root/'spec.json',spec)
            paired_input(SimpleNamespace(spec=root/'spec.json',out=root/'out.json'))
            return io_utils.load(root/'out.json')

    def test_mate_suffixes_and_counts(self):
        result = self.check_pair('@r/1\nACGT\n+\nIIII\n', '@r/2\nTGCA\n+\nIIII\n')
        self.assertEqual(result['num_pairs'],1)

    def test_wrong_pair_ids_fail(self):
        with self.assertRaisesRegex(ValueError,'IDs differ'):
            self.check_pair('@r1\nACGT\n+\nIIII\n','@r2\nACGT\n+\nIIII\n')

    def test_missing_mate_fails(self):
        with self.assertRaisesRegex(ValueError,'unequal'):
            self.check_pair('@r1\nACGT\n+\nIIII\n','')


if __name__ == '__main__':
    unittest.main()
