"""Biologically material delivery gates: exact coverage, haplotypes and stale reviews."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from assembly.lib import sequences as agp_utils
from assembly.lib import io as io_utils
from assembly.lib.postassembly import fastq_statistics


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
            return fastq_statistics([a,b])

    def test_mate_suffixes_and_counts(self):
        result = self.check_pair('@r/1\nACGT\n+\nIIII\n', '@r/2\nTGCA\n+\nIIII\n')
        self.assertEqual(result['num_pairs'],1)

    def test_wrong_pair_ids_fail(self):
        with self.assertRaisesRegex(ValueError,'IDs/order differ'):
            self.check_pair('@r1\nACGT\n+\nIIII\n','@r2\nACGT\n+\nIIII\n')

    def test_missing_mate_fails(self):
        with self.assertRaisesRegex(ValueError,'counts differ'):
            self.check_pair('@r1\nACGT\n+\nIIII\n','')


if __name__ == '__main__':
    unittest.main()
