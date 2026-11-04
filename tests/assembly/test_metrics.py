"""Small, independently calculable statistics and empty/invalid set contracts."""
from pathlib import Path
import tempfile
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'assembly/scripts'))
import fasta_utils as tasks


class MetricTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.fasta = Path(self.temp.name) / 'assembly.fa'

    def test_known_lengths_and_bases(self):
        self.fasta.write_text('>long\nAC\nGN\n>short\nAA\n>tiny\nT\n')
        stat = tasks.fasta_statistics(self.fasta)
        self.assertEqual([stat[k] for k in ['num_contigs', 'total_bases', 'N50', 'L50', 'N_bases']],
                         [3, 7, 4, 1, 1])
        self.assertAlmostEqual(stat['GC_percent'], 200 / 7)

    def test_empty_alternate_set(self):
        self.fasta.write_text('')
        stat = tasks.fasta_statistics(self.fasta)
        self.assertEqual(stat['num_contigs'], 0)
        self.assertEqual(stat['N50'], 0)
        self.assertIsNone(stat['GC_percent'])

    def test_duplicate_contig_ids(self):
        self.fasta.write_text('>same\nACGT\n>same\nACGT\n')
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            tasks.fasta_statistics(self.fasta)

    def test_sequence_without_header(self):
        self.fasta.write_text('ACGT\n')
        with self.assertRaisesRegex(ValueError, 'without a header'):
            tasks.fasta_statistics(self.fasta)


if __name__ == '__main__':
    unittest.main()
