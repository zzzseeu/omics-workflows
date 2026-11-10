"""BAM dictionary translation must preserve genomic and mate coordinates."""
from pathlib import Path
import sys
import tempfile
import unittest
import pysam

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from assembly.lib.hic import bam_pairs, combine_haplotype_bams


class HaplotypeBamTests(unittest.TestCase):
    def write_pair(self, path, hap):
        header = {'HD': {'VN': '1.6', 'SO': 'unsorted'},
                  'SQ': [{'SN': f'{hap}_a', 'LN': 1000}, {'SN': f'{hap}_b', 'LN': 2000}]}
        with pysam.AlignmentFile(str(path), 'wb', header=header) as bam:
            for mate in (0, 1):
                read = pysam.AlignedSegment(bam.header)
                read.query_name = 'same_molecule'
                read.flag = 65 if mate == 0 else 129
                read.reference_id = mate
                read.reference_start = 100 + mate
                read.next_reference_id = 1-mate
                read.next_reference_start = 101-mate
                read.mapping_quality = 60
                read.cigarstring = '10M'
                read.query_sequence = 'ACGTACGTAC'
                read.query_qualities = pysam.qualitystring_to_array('I'*10)
                read.set_tag('NM', 0)
                bam.write(read)

    def test_translate_reference_and_mate_dictionaries(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = [Path(tmp)/f'{i}.bam' for i in (1,2)]
            for i,path in enumerate(paths,1): self.write_pair(path, f'hap{i}')
            out = Path(tmp)/'union.bam'
            combine_haplotype_bams(paths, out)
            self.assertEqual(bam_pairs(out), 2)
            with pysam.AlignmentFile(str(out)) as bam:
                reads = list(bam)
            self.assertEqual([r.reference_name for r in reads], ['hap1_a','hap1_b','hap2_a','hap2_b'])
            self.assertEqual([r.next_reference_name for r in reads], ['hap1_b','hap1_a','hap2_b','hap2_a'])
            self.assertEqual([r.reference_start for r in reads], [100,101,100,101])
            self.assertEqual([r.next_reference_start for r in reads], [101,100,101,100])
            self.assertEqual([r.query_name for r in reads], ['hap1:same_molecule']*2+['hap2:same_molecule']*2)
            self.assertTrue(all(r.get_tag('NM') == 0 for r in reads))

    def test_overlapping_reference_names_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = [Path(tmp)/f'{i}.bam' for i in (1,2)]
            for path in paths: self.write_pair(path, 'hap1')
            with self.assertRaisesRegex(ValueError, 'disjoint'):
                combine_haplotype_bams(paths, Path(tmp)/'bad.bam')


if __name__ == '__main__': unittest.main()
