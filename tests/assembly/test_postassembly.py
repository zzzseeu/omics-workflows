"""Post-assembly contracts that protect copy identity and evidence integrity."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from assembly.lib.sequences import fasta_statistics
from assembly.lib.io import checksum
from assembly.lib.postassembly import inspect_sequences, rebase_symlinks, supported_alignments, terminal_repeats, assembly_bundle, write_sequences


class PostassemblyTests(unittest.TestCase):
    def test_internal_telomeric_repeat_is_not_terminal(self):
        seq='ACGT'*2000+'TTTAGGG'*40+'ACGT'*2000
        self.assertFalse(any(terminal_repeats(seq,'TTTAGGG',10000,100,10).values()))

    def test_both_terminal_orientations_and_rotations(self):
        seq='TAGGGTT'*20+'ACGT'*2000+'CCCTAAA'*20
        self.assertTrue(all(terminal_repeats(seq,'TTTAGGG',10000,10,10).values()))

    def test_ambiguity_blocks_gap_free_screening(self):
        rows=inspect_sequences({'chr':'TTTAGGG'*20+'ACGTR'+'CCCTAAA'*20},
            {'terminal_window':10000,'terminal_max_distance':10,'terminal_min_repeats':10},'TTTAGGG')
        self.assertEqual(rows[0]['N_bases'],0)
        self.assertEqual(rows[0]['ambiguous_bases'],1)
        self.assertFalse(rows[0]['gap_free_two_telomeres'])

    def test_changed_copy_fasta_is_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)
            for name,seq in [('chromosomes','ACGT'),('hap1.chromosomes','ACGA'),('unplaced','')]:
                (p/(name+'.fa')).write_text('>chr\n'+seq+'\n' if seq else '')
            manifest={'sample_id':'S','quality_level':'CHROMOSOME_CANDIDATE',
                      'sets':{label:fasta_statistics(p/(label+'.fa')) for label in ['chromosomes','hap1.chromosomes','unplaced']}}
            (p/'manifest.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,'inconsistent sequence sets'):assembly_bundle(p/'manifest.json')

    def test_stale_source_hash_is_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);(p/'chromosomes.fa').write_text('>chr\nACGT\n')
            (p/'manifest.json').write_text(json.dumps({'sample_id':'S','sets':{'chromosomes':{'sha256':'bad'}}}))
            with self.assertRaisesRegex(ValueError,'changed after delivery'):assembly_bundle(p/'manifest.json')

    def test_native_evidence_links_survive_publication(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);stage=root/'stage';stage.mkdir();(stage/'data').write_text('evidence')
            (stage/'link').symlink_to(stage/'data');dest=root/'final'
            rebase_symlinks(stage,dest);stage.replace(dest)
            self.assertEqual((dest/'link').read_text(),'evidence')

    def test_ambiguous_and_secondary_reads_cannot_support_copy_specific_edits(self):
        import pysam
        records=[]
        for name,mapq,flag in [('ambiguous',0,0),('secondary',60,256),('unique',60,0)]:
            r=pysam.AlignedSegment();r.query_name=name;r.query_sequence='A'*1000
            r.flag=flag;r.reference_start=0;r.mapping_quality=mapq;r.cigarstring='1000M';r.set_tag('NM',0)
            records.append(r)
        class Bam:
            def fetch(self,*args):return iter(records)
        self.assertEqual([r.query_name for r in supported_alignments(Bam(),'chr',0,1000,30,.99)],['unique'])

    def test_duplicate_primary_molecule_id_is_rejected(self):
        import pysam
        r=pysam.AlignedSegment();r.query_name='same_molecule';r.query_sequence='A'*1000
        r.flag=0;r.reference_start=0;r.mapping_quality=60;r.cigarstring='1000M';r.set_tag('NM',0)
        class Bam:
            def fetch(self,*args):return iter([r,r])
        with self.assertRaisesRegex(ValueError,'Duplicate primary'):
            list(supported_alignments(Bam(),'chr',0,1000,30,.99))


if __name__=='__main__':unittest.main()
