"""Polyploid strategy contracts and complete-copy delivery gates."""
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
import yaml
from tests.assembly import test_config as fixtures

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from assembly.lib.scaffolding import export_polyploid as export
from assembly.lib.io import checksum, load, publish, save
from assembly.lib.assembly import relocate_internal_links
from assembly.lib.assembly import mapping_commands


class CopyNumberInputTests(unittest.TestCase):
    def test_multiple_hifi_readsets_form_one_unpaired_stream(self):
        import gzip
        import subprocess
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first, second = root/'cell1.fq', root/'cell2.fq.gz'
            a, b = '@a\nACGT\n+\nIIII\n', '@b\nTTAA\n+\nIIII\n'
            first.write_text(a)
            with gzip.open(second, 'wt') as handle:
                handle.write(b)
            commands = mapping_commands(root/'unitigs.fa', [first, second], root/'reads.bam', 8)
            streamed = subprocess.run(list(map(str, commands[0])), check=True, capture_output=True).stdout
            self.assertEqual(streamed.decode(), a+b)
            self.assertEqual(commands[1][-1], '-')
            self.assertNotIn(first, commands[1])
            self.assertNotIn(second, commands[1])

    def test_one_hifi_readset_is_mapped_directly(self):
        path = Path('/tmp/one-cell.fq')
        commands = mapping_commands('/tmp/unitigs.fa', [path], '/tmp/reads.bam', 8)
        self.assertEqual(len(commands), 2)
        self.assertEqual(commands[0][-1], path)


class PublicationTests(unittest.TestCase):
    def test_published_graph_links_keep_their_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);stage=root/'new';(stage/'preprocessing').mkdir(parents=True)
            (stage/'preprocessing/counts').write_text('known')
            (stage/'count_link').symlink_to(stage/'preprocessing/counts')
            external=root/'raw_input';external.write_text('input')
            (stage/'source_link').symlink_to(external)
            destination=root/'published'
            relocate_internal_links(stage,destination);publish(stage,destination,root)
            self.assertEqual((destination/'count_link').read_text(),'known')
            self.assertEqual((destination/'source_link').resolve(),external)

    def test_nested_empty_output_parents_are_valid_placeholders(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);stage=root/'new';stage.mkdir();(stage/'result').write_text('valid')
            destination=root/'published';(destination/'nested/output').mkdir(parents=True)
            publish(stage,destination,root)
            self.assertEqual((destination/'result').read_text(),'valid')

    def test_unmanaged_file_is_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);stage=root/'new';stage.mkdir()
            destination=root/'published';(destination/'nested').mkdir(parents=True)
            protected=destination/'nested/user_file';protected.write_text('keep')
            with self.assertRaisesRegex(ValueError,'unmanaged'):publish(stage,destination,root)
            self.assertEqual(protected.read_text(),'keep')


class PolyploidConfigTests(unittest.TestCase):
    def setUp(self):
        self.fixture=fixtures.ConfigTests();self.fixture.setUp();self.addCleanup(self.fixture.doCleanups)
        self.fixture.joint()
        self.fixture.config['defaults'].update(ploidy=4)
        self.fixture.config['defaults']['assembly']['representation']='unitigs'
        self.fixture.config['defaults']['scaffolding']['backend']='gphase'

    def resolve(self):
        f=self.fixture;p=f.path/'config.yaml';p.write_text(yaml.safe_dump(f.config))
        return fixtures.module.resolve_config(p,ROOT,f.path/'run',['A'],16,64000)

    def test_explicit_graph_route_preserves_real_ploidy(self):
        plan=self.resolve()
        self.assertEqual(plan['workflow_route'],'polyploid')
        self.assertEqual(plan['samples']['A']['settings']['ploidy'],4)

    def test_no_diploid_contig_fallback(self):
        self.fixture.config['defaults']['assembly']['representation']='phased'
        with self.assertRaisesRegex(fixtures.module.ConfigError,'polyploid assembly strategy'):self.resolve()

    def test_diploid_gphase_not_silently_enabled(self):
        self.fixture.config['defaults']['ploidy']=2
        with self.assertRaisesRegex(fixtures.module.ConfigError,'GPhase requires'):self.resolve()

    def test_unimplemented_polyploid_review_rejected(self):
        self.fixture.config['defaults']['scaffolding'].update(review_agp='anything',review_manifest='anything')
        with self.assertRaisesRegex(fixtures.module.ConfigError,'polyploid AGP import'):self.resolve()


class PolyploidDeliveryTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.native=self.root/'gphase/gphase_output/scaffold_hap';self.native.mkdir(parents=True)
        raw=''.join(f'>utg{i}\nACGTAC\n' for i in range(1,6))
        (self.root/'raw.fa').write_text(raw)
        save(self.root/'assembly.json',{'sample_id':'P','assembly_id':'known','fasta_sha256':checksum(self.root/'raw.fa')})
        save(self.root/'meta.json',{'sample_id':'P','settings':{'ploidy':4,'chromosomes_per_haplotype':1}})
        self.set_outputs()

    def set_outputs(self, extra=False, wrong=False):
        components=''.join(f'>utg{i}\nACGTAC\n' for i in range(1,5))
        rows=''.join(f'chr1g{i}_scaffold_1\t1\t6\t1\tW\tutg{i}\t1\t6\t+\n' for i in range(1,5))
        fasta=''.join(f'>chr1g{i}_scaffold_1\nACGTAC\n' for i in range(1,5))
        if extra:
            components+='>utg1_dup1\nACGTAC\n'
            rows+='chr1g1_scaffold_2\t1\t6\t1\tW\tutg1_dup1\t1\t6\t+\n'
            fasta+='>chr1g1_scaffold_2\nACGTAC\n'
        for kind in ('unitig','contig'):
            (self.native/f'gphase_final.{kind}.fasta').write_text(components)
            (self.native/f'gphase_final.{kind}.scaffold.agp').write_text(rows)
            (self.native/f'gphase_final.{kind}.scaffold.fasta').write_text(fasta.replace('ACGTAC','ACGTAA',1) if wrong else fasta)
        (self.native/'gphase_final.ctg2utg.txt').write_text('known fixture\n')

    def deliver(self):
        export(SimpleNamespace(metadata=self.root/'meta.json',assembly=self.root/'assembly.json',
               raw_unitigs=self.root/'raw.fa',gphase=self.root/'gphase',out_dir=self.root/'delivery'))
        return load(self.root/'delivery/manifest.json')

    def test_all_four_copies_and_unassigned_sequence_retained(self):
        result=self.deliver()
        self.assertEqual(result['quality_level'],'POLYPLOID_CHROMOSOME_CANDIDATE')
        self.assertEqual(result['sets']['chromosomes']['num_contigs'],4)
        self.assertEqual(result['unassigned_raw_unitigs'],1)
        self.assertIn('unassigned_unitig_utg5',(self.root/'delivery/unplaced.fa').read_text())
        self.assertEqual(result['T2T_status'],'NOT_ASSESSED')

    def test_extra_scaffold_is_not_discarded_to_force_target(self):
        self.set_outputs(extra=True);result=self.deliver()
        self.assertEqual(result['quality_level'],'INCOMPLETE_SCAFFOLDING')
        self.assertEqual(result['sets']['chromosomes']['num_contigs'],5)

    def test_native_fasta_agp_mismatch_blocks_delivery(self):
        self.set_outputs(wrong=True)
        with self.assertRaisesRegex(ValueError,'differs from AGP'):self.deliver()

    def test_raw_unitig_copy_mismatch_blocks_delivery(self):
        p=self.root/'raw.fa';p.write_text(p.read_text().replace('ACGTAC','ACGTAA',1))
        save(self.root/'assembly.json',{'assembly_id':'known','fasta_sha256':checksum(p)})
        with self.assertRaisesRegex(ValueError,'differs from its raw sequence'):self.deliver()


if __name__=='__main__':unittest.main()
