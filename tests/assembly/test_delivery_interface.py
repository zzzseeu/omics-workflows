"""Both routes deliver self-contained AGP sequences and the same public fields."""
import csv
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
import omics_workflows

ROOT = Path(omics_workflows.__file__).resolve().parent
from omics_workflows.workflows.assembly.lib.sequences import parse_agp, reconstruct
from omics_workflows.workflows.assembly.lib.scaffolding import deliver
from omics_workflows.workflows.assembly.lib.scaffolding import export_polyploid as export
from omics_workflows.workflows.assembly.lib.sequences import fasta_dict, write_fasta
from omics_workflows.workflows.assembly.lib.io import checksum, load, save
from omics_workflows.workflows.assembly.lib.postassembly import assembly_bundle


class DeliveryInterfaceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def diploid(self):
        base = self.root / 'diploid'
        native = base / 'haphic/04.build'
        native.mkdir(parents=True)
        tours = base / 'haphic/03.sort/final_tours'
        tours.mkdir(parents=True)
        for group in ('group1', 'group2'):
            (tours / f'{group}_1.tour').write_text('known fixture\n')
        sequences = {'hap1_a': 'AACCGG', 'hap2_b': 'TTGGAA', 'hap1_u': 'CCAA'}
        with (base / 'combined.fa').open('w') as out:
            for name, seq in sequences.items():
                write_fasta(out, name, seq)
        (base / 'mapping.tsv').write_text('contig_id\thaplotype\toriginal_id\n'
                                         'hap1_a\thap1\ta\nhap2_b\thap2\tb\nhap1_u\thap1\tu\n')
        agp = native / 'scaffolds.raw.agp'
        agp.write_text('group1\t1\t3\t1\tW\thap1_a\t1\t3\t+\n'
                       'group1\t4\t6\t2\tU\t3\tscaffold\tyes\tproximity_ligation\n'
                       'group1\t7\t9\t3\tW\thap1_a\t4\t6\t-\n'
                       'group2\t1\t6\t1\tW\thap2_b\t1\t6\t+\n'
                       'unused\t1\t4\t1\tW\thap1_u\t1\t4\t+\n')
        with (native / 'scaffolds.fa').open('w') as out:
            for name, rows in parse_agp(agp, sequences).items():
                write_fasta(out, name, reconstruct(rows, sequences))
        save(base / 'meta.json', {'sample_id': 'S', 'settings': {'chromosomes_per_haplotype': 1,
             'scaffolding': {'review_agp': None, 'review_manifest': None}}})
        save(base / 'assembly.json', {'assembly_id': 'diploid', 'phase_level': 'HIC_PHASED_UNVALIDATED',
             'combined_sha256': checksum(base / 'combined.fa')})
        save(base / 'haphic/manifest.json', {'hic_reference': 'combined'})
        deliver(SimpleNamespace(metadata=base / 'meta.json', assembly=base / 'assembly.json',
                fasta=base / 'combined.fa', mapping=base / 'mapping.tsv',
                haphic=base / 'haphic', out_dir=base / 'delivery'))
        return base / 'delivery'

    def polyploid(self, unplaced=True):
        base = self.root / 'polyploid'
        native = base / 'gphase/gphase_output/scaffold_hap'
        native.mkdir(parents=True)
        raw = '>utg1\nAACCGG\n' + ('>unused\nTTAA\n' if unplaced else '')
        (base / 'raw.fa').write_text(raw)
        # One collapsed raw unitig is assigned to all four copies. The graph
        # components use only part of it after overlap trimming, unlike diploids.
        unitigs = ''.join(f'>utg1{("_dup"+str(i)) if i else ""}\nAACCGG\n' for i in range(4))
        (native / 'gphase_final.unitig.fasta').write_text(unitigs)
        (native / 'gphase_final.contig.fasta').write_text(''.join(f'>ctg{i}\nAACC\n' for i in range(4)))
        for kind in ('unitig', 'contig'):
            agp = native / f'gphase_final.{kind}.scaffold.agp'
            rows = []
            for i in range(4):
                component = 'utg1' + (f'_dup{i}' if i else '') if kind == 'unitig' else f'ctg{i}'
                rows.append(f'chr1g{i+1}_scaffold_1\t1\t4\t1\tW\t{component}\t1\t4\t+\n')
            agp.write_text(''.join(rows))
            (native / f'gphase_final.{kind}.scaffold.fasta').write_text(
                ''.join(f'>chr1g{i+1}_scaffold_1\nAACC\n' for i in range(4)))
        (native / 'gphase_final.ctg2utg.txt').write_text('known graph fixture\n')
        save(base / 'meta.json', {'sample_id': 'S', 'settings': {'ploidy': 4, 'chromosomes_per_haplotype': 1}})
        save(base / 'assembly.json', {'assembly_id': 'polyploid', 'fasta_sha256': checksum(base / 'raw.fa')})
        export(SimpleNamespace(metadata=base / 'meta.json', assembly=base / 'assembly.json',
               raw_unitigs=base / 'raw.fa', gphase=base / 'gphase', out_dir=base / 'delivery'))
        return base / 'delivery'

    def test_common_files_reconstruct_both_routes_without_original_directories(self):
        for base, ploidy in ((self.diploid(), 2), (self.polyploid(), 4)):
            with self.subTest(ploidy=ploidy):
                manifest, chromosomes, groups = assembly_bundle(base / 'manifest.json')
                for name in manifest['artifacts']['sequence_sets'].values():
                    self.assertTrue((base / name).is_file(), name)
                self.assertEqual(manifest['analysis_stage'], 'phasing_scaffolding')
                self.assertEqual(manifest['ploidy'], ploidy)
                self.assertEqual(manifest['phase_validation'], 'NOT_ASSESSED')
                self.assertEqual(len(set(groups.values())), ploidy)
                components = fasta_dict(base / 'components.fa')
                objects = parse_agp(base / 'assembly.agp', components, require_complete=False)
                chr_objects = parse_agp(base / 'chromosomes.agp', components, require_complete=False)
                self.assertEqual(set(chr_objects), set(chromosomes))
                delivered = {**chromosomes, **fasta_dict(base / 'unplaced.fa')}
                for name, rows in objects.items():
                    self.assertEqual(reconstruct(rows, components), delivered[name])
                with (base / 'chromosome_names.tsv').open() as handle:
                    reader = csv.DictReader(handle, delimiter='\t')
                    self.assertEqual(reader.fieldnames,
                                     ['sequence_id', 'source_object', 'haplotype', 'status', 'length'])
                    table = {row['sequence_id']: row for row in reader}
                self.assertEqual(set(table), set(delivered))
                for name, seq in delivered.items():
                    self.assertEqual(int(table[name]['length']), len(seq))
                for label, fasta in manifest['artifacts']['sequence_sets'].items():
                    expected = {name: len(seq) for name, seq in fasta_dict(base / fasta).items()}
                    sizes = {name: int(size) for name, size in
                             (line.split('\t') for line in (base / f'{label}.chrom.sizes').read_text().splitlines())}
                    self.assertEqual(sizes, expected)
                for line in (base / 'SHA256SUMS').read_text().splitlines():
                    digest, name = line.split('  ', 1)
                    self.assertEqual(checksum(base / name), digest)

    def test_empty_unplaced_is_a_valid_sequence_set(self):
        base = self.polyploid(unplaced=False)
        for name in ('unplaced.fa', 'unplaced.fa.fai', 'unplaced.chrom.sizes'):
            self.assertEqual((base / name).read_text(), '')
        manifest, _, _ = assembly_bundle(base / 'manifest.json')
        self.assertEqual(manifest['sets']['unplaced']['num_contigs'], 0)


if __name__ == '__main__':
    unittest.main()
