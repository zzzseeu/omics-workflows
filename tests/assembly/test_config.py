"""Contract failures that must be caught before starting compute rules."""
import copy
from pathlib import Path
import tempfile
import unittest
import omics_workflows
import yaml

ROOT = Path(omics_workflows.__file__).resolve().parent
from omics_workflows.assembly import runner as module


class ConfigTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)
        (self.path / 'reads.fq').write_text('@r\nACGT\n+\nIIII\n')
        (self.path / 'samples.tsv').write_text('sample_id\tspecies\tprofile\nA\tSpecies\t\nB\tOther\t\n')
        (self.path / 'readsets.tsv').write_text(
            'sample_id\treadset_id\tlibrary_id\ttechnology\tformat\tfile1\n'
            'A\tr1\tl1\thifi\tfastq\treads.fq\n'
            'B\tr1\tl1\thifi\tfastq\tmissing.fq\n')
        self.config = {'samples_table': 'samples.tsv',
                       'readsets_table': 'readsets.tsv',
                       'defaults': {'assembly': {'mem_mb': 2048, 'bloom_filter_bits': 0}}}

    def resolve(self, config=None, selected=None, budget=4096):
        (self.path / 'config.yaml').write_text(yaml.safe_dump(config or self.config))
        return module.resolve_config(self.path / 'config.yaml', ROOT, self.path / 'run',
                                     selected if selected is not None else ['A'], 2, budget)

    def test_selection_does_not_require_unselected_files(self):
        result = self.resolve()
        self.assertEqual(result['sample_order'], ['A'])
        self.assertIsNone(result['samples']['A']['settings']['ploidy'])

    def test_plan_contains_only_current_run_fields(self):
        self.assertEqual(set(self.resolve()), {
            'layout_version', 'config_path', 'repository', 'framework', 'run_dir',
            'budget', 'sample_order', 'samples', 'workflow_route',
        })

    def test_genome_profiling_is_not_an_assembly_stage(self):
        self.assertNotIn('genome_profile', self.resolve()['samples']['A']['settings'])

    def test_old_disabled_profiling_placeholder_preserves_settings(self):
        baseline = self.resolve()['samples']['A']['settings']
        for scope in ('defaults', 'analysis_profiles', 'samples'):
            with self.subTest(scope=scope):
                config = copy.deepcopy(self.config)
                if scope == 'defaults':
                    config['defaults']['genome_profile'] = {'enabled': False}
                elif scope == 'samples':
                    config['samples'] = {'A': {'genome_profile': {'enabled': False}}}
                else:
                    config['analysis_profiles'] = {'unused': {'genome_profile': {'enabled': False}}}
                self.assertEqual(self.resolve(config)['samples']['A']['settings'], baseline)
                self.assertIn('genome_profile', config[scope] if scope == 'defaults' else
                              next(iter(config[scope].values())))

    def test_enabled_or_malformed_profiling_is_not_silently_ignored(self):
        for value in ({'enabled': True}, {'enabled': False, 'enabeld': False}, {'enabled': 0}, None):
            with self.subTest(value=value):
                config = copy.deepcopy(self.config)
                config['defaults']['genome_profile'] = value
                with self.assertRaisesRegex(module.ConfigError, 'outside the assembly workflow'):
                    self.resolve(config)

    def test_missing_selected_input(self):
        with self.assertRaisesRegex(module.ConfigError, 'missing or empty'):
            self.resolve(selected=['B'])

    def test_duplicate_selected_sample(self):
        with self.assertRaisesRegex(module.ConfigError, 'Duplicate or unknown'):
            self.resolve(selected=['A', 'A'])

    def test_known_polyploid_requires_its_own_strategy(self):
        self.config['samples'] = {'A': {'ploidy': 4}}
        with self.assertRaisesRegex(module.ConfigError, 'polyploid assembly strategy'):
            self.resolve()

    def test_unknown_mode_does_not_fallback(self):
        self.config['defaults']['assembly']['mode'] = 'hifi_ul_hic'
        with self.assertRaisesRegex(module.ConfigError, 'no automatic fallback'):
            self.resolve()

    def test_optional_module_not_silently_ignored(self):
        self.config['defaults']['finishing'] = {'enabled': True}
        with self.assertRaisesRegex(module.ConfigError, 'not implemented'):
            self.resolve()

    def test_unavailable_evaluation_not_silently_passed(self):
        self.config['defaults']['assessment'] = {'busco': True}
        with self.assertRaisesRegex(module.ConfigError, 'basic contig statistics only'):
            self.resolve()

    def test_budget_failure(self):
        with self.assertRaisesRegex(module.ConfigError, 'task budget'):
            self.resolve(budget=1024)

    def test_bloom_memory_failure(self):
        self.config['defaults']['assembly']['bloom_filter_bits'] = 37
        with self.assertRaisesRegex(module.ConfigError, 'bloom_filter_bits'):
            self.resolve()

    def test_typo_is_not_ignored(self):
        self.config['defaults']['assembly']['thread'] = 2
        with self.assertRaises(module.ConfigError):
            self.resolve()

    def test_override_precedence(self):
        self.config['analysis_profiles'] = {'p': {'assembly': {'threads': 1}}}
        (self.path / 'samples.tsv').write_text('sample_id\tspecies\tprofile\nA\tSpecies\tp\nB\tOther\t\n')
        self.config['samples'] = {'A': {'assembly': {'threads': 2}}}
        result = self.resolve()
        self.assertEqual(result['samples']['A']['settings']['assembly']['threads'], 2)

    def test_duplicate_readset_rejected(self):
        with (self.path / 'readsets.tsv').open('a') as f:
            f.write('A\tr1\tl1\thifi\tfastq\treads.fq\n')
        with self.assertRaisesRegex(module.ConfigError, 'Duplicate readset_id'):
            self.resolve()

    def test_repository_output_rejected(self):
        (self.path / 'config.yaml').write_text(yaml.safe_dump(self.config))
        with self.assertRaisesRegex(module.ConfigError, 'outside'):
            module.resolve_config(self.path / 'config.yaml', ROOT, ROOT / 'results', ['A'], 2, 4096)

    def test_legacy_run_is_not_silently_reassembled(self):
        import json
        run = self.path / 'run'
        (run / 'config').mkdir(parents=True)
        (run / 'samples').mkdir()
        (run / 'config/plan.json').write_text(json.dumps({}))
        with self.assertRaisesRegex(module.ConfigError, 'Legacy output layout'):
            self.resolve()

    def test_migrated_completed_run_is_historical(self):
        run = self.path / 'run'
        (run / 'layout_migration').mkdir(parents=True)
        (run / 'layout_migration/status.json').write_text('{}')
        with self.assertRaisesRegex(module.ConfigError, 'historical result'):
            self.resolve()

    def joint(self, mode='hifi_ul_hic'):
        self.config['defaults'].update(ploidy=2, chromosomes_per_haplotype=12,
                                      scaffolding={'enabled': True, 'mem_mb': 32000})
        self.config['defaults']['assembly'].update(mode=mode, representation='phased')
        for name in ('ont.fq', 'R1.fq', 'R2.fq'):
            (self.path/name).write_text('@r\nACGT\n+\nIIII\n')
        (self.path/'readsets.tsv').write_text(
            'sample_id\treadset_id\tlibrary_id\ttechnology\tformat\tfile1\tfile2\n'
            'A\th\tH\thifi\tfastq\treads.fq\t\n'
            'A\tu\tU\tont\tfastq\tont.fq\t\n'
            'A\tc\tC\thic\tfastq\tR1.fq\tR2.fq\n')

    def resolve_joint(self):
        (self.path/'config.yaml').write_text(yaml.safe_dump(self.config))
        return module.resolve_config(self.path/'config.yaml',ROOT,self.path/'run',['A'],8,64000)

    def test_joint_activates_all_required_technologies(self):
        self.joint()
        plan=self.resolve_joint()
        self.assertEqual({r['technology'] for r in plan['samples']['A']['readsets'].values()},
                         {'hifi','ont','hic'})
        self.assertEqual(plan['samples']['A']['settings']['scaffolding']['processes'],4)

    def test_independent_haplotype_mapping_is_explicit(self):
        self.joint()
        self.config['defaults']['scaffolding']['hic_reference']='per_haplotype'
        self.assertEqual(self.resolve_joint()['samples']['A']['settings']['scaffolding']['hic_reference'],
                         'per_haplotype')

    def test_independent_haplotype_mapping_rejects_contig_only_route(self):
        self.config['defaults']['scaffolding']={'hic_reference':'per_haplotype'}
        with self.assertRaisesRegex(module.ConfigError,'diploid phased HapHiC'):
            self.resolve()

    def test_hifi_hic_does_not_activate_ont(self):
        self.joint('hifi_hic')
        (self.path/'ont.fq').unlink()
        self.assertEqual(self.resolve_joint()['samples']['A']['unused_readsets'],['u'])

    def test_joint_requires_explicit_karyotype(self):
        self.joint()
        self.config['defaults']['ploidy']=None
        with self.assertRaisesRegex(module.ConfigError,'explicit ploidy=2'):
            self.resolve_joint()

    def test_joint_missing_technology_never_falls_back(self):
        self.joint()
        path=self.path/'readsets.tsv'
        path.write_text(''.join(line for line in path.read_text().splitlines(True) if '\tont\t' not in line))
        with self.assertRaisesRegex(module.ConfigError,'required read technologies missing'):
            self.resolve_joint()

    def test_joint_requires_both_mates(self):
        self.joint()
        (self.path/'R2.fq').unlink()
        with self.assertRaisesRegex(module.ConfigError,'missing or empty external mate'):
            self.resolve_joint()

    def test_finishing_requires_explicit_motif(self):
        self.joint();self.config['defaults']['finishing']={'enabled':True}
        with self.assertRaisesRegex(module.ConfigError,'telomere_motif'):self.resolve_joint()

    def test_nextpolish_does_not_accept_hic_as_wgs(self):
        self.joint();self.config['defaults']['finishing']={'enabled':True,'telomere_motif':'TTTAGGG'}
        self.config['defaults']['polishing']={'enabled':True,'on_missing_wgs':'error'}
        with self.assertRaisesRegex(module.ConfigError,'ordinary WGS'):self.resolve_joint()

    def test_missing_wgs_can_be_explicitly_skipped(self):
        self.joint();self.config['defaults']['finishing']={'enabled':True,'telomere_motif':'TTTAGGG'}
        self.config['defaults']['polishing']={'enabled':True,'on_missing_wgs':'skip'}
        settings=self.resolve_joint()['samples']['A']['settings']
        self.assertEqual(settings['polishing']['on_missing_wgs'],'skip')

    def test_busco_requires_existing_offline_lineage(self):
        self.joint();self.config['defaults']['finishing']={'enabled':True,'telomere_motif':'TTTAGGG'}
        self.config['defaults']['assessment']={'busco':True}
        with self.assertRaisesRegex(module.ConfigError,'offline lineage'):self.resolve_joint()

    def test_review_requires_both_files(self):
        self.joint()
        self.config['defaults']['scaffolding']['review_agp']='review.agp'
        with self.assertRaisesRegex(module.ConfigError,'supplied together'):
            self.resolve_joint()


if __name__ == '__main__':
    unittest.main()
