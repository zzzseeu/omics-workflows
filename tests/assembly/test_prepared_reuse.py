"""Prepared-input identity, immutable caching and zero-copy reuse contracts."""
import argparse
import copy
import gzip
import json
from pathlib import Path
import tempfile
import unittest
import omics_workflows
from unittest.mock import patch

ROOT=Path(omics_workflows.__file__).resolve().parent
from omics_workflows.assembly.lib.io import checksum, save, signature
from omics_workflows.assembly.lib.reads import describe_prepared_input
from omics_workflows.assembly.lib.reads import attest, reuse


class PreparedReuseTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        parent=self.root/'old/samples/A/01_inputs';parent.mkdir(parents=True)
        self.raw=self.root/'raw.bam';self.raw.write_bytes(b'raw source provenance')
        self.fastq=parent/'r1.fastq.gz'
        with gzip.open(self.fastq,'wt') as out:out.write('@read\nACGTACGT\n+\nIIIIIIII\n')
        self.manifest=parent/'r1.json';self.qc=parent.parent/'02_read_qc/r1.json'
        self.spec={'sample_id':'A','readset_id':'r1','library_id':'L','technology':'hifi','format':'bam',
                   'file1':str(self.raw),'source_signature':signature(self.raw),'expected_md5':checksum(self.raw,'md5'),
                   'seqkit':{'path':'/must-not-be-called/seqkit'}}
        self.old={'sample_id':'A','readset_id':'r1','library_id':'L','source':str(self.raw),
                  'source_signature':signature(self.raw),'source_md5':checksum(self.raw,'md5'),
                  'source_sha256':checksum(self.raw),'provided_md5_status':'PASS',
                  'conversion':'samtools fastq; unpaired primary reads; original IDs retained',
                  'filtering':'No quality trimming or deduplication performed',
                  'prepared_fastq':'samples/A/01_inputs/r1.fastq.gz'}
        save(self.manifest,self.old)
        save(self.qc,{'status':'PASS','num_reads':1,'num_bases':8,
            'statistics':{'file':'samples/A/01_inputs/r1.fastq.gz','format':'FASTQ','num_seqs':'1','sum_len':'8'}})

    def descriptor(self):
        return describe_prepared_input(self.spec,self.manifest,self.root/'cache',self.root/'new',ROOT)

    def producer(self):
        save(self.manifest,{**self.old,'prepared_sha256':checksum(self.fastq),'prepared_signature':signature(self.fastq)})
        return {**self.spec,'prepared_reuse':self.descriptor()}

    def test_legacy_manifest_is_metadata_only_and_requires_initial_attestation(self):
        r=self.descriptor();self.assertIsNone(r['producer_sha256']);self.assertFalse(Path(r['cache_file']).exists())

    def test_source_sample_cannot_be_substituted(self):
        self.old['sample_id']='B';save(self.manifest,self.old)
        with self.assertRaisesRegex(ValueError,'identity mismatch'):self.descriptor()

    def test_filtering_policy_cannot_be_substituted(self):
        self.old['filtering']='filtered again';save(self.manifest,self.old)
        with self.assertRaisesRegex(ValueError,'filtering policy'):self.descriptor()

    def test_changed_raw_source_is_rejected(self):
        self.raw.write_bytes(b'changed source');self.spec['source_signature']=signature(self.raw)
        with self.assertRaisesRegex(ValueError,'source changed'):self.descriptor()

    def test_nonpass_qc_is_rejected(self):
        q=json.loads(self.qc.read_text());q['status']='FAILED';save(self.qc,q)
        with self.assertRaisesRegex(ValueError,'QC is not complete'):self.descriptor()

    def test_wrong_qc_file_is_rejected(self):
        q=json.loads(self.qc.read_text());q['statistics']['file']='samples/B/01_inputs/r1.fastq.gz';save(self.qc,q)
        with self.assertRaisesRegex(ValueError,'different FASTQ'):self.descriptor()

    def test_producer_hash_uses_metadata_only_cache(self):
        spec=self.producer()
        with patch('omics_workflows.assembly.lib.reads.scan_legacy',side_effect=AssertionError('must not scan')):
            first,hit,sha=attest(spec);self.assertFalse(hit)
            second,hit,sha=attest(spec);self.assertTrue(hit)
        self.assertEqual(first['prepared_sha256'],checksum(self.fastq));self.assertEqual(first,second)

    def test_changed_prepared_file_is_not_silently_recertified(self):
        spec=self.producer();attest(spec)
        self.fastq.touch();changed={**self.spec,'prepared_reuse':self.descriptor()}
        with patch('omics_workflows.assembly.lib.reads.scan_legacy',side_effect=AssertionError('must not rescan')):
            with self.assertRaisesRegex(ValueError,'Cached prepared FASTQ changed'):attest(changed)

    def test_link_does_not_copy_or_touch_original(self):
        spec=self.producer();before=signature(self.fastq);path=self.root/'spec.json';save(path,spec)
        out=self.root/'new/samples/A/01_inputs/r1.fastq.gz'
        reuse(argparse.Namespace(spec=str(path),fastq=str(out),manifest=str(out.with_name('r1.json'))))
        self.assertTrue(out.is_symlink());self.assertEqual(out.resolve(),self.fastq)
        self.assertEqual(signature(self.fastq),before)
        m=json.loads(out.with_name('r1.json').read_text())
        self.assertEqual(m['source_sha256'],checksum(self.raw));self.assertEqual(m['prepared_sha256'],checksum(self.fastq))

    def test_cache_provenance_tampering_is_rejected(self):
        spec=self.producer();attest(spec);p=Path(spec['prepared_reuse']['cache_file'])
        data=json.loads(p.read_text());data['identity']['sample_id']='wrong';save(p,data)
        with self.assertRaisesRegex(ValueError,'provenance is invalid'):attest(spec)

    def test_stage_one_validates_reused_reads_without_seqkit(self):
        import importlib
        from types import SimpleNamespace
        prepare = importlib.import_module("omics_workflows.assembly.scripts.01_prepare_inputs").prepare_sample
        spec=self.producer()
        metadata=self.root/'metadata.json'
        save(metadata, {'sample_id':'A','settings':{'assembly':{'strategy':'single'}},'readsets':{'r1':spec}})
        output=self.root/'new/samples/A/01_inputs'
        prepare(SimpleNamespace(metadata=metadata,out_dir=output,threads=1))
        record=json.loads((output/'manifest.json').read_text())['readsets']['r1']
        self.assertEqual(record['statistics']['num_reads'],1)
        self.assertEqual(record['statistics']['num_bases'],8)
        self.assertIn('reuse',record['sources'])
        self.assertTrue((output/record['files'][0]['file']).is_symlink())


if __name__=='__main__':unittest.main()
