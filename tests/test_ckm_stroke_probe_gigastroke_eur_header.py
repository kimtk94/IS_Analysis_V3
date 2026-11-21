"""Offline tests for publisher-aligned bounded GIGASTROKE outcome header."""
import gzip
import io
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
from scripts.ckm_stroke.probe_gigastroke_eur_header import header_probe, audit

ACC="GCST90104540"
BASE="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/"
URL=BASE+ACC+"/"+ACC+"_buildGRCh37.tsv.gz"
MD5="1bfe4ae8a5042fb24cb0a562b05b0d2f"
ROW={
    "phenotype":"AIS","accession":ACC,"reported_trait":"Ischemic stroke",
    "build":"GRCh37","ancestry":"EUR","sample_total":"1296908",
    "published_raw_md5":MD5,"metadata_url":URL+"-meta.yaml",
    "verification":"METADATA_VERIFIED_RAW_HASH_UNCHECKED",
}
YAML=f"""gwas_id: {ACC}
genome_assembly: GRCh37
data_file_name: {ACC}_buildGRCh37.tsv.gz
data_file_md5sum: {MD5}
trait_description:
  - Ischemic stroke
sample_size: 1296908
"""
MD5_TEXT=f"{MD5}  {ACC}_buildGRCh37.tsv.gz\n"


class Stub(io.BytesIO):
    def __init__(self,raw,url):
        super().__init__(raw)
        self.status=206
        self._url=url
    def geturl(self):return self._url


def fake_fetch(req,timeout):
    return Stub(gzip.compress(b"chromosome\tbase_pair_location\teffect_allele\tother_allele\tp_value\n"),req.full_url)


def fake_text(uri):
    if uri.endswith("md5sum.txt"):return MD5_TEXT
    if uri.endswith("-meta.yaml"):return YAML
    raise ValueError("unknown URL")


class HeaderTests(unittest.TestCase):
    def test_header_with_sidecars_passes_not_full_file(self):
        result=header_probe(ROW,opener=fake_fetch,text_fetcher=fake_text)
        self.assertEqual(result["header_only_status"],"PASS_PUBLISHER_HEADER_ONLY")
        self.assertFalse(result["full_raw_md5_verified"])
        self.assertFalse(result.get("ready_for_ldsc",False))
        self.assertIn("effect_allele",result["gzip_header_columns"])

    def test_metadata_sidecar_mismatch_aborts(self):
        with self.assertRaisesRegex(ValueError,"trait"):
            header_probe(ROW,opener=fake_fetch,
                         text_fetcher=lambda x:fake_text(x).replace("Ischemic stroke","Wrong stroke"))

    def test_bad_gzip_prefix_rejected(self):
        def not_gzip(req,timeout):return Stub(b"not_gzip_data"*30,req.full_url)
        with self.assertRaisesRegex(ValueError,"gzip"):
            header_probe(ROW,opener=not_gzip,text_fetcher=fake_text)

    def test_header_width_invalid_rejected(self):
        def tiny(req,timeout):return Stub(gzip.compress(b"SNP\tP\n"),req.full_url)
        with self.assertRaisesRegex(ValueError,"width"):
            header_probe(ROW,opener=tiny,text_fetcher=fake_text)

    def test_audit_fail_closed_if_no_manifest(self):
        with TemporaryDirectory() as d:
            with self.assertRaises(FileNotFoundError):
                audit(Path(d)/"missing.tsv")


if __name__=="__main__":
    unittest.main()
