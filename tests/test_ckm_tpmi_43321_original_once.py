"""Offline fail-closed checks for single permitted original TPMI 433.21 upload-free ingest."""
import gzip
import io
import tempfile
import unittest
from pathlib import Path
from scripts.ckm_stroke import inspect_tpmi_43321_original_sumstats_once as m

HEADER="#CHROM\tPOS\tID\tREF\tALT\tA1\tOBS_CT\tBETA\tSE\tP\tEAF"
def rows(n=30):
    return "\n".join([HEADER]+
        ["1\t%d\trs%d\tA\tG\tG\t305000\t0.02\t0.05\t0.68\t0.3"%(1000+i,i)
          for i in range(n)])+"\n"
class Fake:
    status=200
    def __init__(self,b,link=m.URL):
        self.fd=io.BytesIO(b)
        self.headers={"Content-Length":str(len(b))}
        self.link=link
    def __enter__(self):return self
    def __exit__(self,*x):return None
    def geturl(self):return self.link
    def read(self,n=-1):return self.fd.read(n)

class Tests(unittest.TestCase):
    def test_original_gzip_header_and_variant_OBS_CT_not_neff(self):
        with tempfile.TemporaryDirectory() as td:
            file=Path(td)/"x.gz"
            file.write_bytes(gzip.compress(rows().encode()))
            r=m.inspect_header_and_sample(file)
        self.assertEqual(r["technical_status"],
                         "PASS_TPMI_43321_ORIGINAL_GWAS_FULL_BYTE_HASH_AND_BOUNDED_HEADER")
        self.assertTrue(r["source_field_coverage"]["has_variantwise_observation_count_OBS_CT_or_N"])
        self.assertTrue(r["OBS_CT_is_NOT_a_case_control_effective_N"])
        self.assertFalse(r["independence_from_GIGASTROKE_participants_established"])
        self.assertEqual(r["first_bounded_source_records_checked"],30)
    def test_one_download_one_checksum_no_local_persistence(self):
        payload=gzip.compress(rows().encode())
        r=m.read_full_original(opener=lambda req,timeout:Fake(payload),
                               expected_bytes=len(payload))
        self.assertEqual(r["original_file_compressed_bytes"],len(payload))
        self.assertTrue(r["file_deleted_after_source_audit"])
        self.assertTrue(r["file_body_downloaded_once_without_retry"])
        self.assertFalse(r["ready_for_ldsc"])
    def test_wrong_file_size_is_never_accepted(self):
        payload=gzip.compress(rows().encode())
        with self.assertRaisesRegex(ValueError,"content length changed"):
            m.read_full_original(opener=lambda req,timeout:Fake(payload),
                                 expected_bytes=len(payload)+1)
    def test_wrong_original_url_blocked_without_network(self):
        with self.assertRaisesRegex(ValueError,"sole authorized"):
            m.read_full_original(source_url="https://attacker.example/gwas")
    def test_gwas_source_redirect_is_rejected(self):
        payload=gzip.compress(rows().encode())
        with self.assertRaisesRegex(ValueError,"host changed"):
            m.read_full_original(
                opener=lambda req,timeout:Fake(payload,"https://example.org/other"),
                expected_bytes=len(payload))
    def test_bad_header_missing_effect_allele_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            file=Path(td)/"x.gz"
            file.write_bytes(gzip.compress(rows().replace("\tA1\t","\tOTHER\t").replace("\tALT\t","\tALLELE_B\t").encode()))
            with self.assertRaisesRegex(ValueError,"variant CHROM POS P and allele"):
                m.inspect_header_and_sample(file)
    def test_missing_se_is_recorded_not_silently_inferred(self):
        with tempfile.TemporaryDirectory() as td:
            file=Path(td)/"x.gz"
            file.write_bytes(gzip.compress(rows().replace("\tSE\t","\tSE_INVALID\t").encode()))
            r=m.inspect_header_and_sample(file)
            self.assertFalse(r["source_field_coverage"]["has_standard_error"])
            self.assertFalse(r["ready_for_ldsc"])
    def test_bad_p_values_cannot_pass_source_QC(self):
        text=rows().replace("\t0.68\t","\t2.0\t")
        with tempfile.TemporaryDirectory() as td:
            file=Path(td)/"x"
            file.write_text(text)
            with self.assertRaisesRegex(ValueError,"insufficient valid"):
                m.inspect_header_and_sample(file)
    def test_original_rare_compressed_binary_without_support_hold(self):
        with tempfile.TemporaryDirectory() as td:
            file=Path(td)/"x"
            file.write_bytes(b"PK\x03\x04"+"x".encode()*100)
            with self.assertRaisesRegex(ValueError,"unsupported archive"):
                m.inspect_header_and_sample(file)

if __name__=="__main__":
    unittest.main()
