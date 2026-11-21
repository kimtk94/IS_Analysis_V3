"""No-network TPMI 433.21 full source GWAS integrity and scientific HOLD tests."""
import gzip
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from scripts.ckm_stroke import audit_tpmi_43321_full_gwas_qc as m

HEADER="\t".join(m.EXPECTED_HEADER)
def variant(i, *, model="ADD", a1="G", ref="A", alt="G", af="0.30",
            obs="305000", odds="1.02", se="0.05", z="0.396", p="0.69", err="NA"):
    return "\t".join([
      "1",str(10000+i),f"1:{10000+i}:{ref}:{alt}",ref,alt,"N",a1,
      "NA",af,model,obs,odds,se,z,p,err,f"rs{i}"])

def payload(lines=None):
    if lines is None:
        lines=[variant(i) for i in range(60)]
    return gzip.compress((HEADER+"\n"+"\n".join(lines)+"\n").encode("utf-8"))

def write_blob(data):
    class Temp:
        def __enter__(self):
            self.d=tempfile.TemporaryDirectory()
            self.path=Path(self.d.name)/"gwas.gz"
            self.path.write_bytes(data)
            return self.path
        def __exit__(self,*args):
            self.d.cleanup()
    return Temp()

class RemoteResponse:
    status=200
    def __init__(self,data,url=m.URL):
        self.buf=io.BytesIO(data)
        self.url=url
        self.headers={"Content-Length":str(len(data))}
    def __enter__(self):return self
    def __exit__(self,*args):return None
    def read(self,n=-1):return self.buf.read(n)
    def geturl(self):return self.url

class AllSourceTests(unittest.TestCase):
    def test_pinned_header_all_rows_and_neff_hold(self):
        with write_blob(payload()) as path:
            r=m.scan_all_rows(path)
        self.assertEqual(r["rows_total_in_original"],60)
        self.assertEqual(r["row_level_QC_counts"]["fully_numeric_ADD_no_error"],60)
        self.assertEqual(r["original_additive_test_categories"]["ADD"],60)
        self.assertEqual(r["value_min_max"]["min_OBS_CT"],305000)
        self.assertTrue(r["A1_is_tested_allele_not_assume_ALT"])
        self.assertFalse(r["GIGASTROKE_and_TPMI_individual_overlap_resolved"])
        self.assertFalse(r["ready_for_ldsc"])

    def test_same_variant_key_detected_but_not_exposed(self):
        lines=[variant(i) for i in range(20)]
        lines.append(lines[0])
        with write_blob(payload(lines)) as f:
            out=m.scan_all_rows(f)
        self.assertEqual(out["duplicate_variant_test_128bit_hashes"],1)
        self.assertNotIn("rs0",json.dumps(out))
        self.assertTrue(out["duplicate_hash_is_not_mathematical_proof_of_duplicate_identity"])

    def test_p_value_underflow_frequency_and_missing_statistic_inventoried(self):
        lines=[variant(i) for i in range(40)]
        lines[5]=variant(5,p="0")
        lines[6]=variant(6,p="NA",odds="NA",err="UNFINISHED")
        lines[7]=variant(7,af="0.00003")
        lines[8]=variant(8,af="1.1")
        lines[9]=variant(9,odds="0")
        lines[10]=variant(10,obs="NA")
        lines[11]=variant(11,z="40")
        with write_blob(payload(lines)) as f:
            out=m.scan_all_rows(f)
        c=out["row_level_QC_counts"]
        self.assertEqual(c["P_reported_zero"],1)
        self.assertEqual(c["P_invalid_or_missing"],1)
        self.assertEqual(c["invalid_or_missing_A1_FREQ"],1)
        self.assertEqual(c["OR_invalid_or_missing"],2)
        self.assertEqual(c["OBS_CT_nonpositive_or_missing"],1)
        self.assertTrue(out["OBS_CT_is_not_SNP_case_control_effective_N"])

    def test_corrupt_column_count_hard_fails(self):
        lines=[variant(i) for i in range(30)]
        lines[12]="1\tbroken\talleles"
        with write_blob(payload(lines)) as f:
            with self.assertRaisesRegex(ValueError,"TSV truncated"):
                m.scan_all_rows(f)

    def test_forged_effect_allele_is_flagged_not_approved(self):
        lines=[variant(i) for i in range(30)]
        lines[0]=variant(0,a1="C")
        with write_blob(payload(lines)) as f:
            r=m.scan_all_rows(f)
        self.assertEqual(r["row_level_QC_counts"]["A1_not_REF_ALT_or_invalid_DNA"],1)
        self.assertEqual(r["row_level_QC_counts"]["fully_numeric_ADD_no_error"],29)
        self.assertFalse(r["all_SNPs_eligible_for_harmonization"])

    def test_wrong_header_hard_fails(self):
        blob=payload()
        uncompressed=gzip.decompress(blob).decode().replace("LOG(OR)_SE","SE") 
        with write_blob(gzip.compress(uncompressed.encode())) as f:
            with self.assertRaisesRegex(ValueError,"header differs"):
                m.scan_all_rows(f)

    def test_same_small_source_downloads_once_only_and_checks_SHA(self):
        p=payload()
        with patch.object(m,"EXPECTED_SHA256",hashlib.sha256(p).hexdigest()):
            d=m.original_single_get(opener=lambda req,timeout:RemoteResponse(p),
                                    expected_size=len(p))
        self.assertEqual(d["original_source_download_attempts"],1)
        self.assertTrue(d["pinned_prior_original_sha256_match"])
        self.assertEqual(d["original_compressed_bytes"],len(p))
        self.assertTrue(d["original_full_file_deleted_after_scan"])
        self.assertFalse(d["source_file_or_variant_effects_uploaded"])

    def test_tampered_publisher_source_hash_rejected(self):
        p=payload()
        with self.assertRaisesRegex(ValueError,"SHA256 drift"):
            m.original_single_get(opener=lambda req,timeout:RemoteResponse(p),
                                  expected_size=len(p))

    def test_wrong_compressed_length_fails_closed(self):
        p=payload()
        with self.assertRaisesRegex(ValueError,"content length"):
            m.original_single_get(opener=lambda req,timeout:RemoteResponse(p),
                                  expected_size=len(p)+1)

    def test_offsite_redirect_fails_closed(self):
        p=payload()
        with self.assertRaisesRegex(ValueError,"off-site redirect"):
            m.original_single_get(opener=lambda req,timeout:RemoteResponse(p,url="https://evil.example/sumstats"),
                                  expected_size=len(p))

    def test_missing_or_invalid_observation_count_never_neff(self):
        lines=[variant(i,obs="99999999") for i in range(30)]
        with write_blob(payload(lines)) as f:
            r=m.scan_all_rows(f)
        self.assertEqual(r["row_level_QC_counts"]["OBS_CT_exceeds_phenotype_page_total"],30)
        self.assertFalse(r["ready_for_ldsc"])

if __name__=="__main__":
    unittest.main()
