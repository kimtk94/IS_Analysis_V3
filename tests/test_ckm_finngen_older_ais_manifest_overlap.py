"""Public FinnGen R9/R10/R13 AIS availability and no-overlap-claim unit gates."""
from __future__ import annotations
import hashlib
import unittest
from unittest.mock import patch
from scripts.ckm_stroke import audit_finngen_older_ais_manifest_overlap as m

def tab(release="R9",more=True):
    first="phenocode\tphenotype\tnum_cases\tnum_controls\tpath_https\n"
    lines=[
      "I9_STR_EXH\tIschaemic stroke excluding haemorrhage\t9000\t220000\thttps://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_I9_STR_EXH.gz\n",
      "I9_STR_EMBOLIC\tEmbolic stroke\t1100\t220000\thttps://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_I9_STR_EMBOLIC.gz\n",
      "I9_STR\tStroke including ICH\t10000\t221000\thttps://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_I9_STR.gz\n"
    ]
    if release=="R13":
        lines=lines[1:]
        lines=[x.replace("r9","r13").replace("R9","R13") for x in lines]
        lines.append("I9_STR_EMBOLIC_WITH_FLIMMER\tEmbolic stroke with AF\t500\t200000\thttps://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_I9_STR_EMBOLIC_WITH_FLIMMER.gz\n")
        lines+=["DUMMY%04d\tother\t100\t1000\t\n"%i for i in range(2751)]
    return (first+"".join(lines)).encode("utf-8")

class TestOlderReleaseSource(unittest.TestCase):
    def test_r9_manifest_exh_is_release_specific(self):
        src=tab()
        result=m.parse_release("R9",src,hashlib.sha256(src).hexdigest())
        self.assertEqual(result["requested_endpoints"]["I9_STR_EXH"]["release_cases"],9000)
        self.assertTrue(result["requested_endpoints"]["I9_STR_EXH"]["present"])
        self.assertFalse(result["publisher_independent_expected_SHA256_known"])
        self.assertFalse(result["per_SNP_Neff_in_manifest"])
    def test_r13_ischaemic_missing_still_not_data(self):
        src=tab("R13")
        with patch.object(m,"R13_MANIFEST_SHA",hashlib.sha256(src).hexdigest()):
            result=m.parse_release("R13",src,hashlib.sha256(src).hexdigest())
        self.assertFalse(result["requested_endpoints"]["I9_STR_EXH"]["present"])
        self.assertTrue(result["requested_endpoints"]["I9_STR_EMBOLIC"]["present"])
    def test_r13_sha_tamper_is_rejected(self):
        src=tab("R13")
        with self.assertRaisesRegex(ValueError,"SHA256"):
            m.parse_release("R13",src,hashlib.sha256(src).hexdigest())
    def test_manifest_duplicate_candidate_rejected(self):
        src=tab().replace(b"I9_STR\t",b"I9_STR_EXH\t")
        with self.assertRaisesRegex(ValueError,"duplicated"):
            m.parse_release("R9",src,hashlib.sha256(src).hexdigest())
    def test_no_variant_neff_in_release_counts(self):
        records={}
        for release in ("R9","R10","R13"):
            raw=tab("R13" if release=="R13" else "R9")
            records[release]=(raw,hashlib.sha256(raw).hexdigest())
        cohort={"participant_ID_overlap_quantified":False}
        def fake_head(url):
            return {"HTTP_status":403,"remote_file_header_or_variant_Neff_verified":False}
        with patch.object(m,"R13_MANIFEST_SHA",records["R13"][1]):
            result=m.audit(records,cohort,head_checker=fake_head)
        self.assertFalse(result["primary_LDSC_SNP_specific_Neff_resolved"])
        self.assertFalse(result["R9_R10_same_FinnGen_project_independent_replication_between_releases"])
        self.assertFalse(result["MEGASTROKE_GIGASTROKE_participant_overlap_resolved"])
        self.assertEqual(len(result["FinnGen_R13_source_GWAS_HEAD_access_only"]),2)
        self.assertFalse(result["ready_for_ldsc"])
    def test_not_all_release_manifest_keys_is_hard_failure(self):
        with self.assertRaisesRegex(ValueError,"all pinned"):
            m.audit({},cohort={})
    def test_cohort_names_not_proof_of_zero_overlap(self):
        table={i:{"B":"Cohort-%s"%i} for i in range(5,49)}
        table[15]["B"]="UKBiobank"
        table[19]["B"]="Helsinki-Stroke"
        with patch.object(m,"workbook_cells",return_value=(table,{})):
            result=m.cohort_review(b"mock")
        self.assertEqual(result["EUR_cohort_rows"],44)
        self.assertEqual(len(result["explicit_Nordic_Finnish_named_cohort_candidates"]),1)
        self.assertFalse(result["absence_of_FinnGen_label_proves_zero_overlap"])
        self.assertFalse(result["participant_ID_overlap_quantified"])
    def test_suspect_official_redirect_is_rejected(self):
        class Fake:
            status=200
            headers={}
            def __enter__(self):return self
            def __exit__(self,*a):return None
            def geturl(self):return "https://example.com/wrong"
            def read(self,n):return b"x"
        with self.assertRaisesRegex(ValueError,"not official"):
            m.get_text(m.MANIFESTS["R9"],opener=lambda req,timeout:Fake())
    def test_not_an_allowed_head_url(self):
        with self.assertRaisesRegex(ValueError,"not original"):
            m.head_check("https://example.org/secret.tsv")
if __name__=="__main__":
    unittest.main()
