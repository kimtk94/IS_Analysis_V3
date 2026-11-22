"""Original 30-region archived source and 24-source lead-signed effect guard tests."""
from __future__ import annotations
import copy
import csv
import gzip
import importlib
import json
from pathlib import Path
import tempfile
import unittest

M=importlib.import_module("scripts.is.audit_gws_eas_24_leads_original")
LEGACY=Path("legacy/docs/data/is/IS_30_REGIONS_SIX_OUTCOMES_GWAS_CONTEXT_500KB.tsv")

def source():
    return M.read_legacy_180(LEGACY)

def official_row(marker,reverse=False):
    alt=marker["historical_alt_allele"]
    ref=marker["historical_ref_allele"]
    eaf=marker["historical_ALT_EAF"]
    return dict(chromosome=marker["chromosome"],
                base_pair_location=str(marker["position"]),
                effect_allele=ref if reverse else alt,
                other_allele=alt if reverse else ref,
                beta=str(-marker["historical_beta_ALT"] if reverse else marker["historical_beta_ALT"]),
                standard_error="0.02",
                p_value=str(marker["historical_p"]),
                effect_allele_frequency=str(1-eaf if reverse else eaf))

def fixture(path,records,reverse_n=0,extra=None):
    rows=[official_row(m,i<reverse_n) for i,m in enumerate(records)]
    rows.extend(extra or [])
    with gzip.open(path,"wt",encoding="utf-8",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=M.HEADER,delimiter="\t")
        wr.writeheader()
        wr.writerows(rows)

class True24NonpalArchiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=source()

    def test_pinned_actual_30_source_regions_partition(self):
        d=self.d
        self.assertEqual(d["historical_rows"],180)
        self.assertEqual(len(d["selected"]),24)
        self.assertEqual(len(d["withheld"]),4)
        self.assertEqual(len(d["missing"]),2)
        self.assertEqual(len({s["region_id"] for s in d["selected"]+
                        d["withheld"]+d["missing"]}),30)
        self.assertTrue(all(x["GRCh37_variant_key"] for x in d["selected"]))
        self.assertTrue(all(x["status"]=="PALINDROMIC_STRAND_HOLD" for x in d["withheld"]))

    def test_real_archived_region_15_palindrome_excluded(self):
        self.assertFalse(any(x["region_id"]=="IS_XDATA_G0015"
                             for x in self.d["selected"]))
        self.assertIn("IS_XDATA_G0015",
            {x["region_id"] for x in self.d["withheld"]})

    def test_source_normalization_effect_allele(self):
        r=self.d["selected"][0]
        first=M.inspect_original_row(official_row(r))
        flip=M.inspect_original_row(official_row(r,reverse=True))
        self.assertEqual(M.compare(r,[first])["status"],
                         "SIGNED_ORIGINAL_EAS_LEAD_REPRODUCED_ARCHIVE")
        self.assertEqual(M.compare(r,[flip])["status"],
                         "SIGNED_ORIGINAL_EAS_LEAD_REPRODUCED_ARCHIVE")
        self.assertFalse(M.compare(r,[flip])["effect_allele_was_historical_ALT"])

    def test_effect_numeric_mismatch_is_visible(self):
        r=self.d["selected"][0]
        wrong=official_row(r);wrong["beta"]=str(r["historical_beta_ALT"]+0.2)
        v=M.compare(r,[M.inspect_original_row(wrong)])
        self.assertEqual(v["status"],"ORIGINAL_EAS_LEAD_NUMERIC_DISCORDANCE_HOLD")
        self.assertFalse(v["joint_effect_p_eaf_reproduced"])

    def test_af_flip_required_for_ref_effect(self):
        r=self.d["selected"][0]
        wrong=official_row(r,reverse=True)
        wrong["effect_allele_frequency"]=str(r["historical_ALT_EAF"])
        v=M.compare(r,[M.inspect_original_row(wrong)])
        if abs(2*r["historical_ALT_EAF"]-1)>M.AF_ABS_TOL:
            self.assertFalse(v["ALT_oriented_eaf_match"])

    def test_allele_not_same_pair_is_hold(self):
        r=self.d["selected"][0]
        wrong=official_row(r)
        wrong["effect_allele"]="A" if "A" not in (r["historical_ref_allele"],r["historical_alt_allele"]) else "G"
        if wrong["effect_allele"] in (r["historical_ref_allele"],r["historical_alt_allele"]):
            wrong["effect_allele"]="T" 
        v=M.compare(r,[M.inspect_original_row(wrong)])
        self.assertEqual(v["status"],"AMBIGUOUS_OR_INCOMPATIBLE_ORIGINAL_ALLELES_AT_POSITION")

    def test_original_indel_not_accepted_as_SNV(self):
        r=official_row(self.d["selected"][0])
        r["effect_allele"]="AC"
        self.assertIsNone(M.inspect_original_row(r))

    def test_duplicated_original_position_is_not_trusted(self):
        r=self.d["selected"][0]
        first=M.inspect_original_row(official_row(r))
        v=M.compare(r,[first,first])
        self.assertEqual(v["status"],"AMBIGUOUS_OR_INCOMPATIBLE_ORIGINAL_ALLELES_AT_POSITION")

    def test_absent_original_position_does_not_fake_zero_p(self):
        a=M.compare(self.d["selected"][0],[])
        self.assertEqual(a["status"],"OFFICIAL_LEAD_POSITION_ABSENT")
        self.assertFalse(a["joint_effect_p_eaf_reproduced"])

    def test_full_synthetic_24_SNV_scan_preserves_36_QTL_hold(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"synthetic.tsv.gz"
            fixture(path,self.d["selected"],reverse_n=4)
            r=M.scan_original_file(path,self.d,min_rows=24)
        self.assertEqual(r["original_nonpalindromic_archived_leads_tested"],24)
        self.assertEqual(r["original_signed_beta_p_eaf_all_three_match_count"],24)
        self.assertEqual(r["original_same_position_allele_pair_match_count"],24)
        self.assertEqual(r["technical_source_archive_signed_reproduction_screen"],
                         "PASS_REVIEW_REQUIRED")
        self.assertFalse(r["published_original_study_Ncase_Ncontrol_independently_approved"])
        self.assertEqual(r["coloc_abf_executed"],0)
        self.assertEqual(r["coloc_susie_executed"],0)

    def test_synthetic_fewer_than_20_matching_leads_cannot_pass(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"synthetic.tsv.gz"
            fixture(path,self.d["selected"],extra=[])
            with gzip.open(path,"rt") as f:rows=list(csv.DictReader(f,delimiter="\t"))
            for row in rows[:7]:row["p_value"]="0.99"
            with gzip.open(path,"wt") as f:
                wr=csv.DictWriter(f,fieldnames=M.HEADER,delimiter="\t")
                wr.writeheader();wr.writerows(rows)
            r=M.scan_original_file(path,self.d,min_rows=24)
        self.assertEqual(r["original_signed_beta_p_eaf_all_three_match_count"],17)
        self.assertEqual(r["technical_source_archive_signed_reproduction_screen"],
                         "HOLD_SOURCE_REPRODUCTION_INSUFFICIENT")

    def test_old_Git_blob_fingerprint_prevents_unreviewed_edits(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"old.tsv"
            f.write_bytes(LEGACY.read_bytes()+b"\n")
            with self.assertRaisesRegex(ValueError,"Git blob SHA drift"):
                M.read_legacy_180(f)

    def test_source_header_missing_SE_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"s.tsv.gz"
            with gzip.open(f,"wt") as z:
                z.write("chromosome\tbase_pair_location\tbeta\tp_value\n")
            with self.assertRaisesRegex(ValueError,"mandatory signed effect"):
                M.scan_original_file(f,self.d,min_rows=1)

if __name__=="__main__":
    unittest.main()
