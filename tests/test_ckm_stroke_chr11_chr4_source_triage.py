"""Chr11 exact variant identity and chr4 P/Wald source precision: no-network tests."""
from __future__ import annotations
import math
import unittest
from scripts.ckm_stroke.audit_gigastroke_chr11_chr4_source_triage import (
    original_chr_target_lines, targeted_records, check_against_original_chr_evidence,
    p_compatible_with_source_rounding, rounded_decimal_interval,
)
from scripts.ckm_stroke.audit_gigastroke_eaf_gc import diagnostic,summarize,DIAGNOSTIC_KEYS
from scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization import COLUMNS

def make(ch,pos,eaf="0.2",beta="0.1",se="0.05",p="0.0455003",ea="A",oa="C"):
    return [str(ch),str(pos),eaf,beta,se,p,
            "1.1051709180756477","1.0","1.4",ea,oa]

def dummy():
    panel={(4,100):("rs-chr4-1","A","C"),
           (11,200):("rs-chr11-1","A","C"),
           (11,300):("rs-chr11-duplicated","A","C")}
    hrs={(4,100):0.2,(11,200):0.85,(11,300):0.5}
    body=[
        make(4,100,p="0.03"),
        make(11,200,eaf="0.20"),
        make(11,300),make(11,300),
        make(1,123456),
    ]
    lines=iter(["\t".join(COLUMNS)+"\n"]+
               ["\t".join(r)+"\n" for r in body])
    n,kept,dup,valid=original_chr_target_lines(lines,panel,hrs)
    return n,kept,dup,valid

def old_reference(retained):
    old={"eaf_gc_source_descriptive":{"by_chr":[{} for _ in range(22)]},
         "chr_attrition_pre_munge":[{} for _ in range(22)]}
    for ch in (4,11):
        cur=summarize({coord:r[0] for coord,r in retained.items() if coord[0]==ch})["by_chr"][ch-1]
        old["eaf_gc_source_descriptive"]["by_chr"][ch-1]=cur
        old["chr_attrition_pre_munge"][ch-1]={
            "allele_aligned_numeric_valid_pre_munge":cur["matched_diagnostic_SNPs"]
        }
    return old

class TargetedDiagnosticTests(unittest.TestCase):
    def test_chr11_outlier_identity_and_chr4_moderate_wald(self):
        n,kept,dups,valid=dummy()
        self.assertEqual((n,dups,valid),(5,1,2))
        source=old_reference(kept)
        check_against_original_chr_evidence("AIS",source,kept)
        f,w=targeted_records(kept,source)
        self.assertEqual(len(f),1)
        self.assertEqual(f[0]["rsid_from_pinned_GIANT"],"rs-chr11-1")
        self.assertEqual(f[0]["position_GRCh37"],200)
        self.assertEqual(f[0]["GIANT_HRS_tested_allele_frequency"],0.85)
        self.assertEqual(f[0]["GIGASTROKE_frequency_oriented_to_GIANT"],0.2)
        self.assertEqual(len(w),1)
        self.assertEqual(w[0]["position_GRCh37"],100)
        self.assertTrue(w[0]["abs_log10P_difference"]>=0.1)
        self.assertFalse(w[0]["source_association_effects_exported"])

    def test_chromosome_exact_recount_fail_closed(self):
        _,records,_,_=dummy()
        previous=old_reference(records)
        previous["eaf_gc_source_descriptive"]["by_chr"][10]["categories"]["hrs_absdiff_ge_020"]+=1
        with self.assertRaisesRegex(ValueError,"diagnostic category drift"):
            check_against_original_chr_evidence("AIS",previous,records)

    def test_candidate_duplicated_position_quarantined(self):
        _,records,duplicates,_=dummy()
        self.assertEqual(duplicates,1)
        self.assertNotIn((11,300),records)

    def test_printed_beta_se_p_rounding_interval(self):
        lo,hi=rounded_decimal_interval("1.234e-3")
        self.assertTrue(lo<0.001234<hi)
        self.assertTrue(p_compatible_with_source_rounding(make(4,100,p="0.0455003")))
        self.assertFalse(p_compatible_with_source_rounding(make(4,100,beta="0.100000000000",se="0.050000000000",p="0.03")))

    def test_p_zero_has_no_rounding_claim(self):
        self.assertIsNone(p_compatible_with_source_rounding(make(4,100,p="0")))

    def test_original_header_must_match(self):
        with self.assertRaisesRegex(ValueError,"header drift"):
            original_chr_target_lines(iter(["chr\tpos\n"]),{},{})

    def test_outlier_counts_match_full_reference_else_hold(self):
        _,retained,_,_=dummy()
        ref=old_reference(retained)
        ref["eaf_gc_source_descriptive"]["by_chr"][3]["categories"]["abs_log10p_delta_010_to_1"]=0
        with self.assertRaisesRegex(ValueError,"outlier count drift"):
            targeted_records(retained,ref)

if __name__=="__main__":
    unittest.main()
