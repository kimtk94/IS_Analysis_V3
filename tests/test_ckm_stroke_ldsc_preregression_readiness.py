"""Fail-closed pre-LDSC descriptors, no LD score regression or participant claims."""
from __future__ import annotations
import gzip
import hashlib
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from scripts.ckm_stroke.audit_ldsc_preregression_readiness import (
    examine_candidate, outcome_gates, EXPECTED_OUTCOMES
)


def fake_candidate(root:Path,name:str,rows:list):
    head=("SNP\tN_REPORTED\tZ_FROM_BETA_SE\tA1\tA2\n"
          if name=="locke_2015" else "SNP\tN\tZ\tA1\tA2\n")
    content=(head+"".join("\t".join(map(str,row))+"\n" for row in rows)).encode()
    path=root/(name+".gz");path.write_bytes(gzip.compress(content,mtime=0))
    return path,{"rows":len(rows),"plaintext_sha256":hashlib.sha256(content).hexdigest()}


def fake_outcomes():
    cohort={"technical_status":"PASS_PRIMARY_CORRECTED_COHORT_AND_TABLE3_ONLY",
            "ready_for_ldsc":False,"outcomes":[]}
    prefix={"technical_status":"PASS_FOUR_BOUNDED_LOGOR_EFFECT_SCALE_ONLY",
            "ready_for_ldsc":False,"prefixes":[]}
    for name,(acc,cases) in EXPECTED_OUTCOMES.items():
        cohort["outcomes"].append({
            "phenotype":name,"accession":acc,
            "EUR_case_count_from_original_table1":cases,
            "EUR_control_count_from_original_table1":1234808,
            "original_CovLDSC_h2_SE_ratio":1.67 if name=="LAA" else 9.8,
        })
        prefix["prefixes"].append({
            "phenotype":name,"accession":acc,
            "source_per_variant_N_available":False,
            "observed":{"parsed_effect_rows":510},
            "maximum_abs_beta_minus_lnOR":5e-15,
        })
    overlap={
        "technical_status":"PASS_FIRSTPARTY_TWO_SUPPLEMENT_COHORT_NAME_RECONCILIATION",
        "exact_individual_overlap_count_identified":False,
        "direct_cohort_label_matches":10,
    }
    return cohort,prefix,overlap


class TestLDSCPreregression(unittest.TestCase):
    def test_2015_fractional_n_preserved_and_positive_z(self):
        with TemporaryDirectory() as d:
            path,info=fake_candidate(Path(d),"locke_2015",[
                ("rs1","234000.5","1.5","A","C"),
                ("rs2","250000","-2","G","T"),
                ("rs3","220000","0","C","T"),
            ])
            r=examine_candidate("locke_2015",path,info,strict=False)
            self.assertEqual(r["SNPs"],3)
            self.assertEqual(r["reported_N"]["fractional_N_rows"],1)
            self.assertEqual(r["descriptive_only"]["mean_chi2"],round((2.25+4)/3,8))
            self.assertFalse(r["ready_for_ldsc"])

    def test_2018_integer_n_and_chisq(self):
        with TemporaryDirectory() as d:
            path,info=fake_candidate(Path(d),"yengo_2018",[
                ("rs1","200000","2","A","C"),
                ("rs2","300000","-2","G","T"),
            ])
            r=examine_candidate("yengo_2018",path,info,strict=False)
            self.assertEqual(r["descriptive_only"]["mean_chi2"],4.0)
            self.assertEqual(r["reported_N"]["fractional_N_rows"],0)

    def test_2018_fractional_n_fails_closed(self):
        with TemporaryDirectory() as d:
            path,info=fake_candidate(Path(d),"yengo_2018",[
                ("rs1","200000.1","2","A","C"),
            ])
            with self.assertRaisesRegex(ValueError,"integral source N"):
                examine_candidate("yengo_2018",path,info,strict=False)

    def test_source_digest_mismatch_rejected_even_in_fixture_mode(self):
        with TemporaryDirectory() as d:
            path,info=fake_candidate(Path(d),"locke_2015",[
                ("rs1","234000","1.5","A","C")
            ])
            info["plaintext_sha256"]="0"*64
            with self.assertRaisesRegex(ValueError,"SHA/row"):
                examine_candidate("locke_2015",path,info,strict=False)

    def test_repeated_rsid_rejected(self):
        with TemporaryDirectory() as d:
            path,info=fake_candidate(Path(d),"locke_2015",[
                ("rs1","234000","2","A","C"),
                ("rs1","234000","2","A","C"),
            ])
            with self.assertRaisesRegex(ValueError,"duplicate"):
                examine_candidate("locke_2015",path,info,strict=False)

    def test_nonfinite_z_rejected(self):
        with TemporaryDirectory() as d:
            path,info=fake_candidate(Path(d),"locke_2015",[
                ("rs1","234000","inf","A","C")
            ])
            with self.assertRaisesRegex(ValueError,"nonfinite"):
                examine_candidate("locke_2015",path,info,strict=False)

    def test_published_h2_benchmark_does_not_authorize_LDSC(self):
        c,p,o=fake_outcomes()
        result=outcome_gates(c,p,o)
        self.assertEqual(len(result),4)
        self.assertEqual(result[1]["phenotype"],"LAA")
        self.assertTrue(result[1]["published_h2_precision_below_precaution_ratio4"])
        self.assertTrue(all(x["input_gate"].startswith("HOLD") for x in result))
        self.assertTrue(all(x["allow_constrain_cross_trait_intercept_to_zero"] is False for x in result))

    def test_swapped_LAA_CE_original_accession_rejected(self):
        c,p,o=fake_outcomes()
        for item in p["prefixes"]:
            if item["phenotype"]=="LAA":item["accession"]="GCST90104541"
        with self.assertRaisesRegex(ValueError,"source phenotype"):
            outcome_gates(c,p,o)

    def test_quantified_overlap_flag_not_permitted_without_new_policy(self):
        c,p,o=fake_outcomes()
        o["exact_individual_overlap_count_identified"]=True
        with self.assertRaisesRegex(ValueError,"cohort-overlap evidence"):
            outcome_gates(c,p,o)

    def test_per_snp_N_claim_not_permitted(self):
        c,p,o=fake_outcomes()
        p["prefixes"][0]["source_per_variant_N_available"]=True
        with self.assertRaisesRegex(ValueError,"source phenotype"):
            outcome_gates(c,p,o)


if __name__=="__main__":
    unittest.main()
