"""Offline source-evidence review tests; no original GWAS or network access."""
from __future__ import annotations
import copy
import unittest
from scripts.ckm_stroke.review_gigastroke_full_evidence import (
    PINNED, OBSERVED_COUNTS, EXPECTED_BMI_CANDIDATE, EXPECTED_BMI_SHA, review
)


def fixture():
    out={
        "schema":"ckm-stroke-gigastroke-eur-full-md5-bmi-harmonization-v1",
        "technical_status":"PASS_ALL_REQUESTED_ORIGINAL_MD5_AND_ALLELE_COUNTS",
        "scientific_status":"HOLD","ready_for_ldsc":False,
        "LDSC_h2_or_rg_executed":False,"source_archives_retained":False,
        "BMI_reference_panel":{
            "candidate_total":EXPECTED_BMI_CANDIDATE,
            "candidate_plaintext_sha256":EXPECTED_BMI_SHA
        },
        "outcomes":[]
    }
    for phenotype,(accession,md5,compressed_bytes) in PINNED.items():
        rows,aligned,unmatched=OBSERVED_COUNTS[phenotype]
        out["outcomes"].append({
            "phenotype":phenotype,"accession":accession,
            "original_compressed_bytes":compressed_bytes,
            "published_md5":md5,"raw_md5_verified":True,
            "publisher_yaml_md5_sidecar_verified":True,
            "technical_status":"PASS_RAW_MD5_AND_ALLELE_CLASSIFICATION_ONLY",
            "source_sha256":"a"*64,
            "outcome_has_rsid_column":False,
            "outcome_per_variant_N_available":False,
            "counts":{
                "original_outcome_rows":rows,
                "outcome_rows_with_BMI_coordinate":aligned,
                "outcome_rows_absent_from_BMI_candidate":rows-aligned,
                "candidate_coordinates_not_in_outcome":unmatched,
                "outcome_collision_coordinates_excluded":0,
                "unique_outcome_coordinate_matches":aligned,
                "unique_harmonizable_BMI_positions":aligned,
                "unique_same":aligned,
                "unique_orientation_sign_positive":aligned,
                "unique_orientation_sign_negative":0,
            }
        })
    return out


class ReviewTests(unittest.TestCase):
    def test_four_original_aggregates_pass_not_scientific_go(self):
        r=review(fixture())
        self.assertEqual(r["technical_status"],"PASS_PUBLISHER_MD5_AND_AGGREGATE_INVARIANTS")
        self.assertEqual(len(r["outcomes"]),4)
        self.assertEqual(r["outcomes"][1]["matched_and_allele_aligned"],956_929)
        self.assertEqual(r["outcomes"][2]["candidate_not_in_outcome"],11_237)
        self.assertTrue(r["uniform_allele_orientation_needs_independent_review"])
        self.assertFalse(r["ready_for_ldsc"])
        self.assertTrue(r["illustrative_Neff_is_NOT_per_variant_or_approved_for_LDSC"])

    def test_duplicate_outcome_rejected(self):
        data=fixture()
        data["outcomes"][-1]["phenotype"]="AIS"
        with self.assertRaisesRegex(ValueError,"duplicate or missing"):
            review(data)

    def test_false_publisher_md5_rejected(self):
        data=fixture()
        data["outcomes"][2]["raw_md5_verified"]=False
        with self.assertRaisesRegex(ValueError,"MD5"):
            review(data)

    def test_nonheld_science_rejected(self):
        data=fixture()
        data["ready_for_ldsc"]=True
        with self.assertRaisesRegex(ValueError,"scientific"):
            review(data)

    def test_broken_row_partition_rejected(self):
        data=fixture()
        data["outcomes"][0]["counts"]["outcome_rows_absent_from_BMI_candidate"]-=1
        with self.assertRaisesRegex(ValueError,"partition"):
            review(data)

    def test_broken_allele_orientation_partition_rejected(self):
        data=fixture()
        data["outcomes"][3]["counts"]["unique_same"]-=1
        with self.assertRaisesRegex(ValueError,"directional"):
            review(data)

    def test_mismatched_candidate_sha_rejected(self):
        data=fixture()
        data["BMI_reference_panel"]["candidate_plaintext_sha256"]="0"*64
        with self.assertRaisesRegex(ValueError,"candidate changed"):
            review(data)

    def test_unexpected_source_row_count_rejected(self):
        data=fixture()
        data["outcomes"][1]["counts"]["original_outcome_rows"]+=1
        data["outcomes"][1]["counts"]["outcome_rows_absent_from_BMI_candidate"]+=1
        with self.assertRaisesRegex(ValueError,"original source counts"):
            review(data)


if __name__=="__main__":
    unittest.main()
