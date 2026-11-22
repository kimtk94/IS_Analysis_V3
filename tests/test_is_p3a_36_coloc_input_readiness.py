"""Synthetic negative controls for the 36-assay source-only readiness gate."""
import copy
import json
import unittest
from pathlib import Path

from scripts.is.audit_p3a_36_coloc_input_readiness import audit, REQUIRED_ABF

MANIFEST = Path("docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json")


def fixtures():
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    source = manifest["groups"]
    recovery_groups = []
    for x in source:
        recovery_groups.append({
            "region_id": x["region_id"], "QTL_dataset": x["QTL_dataset"],
            "gene_id": x["gene_id"],
            "source_BH05_eGene_assay":
                "YES" if x["official_studywide_eGene_BH05_assay"] else
                "NO_SOURCE_BH05_MATCH",
            "complete_cis_tested_SNP_from_BH05_source":
                x["official_full_cis_tested_variant_count"],
            "source_stat_group_membership":
                "UNVERIFIED_25_OF_36_AGGREGATE_ONLY",
        })
    recovery = {
        "schema": "is-p3a-restored-36-pinned-legacy-source-audit-v1",
        "technical_status": "PASS_36_REAL_GROUP_IDS_CROSSCHECKED_TO_THREE_PINNED_ARCHIVED_TABLES",
        "scientific_status": "HOLD", "restored_group_key_count": 36,
        "distinct_gene_id_count": 20, "public_BH05_eGene_assay_count": 25,
        "fully_approved_new_P3A_coloc_groups": 0,
        "new_P3A_coloc_computations": 0, "groups": recovery_groups,
        "source_only_25_group_fullcis_coverage": {
            "sum_tested_fullcis_rows": 142196,
            "sum_BBJ_EAS_nonpal_exact_gene_assay_rows": 59090,
            "sum_nonpal_QTL_EAS_AF_delta_le_0p1_gene_assay_rows": 31313,
        },
    }
    original = {
        "schema": "is-p3a-chr10-original-gwas-combined-36-qtl-groups-v1",
        "technical_status": "PASS_BOTH_PUBLISHER_ORIGINAL_GWAS_CHR10_CONTEXT_SCAN",
        "scientific_status": "HOLD",
        "publisher_accessions": {
            "EUR": "GCST90104540", "EAS": "GCST90104545",
        },
        "EUR_chr10_context_valid_original_variant_rows": 4397,
        "EAS_chr10_context_valid_original_variant_rows": 3913,
        "real_group_keys": 36, "real_gene_ids": 20,
        "distinct_analysed_GWAS_region_count": 1,
        "fully_approved_per_group_gwas_qtl_source_matches": 0,
        "case_control_N_approved": False, "coloc_abf_executed": 0,
        "coloc_susie_executed": 0,
        "per_group_status": [
            {"region_id": x["region_id"], "QTL_dataset": x["QTL_dataset"],
             "gene_id": x["gene_id"],
             "QTL_GWAS_per_variant_same_source_effect_verified": False}
            for x in source
        ],
    }
    eas24 = {
        "schema": "is-gws-eas-original-24-variant-cross-region-readback-v1",
        "source_audit_status": "PASS_ORIGINAL_ROWS_RECONCILED_TO_PINNED_SOURCE_JSON",
        "scientific_status": "HOLD", "original_source": "GCST90104545",
        "nonpalindromic_leads_tested": 24,
        "signed_beta_p_eaf_reproduced": 24, "concordant_chromosome_count": 16,
        "P3A_QTL_gene_assay_count": 36,
        "P3A_exact_QTL_GWAS_source_rows_paired": 0,
        "P3A_coloc_abf": 0, "P3A_coloc_susie": 0,
        "source_GWAS_case_fraction_N_approved": False,
        "held_palindromic_leads": [{"region_id": "IS_XDATA_G0015"}],
    }
    return [manifest, recovery, original, eas24]


class P3A36ReadinessTest(unittest.TestCase):
    def test_valid_36_real_groups_remain_hold(self):
        result = audit(*fixtures())
        self.assertEqual((result["assay_groups"],
                          result["candidate_gene_ids"],
                          result["public_BH05_eGene_assays"]), (36, 20, 25))
        self.assertEqual(result["P3A_36_formal_coloc_abf_eligible"], 0)
        self.assertEqual(result["P3A_36_coloc_susie_executed"], 0)
        self.assertEqual(len(result["genes"]), 20)
        self.assertTrue(all(len(x["coloc_abf_missing_requirements"]) ==
                            len(REQUIRED_ABF) for x in result["groups"]))

    def test_source_membership_not_equal_to_public_BH05(self):
        m, r, o, e = fixtures()
        r["groups"][0]["source_stat_group_membership"] = "SOURCE_POSITIVE"
        with self.assertRaises(ValueError):
            audit(m, r, o, e)

    def test_cannot_pretend_h4_or_coloc_performed(self):
        m, r, o, e = fixtures()
        o["coloc_abf_executed"] = 36
        with self.assertRaises(ValueError):
            audit(m, r, o, e)

    def test_no_legacy_source_manifest_promotion(self):
        m, r, o, e = fixtures()
        m["groups"][0]["approved_for_new_coloc"] = True
        with self.assertRaises(ValueError):
            audit(m, r, o, e)

    def test_same_36_group_keys_required_across_ci(self):
        m, r, o, e = fixtures()
        o["per_group_status"][0]["gene_id"] = "ENSG_FAKE"
        with self.assertRaises(ValueError):
            audit(m, r, o, e)

    def test_no_wrong_qtl_coverage(self):
        m, r, o, e = fixtures()
        r["groups"][0]["complete_cis_tested_SNP_from_BH05_source"] = 999
        with self.assertRaises(ValueError):
            audit(m, r, o, e)

    def test_no_24_lead_case_control_approval(self):
        m, r, o, e = fixtures()
        e["source_GWAS_case_fraction_N_approved"] = True
        with self.assertRaises(ValueError):
            audit(m, r, o, e)

    def test_chr10_palindrome_must_remain_held(self):
        m, r, o, e = fixtures()
        e["held_palindromic_leads"] = []
        with self.assertRaises(ValueError):
            audit(m, r, o, e)

    def test_no_fake_independent_loci(self):
        m, r, o, e = fixtures()
        o["distinct_analysed_GWAS_region_count"] = 36
        with self.assertRaises(ValueError):
            audit(m, r, o, e)

    def test_missing_eas24_audit_cannot_pass(self):
        m, r, o, e = fixtures()
        e["signed_beta_p_eaf_reproduced"] = 0
        with self.assertRaises(ValueError):
            audit(m, r, o, e)


if __name__ == "__main__":
    unittest.main()
