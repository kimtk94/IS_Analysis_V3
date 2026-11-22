"""Candidate-wide P3A regression tests using synthetic group fixtures ONLY."""
from __future__ import annotations
import copy
import importlib
import json
from pathlib import Path
import unittest

m=importlib.import_module("scripts.is.candidate_wide_coloc_source_triage")
BASE=Path("docs/is/evidence/IS_CANDIDATE_WIDE_COLOC_BASELINE_20261122.json")
FOUR=Path("docs/is/evidence/IS_FOUR_LOCUS_PRIVATE_INPUT_QC_20261122.json")

def inputs():
    return json.loads(BASE.read_text(encoding="utf-8")),json.loads(FOUR.read_text(encoding="utf-8"))

def mock_original_group_report():
    groups=[]
    for i in range(36):
        valid=(2500 if i < 24 else (8612 if i==24 else 0))
        groups.append({
            "region_id":"BBJ_IS_L" + str(i//10+1).zfill(3),
            "QTL_dataset":"QTD000131",
            "gene_id":"ENSG"+str(i+1).zfill(11),
            "fullcis_archived_rows":valid+100,
            "source_flagged_stat_rows":valid,
            "coloc_abf_eligible":False,
            "scientific_status":"ABF_FOLLOWUP_FLAG_ONLY_NOT_COLOC_VALIDATION",
            "blocking_reasons":[
                "QTL_GWAS_SAMPLE_SIZES_UNVERIFIED",
                "INDEPENDENT_ALLELE_REFERENCE_PENDING",
                "GWAS_EFFECT_VARIANCE_UNVERIFIED",
            ]
        })
    return {
        "schema":"is-archived-fasta-abf-audit-v1",
        "integrity":"PASS",
        "scientific_status":"ARCHIVED_SOURCE_FLAGS_AND_ABF_READINESS_ONLY",
        "readiness":{"groups":393,"abf_candidate_rows":36},
        "abf_gate_summary":{
            "eligible_for_coloc_abf":0,
            "groups_with_source_stat_rows":25,
            "groups_with_no_source_stat_rows":11,
            "source_stat_rows_total":68612
        },
        "abf_screen":{"screened":36,"coloc_abf_executed":0},
        "p3a_gwas_source_identity":{"accepted_sources":[],"registry_status":"MISSING_APPROVED_REGISTRY"},
        "p3a_source_header_audit":{"field_presence":{
            "bbj_effect_se":False,
            "eas_effect_se":False,
            "bbj_sample_size":False,
            "eas_sample_size":False,
            "stroke_case_fraction":False,
        }},
        "abf_flagged_groups":groups
    }

class BatchCandidateScreenTests(unittest.TestCase):
    def test_2225_and_393_denominators_separate(self):
        b,f=inputs()
        s,rows=m.generate(b,f,None)
        self.assertEqual(s["candidate_total_positional_genes"],2225)
        self.assertEqual(s["readiness_assay_groups"],393)
        self.assertEqual(s["p3a_source_flagged_assay_groups"],36)
        self.assertEqual(s["p3a_groups_with_basic_stats"],25)
        self.assertEqual(s["technical_status"],"PASS_SOURCE_AGGREGATE_ONLY_NO_GROUP_IDS")
        self.assertFalse(s["all_36_original_ids_loaded"])
        self.assertEqual(rows,[])
        self.assertIsNone(s["unique_genes_in_36_assay_groups"])
        self.assertEqual(s["p3a_science_approved_groups"],0)

    def test_synthetic_full_panel_still_no_coloc(self):
        b,f=inputs()
        s,rows=m.generate(b,f,mock_original_group_report())
        self.assertEqual(len(rows),36)
        self.assertEqual(sum(r["source_flagged_stat_rows"] for r in rows),68612)
        self.assertEqual(sum(r["source_flagged_stat_rows"]>0 for r in rows),25)
        self.assertEqual(s["unique_genes_in_36_assay_groups"],36)
        self.assertTrue(all(r["coloc_abf_eligible"] is False for r in rows))
        self.assertTrue(all(r["scientific_status"] == "HOLD" for r in rows))
        self.assertEqual(s["p3a_coloc_executed_in_this_run"],0)

    def test_artificial_group_with_missing_metadata_cannot_be_approved(self):
        b,f=inputs()
        g=mock_original_group_report()
        g["abf_flagged_groups"][0]["coloc_abf_eligible"]=True
        with self.assertRaisesRegex(ValueError,"falsely called eligible"):
            m.generate(b,f,g)

    def test_duplicates_in_gene_tissue_region_not_allowed(self):
        b,f=inputs()
        g=mock_original_group_report()
        g["abf_flagged_groups"][1]["gene_id"]=g["abf_flagged_groups"][0]["gene_id"]
        g["abf_flagged_groups"][1]["region_id"]=g["abf_flagged_groups"][0]["region_id"]
        with self.assertRaisesRegex(ValueError,"duplicate"):
            m.generate(b,f,g)

    def test_valid_source_stats_does_not_equal_25_unique_genes(self):
        b,f=inputs()
        g=mock_original_group_report()
        g["abf_flagged_groups"][25]["gene_id"]=g["abf_flagged_groups"][0]["gene_id"]
        s,_=m.generate(b,f,g)
        self.assertEqual(s["unique_genes_in_36_assay_groups"],35)
        self.assertEqual(s["p3a_groups_with_basic_stats"],25)

    def test_source_counts_are_not_gene_counts(self):
        b,f=inputs()
        b["candidate_universe"]["public_qtl_nominal_available"]=250
        with self.assertRaisesRegex(ValueError,"partition invalid"):
            m.generate(b,f,None)

    def test_unknown_4_legacy_must_not_join_new_p3a(self):
        b,f=inputs()
        f["scientific_status"]="PASS"
        with self.assertRaisesRegex(ValueError,"four-locus"):
            m.generate(b,f,None)

    def test_no_manually_adopted_gwas_n(self):
        b,f=inputs()
        g=mock_original_group_report()
        g["p3a_gwas_source_identity"]["accepted_sources"]=["BBJ"]
        with self.assertRaisesRegex(ValueError,"GWAS source"):
            m.generate(b,f,g)

    def test_qtl_sample_count_in_source_missing_gates(self):
        b,f=inputs()
        g=mock_original_group_report()
        g["p3a_source_header_audit"]["field_presence"]["bbj_sample_size"]=True
        with self.assertRaisesRegex(ValueError,"source lacks"):
            m.generate(b,f,g)

    def test_no_fake_0_when_real_group_ids_missing(self):
        b,f=inputs()
        s,_=m.generate(b,f,None)
        self.assertIsNone(s["unique_genes_in_36_assay_groups"])
        self.assertEqual(s["per_group_rows_output"],0)

    def test_missing_original_archive_report_is_not_silent_36_group_success(self):
        b,f=inputs()
        g=mock_original_group_report()
        g["abf_flagged_groups"]=g["abf_flagged_groups"][:-1]
        with self.assertRaisesRegex(ValueError,"36 actual"):
            m.generate(b,f,g)

    def test_source_candidate_schema_and_baseline_consistent(self):
        b,f=inputs()
        g=mock_original_group_report()
        g["abf_flagged_groups"][0]["source_flagged_stat_rows"] += 1
        with self.assertRaisesRegex(ValueError,"68,612"):
            m.generate(b,f,g)


if __name__=="__main__":
    unittest.main()
