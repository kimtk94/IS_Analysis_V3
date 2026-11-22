"""Synthetic source-join tests with immutable original historical GWAS region table."""
from __future__ import annotations
import copy
import importlib
import json
from pathlib import Path
import tempfile
import unittest

m=importlib.import_module("scripts.is.reconcile_p3a_chr10_actual_gwas_to_legacy")
s=importlib.import_module("scripts.is.audit_p3a_chr10_original_publisher_gwas")
OLD=Path("legacy/docs/data/is/IS_30_REGIONS_SIX_OUTCOMES_GWAS_CONTEXT_500KB.tsv")
MAN=Path("docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json")

def fake(ancestry):
    c={"scientific_status":"HOLD","source_GWAS_effect_and_SE_original_rows_checked":True,
       "source_GWAS_effects_reconciled_to_archived_fullcis_36_groups":False,
       "source_context_GRCh37":dict(s.REGION),
       "source_context_valid_variant_rows":3913 if ancestry=="EAS" else 4397,
       "original_region_valid_variant_rows":2110,
       "source_context_p_lt_5e8_rows":69 if ancestry=="EAS" else 8,
       "source_context_min_p_descriptive":5.584e-11 if ancestry=="EAS" else 1e-9,
       "legacy_GWS_EAS_AIS_lead_comparison":{
           "status":"PALINDROMIC_LEAD_NUMERIC_CONCORDANCE_STRAND_HOLD",
           "source_identity_approved":False,
           "n_independent_original_variant_matches_for_source_approval":0} if ancestry=="EAS" else None}
    groups=s.build_group_rows(s.checked_manifest(MAN),c,ancestry)
    return {
        "schema":"is-p3a-chr10-original-publisher-gwas-all-36-v1",
        "technical_status":"PASS_PUBLISHER_WHOLE_FILE_AND_CHR10_ORIGINAL_ROWS",
        "scientific_status":"HOLD","ancestry":ancestry,
        "source":{"publisher_md5_and_prior_sha256_verified":True,
                  "accession":s.GATE.SOURCES[ancestry]["accession"]},
        "context":c,"actual_36_assay_groups_considered":36,
        "unique_candidate_gene_ids":20,
        "new_coloc_performed":0,
        "independent_GWAS_genomic_loci_replicated":0,
        "all_group_input_statuses":groups
    }

class OriginalGwasArchiveRecheck(unittest.TestCase):
    def test_original_thirty_region_archive_pinned_and_chr10(self):
        z=m.source_table(OLD)
        self.assertEqual(z["region_id"],"IS_XDATA_G0015")
        self.assertEqual(z["source_dataset"],"GWS_EAS_AIS")
        self.assertEqual(int(z["n_rows_in_context"]),3913)
        self.assertEqual(int(z["n_genomewide_p_lt_5e8"]),69)

    def test_numbers_alone_not_causal(self):
        z=m.reconcile(fake("EUR"),fake("EAS"),m.source_table(OLD))
        self.assertTrue(z["descriptive_region_signature_match"])
        self.assertEqual(z["groups"],36)
        self.assertFalse(z["twenty_independent_nonpal_source_SNPs_matched"])
        self.assertFalse(z["sample_size_and_BBJ_study_source_approved"])
        self.assertEqual(z["new_P3A_coloc_executed"],0)

    def test_different_source_counts_flagged_not_silently_rewritten(self):
        eas=fake("EAS");eas["context"]["source_context_p_lt_5e8_rows"]=50
        z=m.reconcile(fake("EUR"),eas,m.source_table(OLD))
        self.assertFalse(z["regional_genomewide_significance_count_equal"])
        self.assertFalse(z["descriptive_region_signature_match"])
        self.assertEqual(z["scientific_status"],"HOLD")

    def test_gwas_source_allele_approval_not_granted_on_one_palindromic_lead(self):
        eas=fake("EAS")
        eas["context"]["legacy_GWS_EAS_AIS_lead_comparison"]["n_independent_original_variant_matches_for_source_approval"]=1
        with self.assertRaisesRegex(ValueError,"palindromic single lead"):
            m.reconcile(fake("EUR"),eas,m.source_table(OLD))

    def test_real_original_digest_tamper_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/"historic.tsv"
            p.write_bytes(OLD.read_bytes()+b"missing\n")
            with self.assertRaisesRegex(ValueError,"Git blob drift"):
                m.source_table(p)

    def test_36_locus_conflation_rejected(self):
        eas=fake("EAS")
        eas["actual_36_assay_groups_considered"]=35
        with self.assertRaisesRegex(ValueError,"scientific STOP"):
            m.reconcile(fake("EUR"),eas,m.source_table(OLD))

if __name__=="__main__":
    unittest.main()
