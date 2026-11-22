"""Synthetic fixtures and archived 36-group ID tests. No GWAS network fetch in unit tests."""
from __future__ import annotations
import copy
import csv
import gzip
import importlib
import json
import tempfile
from pathlib import Path
import unittest

s=importlib.import_module("scripts.is.audit_p3a_chr10_original_publisher_gwas")
j=importlib.import_module("scripts.is.join_p3a_chr10_publisher_source")
MAN=Path("docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json")

def row(pos=105596641,beta="-0.0762",p="6.214e-11",se="0.0115",ea="C",oa="G"):
    return {"chromosome":"10","base_pair_location":str(pos),"effect_allele":ea,
            "other_allele":oa,"effect_allele_frequency":"0.64",
            "beta":beta,"standard_error":se,"p_value":p}

def example_source(ancestry):
    x=s.checked_manifest(MAN)
    ctx={"scientific_status":"HOLD",
         "source_GWAS_effect_and_SE_original_rows_checked":True,
         "source_GWAS_effects_reconciled_to_archived_fullcis_36_groups":False,
         "source_context_GRCh37":dict(s.REGION),
         "source_context_valid_variant_rows":2000,
         "original_region_valid_variant_rows":900,
         "legacy_GWS_EAS_AIS_lead_comparison":None}
    rows=s.build_group_rows(x,ctx,ancestry)
    return {"schema":"is-p3a-chr10-original-publisher-gwas-all-36-v1",
            "technical_status":"PASS_PUBLISHER_WHOLE_FILE_AND_CHR10_ORIGINAL_ROWS",
            "scientific_status":"HOLD","ancestry":ancestry,
            "source":{"publisher_md5_and_prior_sha256_verified":True,
                      "accession":s.GATE.SOURCES[ancestry]["accession"]},
            "context":ctx,"actual_36_assay_groups_considered":36,
            "unique_candidate_gene_ids":20,
            "new_coloc_performed":0,"independent_GWAS_genomic_loci_replicated":0,
            "all_group_input_statuses":rows}

class ThirtySixPanelRealPublisherSourceTests(unittest.TestCase):
    def test_all_36_authentic_qtl_keys_loaded(self):
        groups=s.checked_manifest(MAN)
        self.assertEqual(len(groups),36)
        self.assertEqual(len({g["gene_id"] for g in groups}),20)
        self.assertEqual(len({g["region_id"] for g in groups}),1)

    def test_actual_original_effect_se_and_frequency(self):
        val=s.inspect_row(row())
        self.assertEqual(val["pos"],105596641)
        self.assertEqual(val["effect_allele"],"C")
        self.assertGreater(val["se"],0)
        self.assertTrue(val["palindrome"])

    def test_same_effect_allele_lead_match_but_source_unapproved(self):
        comp=s.compare_prior_lead(s.inspect_row(row()))
        self.assertEqual(comp["status"],"PALINDROMIC_LEAD_NUMERIC_CONCORDANCE_STRAND_HOLD")
        self.assertFalse(comp["source_identity_approved"])
        self.assertEqual(comp["n_independent_original_variant_matches_for_source_approval"],0)
        self.assertTrue(comp["strand_palindromic_unresolved"])

    def test_reversed_effect_allele_sign_is_harmonized_only_as_test(self):
        comp=s.compare_prior_lead(s.inspect_row(row(beta="0.0762",ea="G",oa="C")))
        self.assertTrue(comp["sign_and_p_concordant"])

    def test_opposite_allele_direction_without_sign_flip_is_detected(self):
        comp=s.compare_prior_lead(s.inspect_row(row(beta="0.0762",ea="C",oa="G")))
        self.assertEqual(comp["status"],"PALINDROMIC_LEAD_NUMERIC_DISCORDANCE_STRAND_HOLD")

    def test_wrong_alleles_rejected_as_matching_lead(self):
        comp=s.compare_prior_lead(s.inspect_row(row(ea="A",oa="T")))
        self.assertEqual(comp["status"],"ALLELES_DISCORDANT")

    def test_unusable_source_se_refused(self):
        self.assertIn("error",s.inspect_row(row(se="0")))
        self.assertIn("error",s.inspect_row(row(se="-2")))
        self.assertIn("error",s.inspect_row(row(se="NA")))

    def test_unusable_source_p_eaf_refused(self):
        self.assertIn("error",s.inspect_row(row(p="1.1")))
        e=row();e["effect_allele_frequency"]="-0.1"
        self.assertIn("error",s.inspect_row(e))

    def test_no_region_guess_outside_chr10(self):
        self.assertIsNone(s.inspect_row(row(pos=45000000)))
        e=row();e["chromosome"]="11"
        self.assertIsNone(s.inspect_row(e))

    def test_original_gzip_fixture_gwas_se_scanned_and_not_promoted(self):
        rows=[row(),row(pos=105500001,beta="0.1",p="0.02"),row(pos=105600000,beta="0.2",p="0.1"),
              row(pos=99999999),dict(row(pos=105600010),standard_error="NA")]
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/"fixture.tsv.gz"
            with gzip.open(path,"wt") as f:
                writer=csv.DictWriter(f,fieldnames=s.HEADER,delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            parsed=s.scan_original(path,"EAS",min_total_rows=5,min_context_rows=3)
        self.assertEqual(parsed["source_context_valid_variant_rows"],3)
        self.assertEqual(parsed["source_context_failed_numeric_or_allele_rows"],1)
        self.assertIsNotNone(parsed["legacy_GWS_EAS_AIS_lead_comparison"])
        self.assertFalse(parsed["source_GWAS_effects_reconciled_to_archived_fullcis_36_groups"])

    def test_gzip_missing_se_header_denied(self):
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/"bad.tsv.gz"
            with gzip.open(path,"wt") as f:f.write("chromosome\tbase_pair_location\tp_value\n")
            with self.assertRaisesRegex(ValueError,"missing beta/SE"):
                s.scan_original(path,"EAS",min_total_rows=1,min_context_rows=1)

    def test_no_group_level_evidence_is_inferred(self):
        source=example_source("EAS")
        groups=source["all_group_input_statuses"]
        self.assertEqual(len(groups),36)
        self.assertTrue(all(not x["variant_level_QTL_GWAS_signed_effect_match_proven"] for x in groups))
        self.assertTrue(all(not x["GWAS_exact_case_fraction_and_group_N_approved"] for x in groups))

    def test_merge_two_original_sources_is_single_locus(self):
        out=j.audit(example_source("EUR"),example_source("EAS"))
        self.assertEqual(out["real_group_keys"],36)
        self.assertEqual(out["real_gene_ids"],20)
        self.assertEqual(out["distinct_analysed_GWAS_region_count"],1)
        self.assertEqual(out["coloc_abf_executed"],0)
        self.assertEqual(len(out["per_group_status"]),36)

    def test_merge_fails_if_fake_group_approved(self):
        eur=example_source("EUR");eas=example_source("EAS")
        eas["all_group_input_statuses"][0]["GWAS_exact_case_fraction_and_group_N_approved"]=True
        with self.assertRaisesRegex(ValueError,"incorrectly promoted"):
            j.audit(eur,eas)

    def test_merge_fails_if_gwas_source_not_publisher_verified(self):
        eur=example_source("EUR");eas=example_source("EAS")
        eur["source"]["publisher_md5_and_prior_sha256_verified"]=False
        with self.assertRaisesRegex(ValueError,"publisher source"):
            j.audit(eur,eas)

    def test_merge_fails_if_fake_second_locus(self):
        eur=example_source("EUR");eas=example_source("EAS")
        eas["context"]["source_context_GRCh37"]["original_region_start"]=104663677
        with self.assertRaisesRegex(ValueError,"region source span"):
            j.audit(eur,eas)

    def test_merge_fails_if_original_p3a_group_missing(self):
        eur=example_source("EUR");eas=example_source("EAS")
        eas["all_group_input_statuses"].pop()
        with self.assertRaisesRegex(ValueError,"scientific eligibility"):
            j.audit(eur,eas)

if __name__=="__main__":
    unittest.main()
