"""Negative tests for real old-branch P3A 36-group identity recovery."""
from __future__ import annotations
import copy
import importlib
import json
from pathlib import Path
import unittest

m=importlib.import_module("scripts.is.audit_recovered_p3a_36_groups")
ROOT=Path("legacy/docs/data/is")
MANIFEST=Path("docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json")

def real():
    return (
        json.loads(MANIFEST.read_text(encoding="utf-8")),
        m.get_source(ROOT,"readiness"),
        m.get_source(ROOT,"egene"),
        m.get_source(ROOT,"samesnp"),
    )

class Restored36ActualSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=real()

    def test_real_original_source_is_36_assays_20_genes(self):
        v=m.validate(*copy.deepcopy(self.data))
        self.assertEqual(v["restored_group_key_count"],36)
        self.assertEqual(v["distinct_gene_id_count"],20)
        self.assertEqual(len(v["groups"]),36)
        self.assertEqual(sum(x["source_BH05_eGene_assay"]=="YES" for x in v["groups"]),25)
        self.assertEqual(v["by_QTL_dataset"],{"QTD000131":10,"QTD000171":11,"QTD000434":15})
        self.assertEqual(v["scientific_status"],"HOLD")
        self.assertEqual(len(v["gene_summary"]),20)
        self.assertEqual(sum(x["assay_groups"] for x in v["gene_summary"]),36)
        self.assertEqual(sum(x["public_BH05_eGene_assays"] for x in v["gene_summary"]),25)
        self.assertTrue(all(not x["new_coloc_approved"] for x in v["gene_summary"]))

    def test_no_new_coloc_or_aligned_ld_claim(self):
        v=m.validate(*copy.deepcopy(self.data))
        self.assertEqual(v["new_P3A_coloc_computations"],0)
        self.assertFalse(v["original_server_per_group_source_stats_positive_membership_restored"])
        self.assertTrue(all(x["QTL_matched_LD_verified"] is False for x in v["groups"]))

    def test_p3a_group_cannot_be_synthetic(self):
        v=list(copy.deepcopy(self.data))
        v[0]["groups"][0]["gene_id"]="ENSG12345678900"
        with self.assertRaisesRegex(ValueError,"identifiers"):
            m.validate(*v)

    def test_duplicate_original_gene_tissue_assay_rejected(self):
        v=list(copy.deepcopy(self.data))
        v[0]["groups"][0]=copy.deepcopy(v[0]["groups"][1])
        with self.assertRaisesRegex(ValueError,"missing/duplicated"):
            m.validate(*v)

    def test_bh05_studywide_q_is_read_from_original(self):
        v=list(copy.deepcopy(self.data))
        i=next(i for i,x in enumerate(v[0]["groups"])
               if x["official_studywide_eGene_BH05_assay"])
        v[0]["groups"][i]["official_QTL_studywide_BH_q"]=0.15
        with self.assertRaisesRegex(ValueError,"BH05"):
            m.validate(*v)

    def test_legacy_25_source_stat_membership_not_assumed(self):
        v=list(copy.deepcopy(self.data))
        v[0]["groups"][0]["original_P3A_source_stat_positive_proven"]=True
        with self.assertRaisesRegex(ValueError,"scientific stop"):
            m.validate(*v)

    def test_old_source_commit_pin_enforced(self):
        v=list(copy.deepcopy(self.data))
        v[0]["origin"]["commit_sha"]="0"*40
        with self.assertRaisesRegex(ValueError,"provenance"):
            m.validate(*v)

    def test_source_SNP_count_tamper_rejected(self):
        v=list(copy.deepcopy(self.data))
        v[0]["groups"][0]["archived_nonpalindromic_BBJ_EAS_SNP_overlap"]+=1
        with self.assertRaisesRegex(ValueError,"SNP overlap"):
            m.validate(*v)

    def test_sample_ld_not_approved_by_bh05(self):
        v=list(copy.deepcopy(self.data))
        v[0]["matched_LD_approved"]=True
        with self.assertRaisesRegex(ValueError,"biological"):
            m.validate(*v)

    def test_same_snp_source_no_causality(self):
        v=list(copy.deepcopy(self.data))
        match=next(x for x in v[3] if m.key(x)==m.key(v[0]["groups"][0]))
        match["colocalization_or_causal_gene_proven"]="True"
        with self.assertRaisesRegex(ValueError,"scientific stop"):
            m.validate(*v)

    def test_wrong_public_qtl_flag_rejected(self):
        v=list(copy.deepcopy(self.data))
        v[0]["groups"][0]["official_studywide_eGene_BH05_assay"]=False
        with self.assertRaisesRegex(ValueError,"membership"):
            m.validate(*v)

    def test_wrong_three_QTL_distribution_rejected(self):
        v=list(copy.deepcopy(self.data))
        match=next(x for x in v[1] if x["possible_single_signal_ABF_followup"]=="True")
        match["QTL_dataset"]="QTD000999"
        with self.assertRaisesRegex(ValueError,"distribution"):
            m.validate(*v)

if __name__=="__main__":
    unittest.main()
