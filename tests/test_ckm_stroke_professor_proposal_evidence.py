"""Prevent unverified scientific claims or model/count drift in IS proposal."""
import json
import pathlib
import re
import unittest

ROOT=pathlib.Path(__file__).resolve().parents[1]
FILE=ROOT/"config/ckm_stroke/proposal_evidence_20261122.json"
class ProposalEvidence(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.d=json.loads(FILE.read_text(encoding="utf-8"))
  cls.doc=(ROOT/cls.d["proposal_path"]).read_text(encoding="utf-8")
 def test_evidence_schema_and_proposal_exist(self):
  self.assertEqual(self.d["schema"],"CKM_IS_professor_proposal_gate_v1")
  self.assertIn("# 교수님 연구제안서",self.doc)
  self.assertFalse(self.d["replaces_prior_IS_master"])
 def test_proposal_only_no_science_results(self):
  self.assertEqual(self.d["proposal_status"],"GO_FOR_SUPERVISOR_DISCUSSION_ONLY")
  self.assertEqual(self.d["scientific_status"],"HOLD")
  for key,val in self.d["scientific_gates"].items():
   if key=="original_eas_source_alignment":
    self.assertEqual(val,"PASS_TECHNICAL_ONLY")
   else:
    self.assertIn(val,("HOLD","NOT_RUN"))
 def test_core_stroke_accession_has_ce_laa_no_swap(self):
  self.assertEqual(self.d["eur_outcomes"],{
   "AIS":"GCST90104540","LAA":"GCST90104542",
   "CE":"GCST90104541","SVD":"GCST90104543"})
  for x in self.d["eur_outcomes"].values():self.assertIn(x,self.doc)
 def test_ancestry_and_models_remain_separate(self):
  e=self.d["eas_context"]
  self.assertFalse(e["phenotype_equivalence_proven"])
  self.assertFalse(e["participant_nonoverlap_proven"])
  self.assertEqual((e["TPMI_SAIGE_cases"],e["TPMI_SAIGE_controls"]),(9249,304660))
  self.assertEqual((e["EAS_AIS_discovery_cases"],e["EAS_AIS_discovery_controls"]),(19032,237242))
  self.assertTrue(self.d["methodology_rules"]["no_cross_model_PLINK2_SAIGE_N_transfer"])
 def test_full_qc_evidence_arithmetic(self):
  e=self.d["eas_context"]
  self.assertEqual(e["SAIGE_original_rows"],15588353)
  self.assertEqual(e["SAIGE_source_QC_PASS"],8042797)
  self.assertGreater(e["SAIGE_original_rows"],e["SAIGE_source_QC_PASS"])
  self.assertGreater(e["SAIGE_source_QC_PASS"],e["SAIGE_EAS_allele_matched"])
  self.assertGreater(e["SAIGE_EAS_allele_matched"],e["SAIGE_EAS_strict_joint_INFO_0_8_MAF_0_01_AF_gap_0_05"])
  self.assertEqual(e["SAIGE_EAS_strict_joint_INFO_0_8_MAF_0_01_AF_gap_0_05"],3745225)
  for n in ("15,588,353","8,042,797","5,228,340","3,745,225"):self.assertIn(n,self.doc)
 def test_original_ci_evidence_links_present(self):
  e=self.d["eas_context"]
  for s in (e["action_source"],e["action_audit"]):self.assertIn(s,self.doc)
  self.assertTrue(e["independently_audited"])
 def test_prespecified_formal_subtype_difference_not_naive_sig(self):
  self.assertEqual(self.d["primary_exposures"],["LDL-C","SBP"])
  self.assertEqual(self.d["secondary_exposures"],["eGFR","UACR"])
  self.assertIn("공분산",self.doc)
  self.assertTrue(self.d["methodology_rules"]["no_subtype_effect_difference_from_significance_discordance"])
 def test_prior_art_and_proposal_research_logic(self):
  self.assertGreaterEqual(len(self.d["prior_art"]),5)
  self.assertIn("colocalization",self.doc)
  self.assertIn("NOT RUN",self.doc)
  self.assertTrue(self.d["methodology_rules"]["no_unverified_coloc_MR_LDSC_findings"])
 def test_scientific_hold_explicit(self):
  self.assertIn("SCIENTIFIC_RESULTS_HOLD",self.doc)
  self.assertIn("GO",self.doc)
  self.assertIn("HOLD",self.doc)
 def test_specimen_contact_deferred(self):
  self.assertTrue(self.d["methodology_rules"]["sample_or_specimen_contact_deferred"])
  self.assertIn("검체",self.doc)

if __name__=="__main__":unittest.main()
