import csv
import json
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_priority4_ebi_bbj_abf_grid import analyze
from render_is_priority4_ebi_bbj_abf_prior_n_figure import render
from prepare_is_priority4_ebi_bbj_abf_pilot import FIELDS

ROOT=Path("/srv/is-analysis/results/is/audits")
class PriorityFourQTDColocPilotContract(unittest.TestCase):
 def sources(self):
  grid=ROOT/"is_priority4_ebi_bbj_abf_pilot_20261011_v3/IS_PRIORITY4_EBI_BBJ_ABF_60_SCENARIO_SENSITIVITY.tsv"
  qc=ROOT/"is_priority4_ebi_bbj_abf_pilot_20261011_v3/IS_PRIORITY4_EBI_BBJ_ABF_QTL_N_PROXY_QC.tsv"
  ready=ROOT/"is_priority4_ebi_bbj_abf_pilot_20261011_v2/IS_PRIORITY4_EBI_BBJ_ABF_QC_VALIDATED_INPUT.tsv"
  if not all(x.exists() for x in (grid,qc,ready)):
   self.skipTest("Immutable research-server reference source is unavailable")
  return grid,qc,ready
 def test_full_four_gene_grid_and_science_gates(self):
  g,q,i=self.sources()
  with tempfile.TemporaryDirectory() as t:
   r=analyze(g,q,i,Path(t))
   self.assertEqual(r["audited_grid"],60)
   self.assertTrue(r["qtl_effective_N_is_unverified"])
   self.assertFalse(r["original_historical_coloc_H4_modified"])
   self.assertFalse(r["coloc_susie_performed"])
   v={x["gene_symbol"]:x for x in r["summary"]}
   self.assertEqual({x:x2["matched_nonpal_SNPs"] for x,x2 in v.items()},
      {"FGF5":2773,"C4orf22":2772,"CALHM2":2766,"NEURL":2739})
   self.assertAlmostEqual(v["FGF5"]["H4_reference_N175_p12_1e5"],0.913169913672319,places=10)
   self.assertAlmostEqual(v["CALHM2"]["H4_reference_N175_p12_1e5"],0.905405114815332,places=10)
   self.assertAlmostEqual(v["NEURL"]["H4_reference_N175_p12_1e5"],0.0605188397,places=6)
   self.assertFalse(v["NEURL"]["H4_crosses_0p8_among_scenarios"])
   self.assertTrue(v["FGF5"]["H4_crosses_0p8_among_scenarios"])
 def test_preparation_preserves_original_refalt_and_unattested_qtl_n(self):
  _,_,ready=self.sources()
  with ready.open(newline="") as f:
   rd=csv.DictReader(f,delimiter="\t")
   self.assertEqual(list(rd.fieldnames),FIELDS)
   data=list(rd)
  self.assertEqual(len(data),11050)
  self.assertEqual(len({r["variant_GRCh37_REF_ALT"] for r in data}),5666)
  for r in data:
   self.assertEqual(r["QTD_nominal_full_cis_tested"],"True")
   self.assertEqual(r["BBJ_GWAS_exact_FASTA_ref"],"True")
   self.assertEqual(r["harmonization"],"ALT_AFTER_CHAIN_AND_FASTA_REF_EQUIVALENCE")
   self.assertGreater(float(r["BBJ_se"]),0)
   self.assertGreater(float(r["EBI_QTD_se"]),0)
   self.assertLessEqual(float(r["ebi_genotype_AN_over_2_not_QTL_effective_N"]),175)
 def test_figure_contains_sensitivity_boundary_labels(self):
  grid,_,_=self.sources()
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/"fig.svg"
   self.assertEqual(render(grid,out),60)
   root=ET.parse(out).getroot()
   texts=" ".join((e.text or "") for e in root.iter())
   self.assertIn("Exploratory single-signal ABF",texts)
   self.assertIn("FGF5",texts)
   self.assertIn("NEURL",texts)
if __name__=="__main__":unittest.main()
