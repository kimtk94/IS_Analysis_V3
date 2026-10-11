import csv
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_priority4_ebi_bbj_nested_abf import run
from render_is_priority4_nested_abf_control_figure import render
from prepare_is_priority4_ebi_bbj_nested_snps import EXPECTED

BASE=Path("/srv/is-analysis/results/is/audits")
class NestedSameSourceABFTest(unittest.TestCase):
 def sources(self):
  full=BASE/"is_priority4_ebi_bbj_abf_pilot_20261011_v3/IS_PRIORITY4_EBI_BBJ_ABF_60_SCENARIO_SENSITIVITY.tsv"
  nested=BASE/"is_priority4_ebi_bbj_nested_fixedsource_20261011_v2/IS_PRIORITY4_EBI_BBJ_ABF_60_SCENARIO_SENSITIVITY.tsv"
  exact=BASE/"is_priority4_ebi_bbj_nested_fixedsource_20261011_v1/IS_PRIORITY4_SAME_EBI_BBJ_ARCHIVED_SNP_NESTED_INPUT.tsv"
  if not all(x.exists() for x in (full,nested,exact)):
   self.skipTest("Frozen server science sources unavailable")
  return full,nested,exact
 def test_same_source_paired_60_controls_and_invariants(self):
  full,nested,exact=self.sources()
  with tempfile.TemporaryDirectory() as td:
   data=run(full,nested,Path(td))
   self.assertEqual(data["paired_scenarios"],60)
   self.assertTrue(data["study_sources_fixed_between_paired_runs"])
   self.assertTrue(data["full_minus_nested_H4_all_negative"])
   self.assertTrue(data["no_SuSiE_finemapping_or_causal_gene_promotion"])
   self.assertTrue(data["cannot_attribute_historical_vs_new_H4_to_variant_coverage_alone"])
   rows={x["gene_symbol"]:x for x in data["summary"]}
   self.assertEqual({g:v["n_nested"] for g,v in rows.items()},EXPECTED)
   self.assertAlmostEqual(rows["FGF5"]["baseline_PP_H4_full_minus_nested"],-0.012829623779341,places=10)
   self.assertAlmostEqual(rows["C4orf22"]["baseline_PP_H4_full_minus_nested"],-0.034703348841776,places=10)
   for x in rows.values():
    self.assertEqual(x["scenarios_full_H4_lower_than_nested"],15)
 def test_source_records_in_nested_equal_a_subset_of_full(self):
  full,nested,exact=self.sources()
  a=BASE/"is_priority4_ebi_bbj_abf_pilot_20261011_v2/IS_PRIORITY4_EBI_BBJ_ABF_QC_VALIDATED_INPUT.tsv"
  with a.open(newline="") as f:
   rr=csv.DictReader(f,delimiter="\t")
   old={(x["gene_symbol"],x["variant_GRCh37_REF_ALT"]):x for x in rr}
  with exact.open(newline="") as f:
   rr=csv.DictReader(f,delimiter="\t")
   subset=list(rr)
  self.assertEqual(len(subset),5256)
  for x in subset:
   key=(x["gene_symbol"],x["variant_GRCh37_REF_ALT"])
   self.assertIn(key,old)
   self.assertEqual(x,old[key])
 def test_figure_renders_without_mutating_statistics(self):
  full,nested,exact=self.sources()
  with tempfile.TemporaryDirectory() as td:
   new=Path(td);r=run(full,nested,new)
   fig=new/"f.svg"
   self.assertEqual(render(new/"IS_PRIORITY4_SAME_SOURCE_FULL_VS_NESTED_SUMMARY.tsv",fig),4)
   text=" ".join(e.text or "" for e in ET.parse(fig).getroot().iter())
   self.assertIn("no validated SuSiE",text)
   self.assertIn("FGF5",text)
if __name__=="__main__":unittest.main()
