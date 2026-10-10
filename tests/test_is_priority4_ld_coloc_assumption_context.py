import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from build_is_priority4_ld_coloc_assumption_context import run,render

class Priority4AssumptionContextTest(unittest.TestCase):
 def test_combined_context_reproduces_archived_h4_without_ld_reestimate(self):
  base=Path("/srv/is-analysis/results/is/audits")
  gene=base/"is_priority4_ld_proxy_ancestry_stability_20261011_v1/IS_PRIORITY4_GENE_NOMINAL_EQTL_TAGGING.tsv"
  prior=base/"is_30_abf_joint_sensitivity_20261011_v1/IS_30_ABF_JOINT_SOURCE_N_AND_PRIOR_SENSITIVITY.tsv"
  replay=base/"legacy_coloc_646_batch_cache30_20261010_v1/IS_LEGACY_646_SNP_REPLAY_STATUS.tsv"
  if not all(x.exists() for x in [gene,prior,replay]):
   self.skipTest("Research server audit sources not present")
  with tempfile.TemporaryDirectory() as t:
   o=Path(t);m=run(gene,prior,replay,o,o/"context.svg")
   self.assertEqual(m["source_genes"],4)
   self.assertEqual(m["distinct_loci"],2)
   self.assertFalse(m["current_analysis_recomputes_PP_H4_using_LD"])
   self.assertFalse(m["original_coloc_ancestry_LD_input_provided"])
   self.assertEqual(m["genes_top_QTL_lead_r2_ge_0p8_in_both"],1)
   self.assertEqual(m["original_candidate_gene_universe"],2225)
   g={x["symbol"]:x for x in m["rows"]}
   self.assertAlmostEqual(g["CALHM2"]["archived_PP_H4"],0.79768355258068)
   self.assertEqual(g["CALHM2"]["n_nominal_QTL_p_lt_1e4"],"0")
   self.assertEqual(g["FGF5"]["n_nominal_QTL_p_lt_1e4"],"19")
   self.assertTrue((o/"context.svg").exists())
   root=ET.parse(o/"context.svg").getroot()
   content=" ".join((x.text or "") for x in root.iter())
   self.assertIn("NOT formal multi-signal colocalization",content)
if __name__=="__main__":unittest.main()
