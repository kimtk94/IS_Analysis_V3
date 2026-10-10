"""Regression: original 646 assays, priors, 2225 gene coverage, no causal promotion."""
import importlib.util,sys,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/is/summarize_is_646_prior_h3h4_2225_genes.py"
spec=importlib.util.spec_from_file_location("is_646_prior_gate",p)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m
spec.loader.exec_module(m)

def fixture():
 genes=[{
  "gene_id_stable":f"ENSG{i:011d}",
  "canonical_GENCODE_v19_gene_symbol":f"GENE{i}",
  "n_positional_regions":"1",
  "human_reference_feature_status":"NOT_TESTED"
 } for i in range(2225)]
 grid=[]
 for i in range(646):
  for pp in m.P12:
   baseline=pp==1e-5
   high=pp==1e-4
   h4=(0.9 if i<12 else 0.55) if high and i<565 else (
       .05 if high else .4 if baseline and i<40 else .05)
   h3=.3 if (high and i<565) or (baseline and i<40) else .6
   grid.append({
    "locus":m.LOCI[i%4],"dataset_key":f"GTEx_V8__Assay{i}",
    "gene_base":genes[i%43]["gene_id_stable"],
    "p1":"0.0001","p2":"0.0001","p12":str(pp),
    "PP_H0":"0.01","PP_H1":"0.05","PP_H2":"0.01",
    "PP_H3":str(h3),"PP_H4":str(h4),"n_snps":"200"
   })
 return grid,genes
class TestPriorSensitivity(unittest.TestCase):
 def test_original_count_and_preserved_universe(self):
  grid,genes=fixture()
  summary,regions,assays,outgenes=m.audit(grid,genes)
  self.assertEqual(len(summary),5)
  self.assertEqual(len(assays),646)
  self.assertEqual(len(outgenes),2225)
  self.assertEqual(sum(x["n_archived_coloc_gene_tissue_assays"]==0 for x in outgenes),2182)
  self.assertEqual(sum(x["original_H4_gt_H3"] for x in assays),40)
  self.assertEqual(sum(x["high_H4_gt_H3"] for x in assays),565)
  self.assertEqual(sum(x["high_H4_ge_0p8"] for x in assays),12)
  self.assertTrue(all(x["causal_status"]=="NOT_ESTABLISHED" for x in outgenes))
 def test_missing_prior_fails_closed(self):
  grid,genes=fixture()
  grid.pop()
  with self.assertRaises(ValueError):m.audit(grid,genes)
 def test_changed_allele_set_across_priors_fails_closed(self):
  grid,genes=fixture()
  grid[2]["n_snps"]="199"
  with self.assertRaises(ValueError):m.audit(grid,genes)
if __name__=="__main__":unittest.main()
