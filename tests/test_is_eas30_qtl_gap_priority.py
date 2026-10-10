"""Regression: 2225 genes, 869 EAS genes/30 regions, 826 molecular-QTL gaps."""
import importlib.util,sys,unittest
from pathlib import Path
src=Path(__file__).resolve().parents[1]/"scripts/is/build_is_eas30_qtl_gap_priority.py"
spec=importlib.util.spec_from_file_location("is_scope_eas30",src)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m
spec.loader.exec_module(m)
def fixture():
 rows=[]
 for i in range(2225):
  eas=i<869
  rid=(i%4 if i<43 else i%30)
  tested=i<43
  n=604 if i==42 else 1 if tested else 0
  rows.append({
    "gene_id_stable":f"ENSG{i:011d}",
    "canonical_GENCODE_v19_gene_symbol":f"TEST{i}",
    "n_EAS_Japanese_regions":"1" if eas else "0",
    "n_EUR_regions":"1" if i>=669 else "0",
    "new_DIRECT_prior_grid_assays":str(n),
    "region_ids":f"IS_XDATA_G{rid+1:04d}" if eas else "IS_EUR_AIS_G0001"
  })
 return rows
class TestEAS30Scope(unittest.TestCase):
 def test_full_coverage_stays_2225_and_826_untested(self):
  genes,regions,groups=m.audit(fixture())
  self.assertEqual(len(genes),2225)
  self.assertEqual(len(regions),30)
  self.assertEqual(groups["EAS_REGION_QTL_NOT_TESTED"],826)
  self.assertEqual(groups["NON_EAS_REGION_OTHER_PIPELINE"],1356)
  self.assertEqual(sum(r["n_archived_gene_tissue_coloc_assays"] for r in regions),646)
  self.assertEqual(sum(r["n_genes_archived_GTEx_ABF_tested"]>0 for r in regions),4)
  self.assertTrue(all(not r["EAS_missing_QTL_is_negative_evidence"] for r in genes))
 def test_mismatched_eas_region_fails_closed(self):
  inp=fixture();inp[868]["region_ids"]="IS_EUR_AIS_G0002"
  with self.assertRaises(ValueError):m.audit(inp)
 def test_no_false_molecular_screen_outside_eas(self):
  inp=fixture();inp[-1]["new_DIRECT_prior_grid_assays"]="1"
  with self.assertRaises(ValueError):m.audit(inp)
if __name__=="__main__":unittest.main()
