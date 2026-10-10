"""Offline regression: 30 GWAS regional categories and 164 T1 unscreened genes."""
import importlib.util,sys,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/is/build_is_eas30_gwas_qtl_expansion_queue.py"
s=importlib.util.spec_from_file_location("eas30queue",p)
m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m)
TRAITS=("BBJ_IS","GWS_EAS_AS","GWS_EAS_AIS","GWS_EAS_CES","GWS_EAS_LAS","GWS_EAS_SVS")
def fixture():
 gwas=[];regions=[]
 for i in range(1,31):
  rid=f"IS_XDATA_G{i:04d}"
  miss=41 if i<=4 else 37 if i==5 else 25
  tested=43 if i==1 else 0
  regions.append({
   "EAS_region_id":rid,"n_positional_genes":str(miss+tested),
   "n_genes_archived_GTEx_ABF_tested":str(tested),
   "n_genes_not_tested":str(miss)})
  for trait in TRAITS:
   p=.001
   if i<=4 and trait in ("BBJ_IS","GWS_EAS_AIS"):p=1e-9
   if i==5 and trait=="GWS_EAS_CES":p=1e-9
   if i in (6,7) and trait=="GWS_EAS_AS":p=1e-9
   gwas.append({
    "region_id":rid,"source_dataset":trait,"best_p":str(p),
    "group_source_lead_present":"True"
   })
 return gwas,regions,[{} for i in range(2225)]
class TestEAS30Queue(unittest.TestCase):
 def test_four_t1_genetic_regions_not_causality(self):
  rows,counts=m.build(*fixture())
  self.assertEqual(len(rows),30)
  self.assertEqual(dict(counts),{"T1":4,"T2":1,"T3":2,"T4":23})
  self.assertEqual(sum(r["n_EAS_positional_genes_without_archived_coloc"] for r in rows),826)
  self.assertEqual(sum(r["n_EAS_positional_genes_without_archived_coloc"] for r in rows if r["primary_GWAS_EAS_BBJ_AIS_tier"].startswith("T1")),164)
  self.assertTrue(all(not r["cross_source_genetic_replication_tested"] for r in rows))
 def test_incomplete_outcome_fails_closed(self):
  gw,regs,genes=fixture()
  gw.pop()
  with self.assertRaises(ValueError):m.build(gw,regs,genes)
 def test_public_qtl_is_not_coloc(self):
  gw,regs,genes=fixture()
  pilot=[{"region_id":"IS_XDATA_G0025","dataset_id":f"QTD{i:06d}",
          "status":"PUBLIC_NOMINAL_QTL_REGION_DOWNLOADED"} for i in range(6)]
  rows,_=m.build(gw,regs,genes,pilot)
  x=next(r for r in rows if r["region_id"]=="IS_XDATA_G0025")
  self.assertTrue(x["new_public_QTL_nominal_row_coverage_not_coloc"])
  self.assertFalse(x["new_cohort_matched_QTL_LD_verified"])
if __name__=="__main__":unittest.main()
