"""Offline tests: 15 public nominal slices, 124 new sourced but not coloc genes."""
import importlib.util,sys,tempfile,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/is/merge_is_eas30_public_nominal_qtl_15_sources.py"
s=importlib.util.spec_from_file_location("merge_eas_15",p)
m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m)
def buildfixture(td):
 base=Path(td);technical=base/"technical";t1=base/"t1"
 technical.mkdir();t1.mkdir()
 genes=[];slices=[];sources=[];pr=[]
 total=0
 for region,(n,covered) in m.EXPECT.items():
  for i in range(n):
   gene=f"ENSG{total:011d}";total+=1
   if region=="IS_XDATA_G0015":old=(i<16 or i==70)
   elif region=="IS_XDATA_G0023":old=i<6
   elif region=="IS_XDATA_G0007":old=(i<3 or i==9)
   else:old=False
   genes.append({"gene_id_stable":gene,"canonical_GENCODE_v19_gene_symbol":f"SYMBOL{total}",
                 "region_ids":region,"new_DIRECT_prior_grid_assays":"1" if old else "0"})
   for did in sorted(m.DATASETS):
    slices.append({"region_id":region,"dataset_id":did,"gene_id":gene,
                   "n_raw_QTL_records":"10" if i<covered else "0",
                   "min_QTL_p":"0.002" if i<covered else ""})
  for did in sorted(m.DATASETS):
   root=t1 if region in ("IS_XDATA_G0015","IS_XDATA_G0023","IS_XDATA_G0007") else technical
   path=root/f"{region}__{did}__GRCh38_1_100_200.tsv.gz"
   path.write_bytes(b"tiny sample source sha")
   sources.append({"region_id":region,"dataset_id":did,
     "status":"PUBLIC_NOMINAL_QTL_REGION_DOWNLOADED",
     "regional_qtl_rows":"1000","target_gene_qtl_rows":"100",
     "queried_region_GRCh38":"1:100-200","raw_sha256":m.sha(path)})
 while total<869:
  i=total;total+=1
  region=f"IS_XDATA_G{(i%25)+1:04d}"
  if region in m.EXPECT:region="IS_XDATA_G0019"
  genes.append({"gene_id_stable":f"ENSG{i:011d}","canonical_GENCODE_v19_gene_symbol":f"SYMBOL{i}",
   "region_ids":region,"new_DIRECT_prior_grid_assays":"0"})
 while total<2225:
  i=total;total+=1
  genes.append({"gene_id_stable":f"ENSG{i:011d}","canonical_GENCODE_v19_gene_symbol":f"SYMBOL{i}",
   "region_ids":"IS_EUR_AIS_G0010","new_DIRECT_prior_grid_assays":"0"})
 for i in range(1,31):
  rid=f"IS_XDATA_G{i:04d}"
  label="T1" if rid in ("IS_XDATA_G0015","IS_XDATA_G0023","IS_XDATA_G0007","IS_XDATA_G0022") else "T4"
  pr.append({"region_id":rid,"primary_GWAS_EAS_BBJ_AIS_tier":label})
 return sources,slices,genes,pr,{"T1":t1,"T4":technical}
class TestSourceMerger(unittest.TestCase):
 def test_all_15_local_sha_and_124_new_gene_readiness(self):
  with tempfile.TemporaryDirectory() as td:
   src,ql,genes,pr,roots=buildfixture(td)
   gr,reg,s=m.build(src,ql,genes,pr,roots)
   self.assertEqual(len(gr),2225)
   self.assertEqual(s["n_new_nominal_QTL_union_gene_ids"],149)
   self.assertEqual(s["n_nominal_QTL_union_gene_ids_without_historic_ABF"],124)
   self.assertEqual(s["n_T1_nominal_QTL_genes_without_historic_ABF"],34)
   self.assertTrue(all(not x["GWAS_and_QTL_alleles_coharmonized"] for x in gr))
 def test_changed_source_sha_fails(self):
  with tempfile.TemporaryDirectory() as td:
   src,ql,genes,pr,roots=buildfixture(td)
   src[0]["raw_sha256"]="INVALID"
   with self.assertRaises(ValueError):m.build(src,ql,genes,pr,roots)
 def test_missing_nominal_source_row_fails(self):
  with tempfile.TemporaryDirectory() as td:
   src,ql,genes,pr,roots=buildfixture(td)
   with self.assertRaises(ValueError):m.build(src,ql[:-1],genes,pr,roots)
if __name__=="__main__":unittest.main()
