"""COL4A1/COL4A2 matched biological SNP and public cis-source invariants."""
import importlib.util,sys,unittest
from pathlib import Path
from collections import Counter
p=Path(__file__).resolve().parents[1]/"scripts/is/prepare_is_col4a1_col4a2_paired_cis_coloc.py"
spec=importlib.util.spec_from_file_location("col4_pair",p)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m
spec.loader.exec_module(m)
class TestCOL4Pair(unittest.TestCase):
 def test_all_tissue_N_from_qtl_original_AN(self):
  self.assertEqual(m.DATASETS,{"QTD000131":387,"QTD000171":205,"QTD000434":560})
  self.assertEqual(m.SOURCE_COUNTS[("COL4A1","QTD000434")],8763)
  self.assertEqual(m.SOURCE_COUNTS[("COL4A2","QTD000434")],8764)
 def test_pair_matched_snp_universe_identical_by_dataset(self):
  qtl={};bbj={};ais={}
  for did in m.DATASETS:
   for loc in range(2):
    snp=f"13:{110954237+loc}:A:G"
    b={"beta_alt":.08,"se":.01,"p":2e-9,"af_alt":.30}
    a={"beta_alt":.07,"se":.02,"p":1e-8,"af_alt":.31}
    bbj[snp]=b;ais[snp]=a
    for gn in m.GENES:
     qtl[(gn,did,snp)]={
      "QTL_p":1e-5 if gn=="COL4A2" else .5,
      "QTL_alt_AF":.30,"QTL_beta_ALT":-.2,"QTL_se":.04,
      "QTL_n_variant_AN_div_2":m.DATASETS[did],
      "source_GRCh38":"chr13_110301890_A_G"}
  data,summ=m.compare(qtl,set(bbj),bbj,ais,[])
  self.assertEqual(len(data),12)
  self.assertEqual(len(summ),6)
  self.assertTrue(all(int(q["paired_nonpal_SNP_count"])==2 for q in summ))
  self.assertEqual(sum(q["both_marginal_sameSNP_QTL1e4_AIS5e8"] for q in summ if q["gene_symbol"]=="COL4A2"),6)
  self.assertEqual(sum(q["both_marginal_sameSNP_QTL1e4_AIS5e8"] for q in summ if q["gene_symbol"]=="COL4A1"),0)
  self.assertTrue(all(r["causal_gene_label"]=="NOT_ESTABLISHED" for r in data))
 def test_do_not_include_palindromic_snp(self):
  q={};b={};a={}
  for did in m.DATASETS:
   for gn in m.GENES:
    v="13:10000000:A:T"
    q[(gn,did,v)]={
     "QTL_p":1e-9,"QTL_alt_AF":.30,"QTL_beta_ALT":.2,"QTL_se":.01,
     "QTL_n_variant_AN_div_2":m.DATASETS[did],"source_GRCh38":"chr13_10000000_A_T"}
   b[v]={"beta_alt":1,"se":1,"p":1e-9,"af_alt":.3}
   a[v]=dict(b[v])
  z,_=m.compare(q,set(b),b,a,[])
  self.assertEqual(len(z),0)
if __name__=="__main__":unittest.main()
