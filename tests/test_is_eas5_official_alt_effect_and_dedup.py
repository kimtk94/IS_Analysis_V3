"""Official eQTL Catalogue ALT beta, AC/AN frequency and duplicate-rsid regression."""
import copy,csv,gzip,importlib.util,sys,tempfile,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/is/audit_is_eas5_official_alt_effect_and_dedup.py"
s=importlib.util.spec_from_file_location("official_qtl_alt",p)
m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m)
def fixture():
 sources=[{"region_id":f"IS_XDATA_G{i:04d}","dataset_id":f"QTD{i+130:06d}"} for i in range(15)]
 sources[0]={"region_id":"IS_XDATA_G0015","dataset_id":"QTD000131"}
 genes=[{"gene_id_stable":f"ENSG{i:011d}",
         "canonical_gene_symbol":"STUDYGENE" if i==0 else f"GENE{i}",
         "archived_646_coloc_assays":0} for i in range(2225)]
 row={"region_id":"IS_XDATA_G0015","QTL_dataset":"QTD000131",
      "gene_id":"ENSG00000000000","mapped_variant_GRCh37":"10:101:A:G",
      "source_variant_GRCh38":"chr10_110_A_G","mapped_strand":"+",
      "QTL_p":"0.000001","QTL_beta_SOURCE_ORIENTATION_UNVERIFIED":"0.4",
      "QTL_se":"0.1","QTL_maf":"0.25","QTL_ac":"50","QTL_an":"200",
      "BBJ_is_present":"True","BBJ_beta_ALT":"-0.2",
      "BBJ_p":"0.00001","BBJ_ALT_EAF":"0.30",
      "EAS_AIS_is_present":"True","EAS_AIS_beta_ALT":"-0.3",
      "EAS_AIS_p":"0.000001","EAS_AIS_ALT_EAF":"0.28",
      "palindromic_SNP_review":"False"}
 return sources,genes,row
def testsource(tmp,rows):
 p=Path(tmp)/"matched.tsv.gz"
 with gzip.open(p,"wt",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
  w.writeheader();w.writerows(rows)
 return p
class TestEffectConvention(unittest.TestCase):
 def test_official_alt_ac_an_and_duplicate_rsid(self):
  with tempfile.TemporaryDirectory() as td:
   src,genes,r=fixture()
   input=testsource(td,[r,r.copy()])
   out=Path(td)/"output";out.mkdir()
   report=m.process(input,src,genes,out)
   self.assertEqual(report["n_unique_harmonized_geneSNP_records"],1)
   self.assertEqual(report["duplicate_source_rsids_removed_per_gene_variant"],1)
   vals=m.rd(out/"IS_EAS5_OFFICIAL_QTL_ALT_GWAS_HARMONIZED_VARIANTS.tsv")
   self.assertEqual(len(vals),1)
   self.assertAlmostEqual(float(vals[0]["QTL_beta_mapped_GRCh37_ALT_OFFICIAL_DEFINITION"]),.4)
   self.assertAlmostEqual(float(vals[0]["QTL_alt_eaf_ac_div_an"]),.25)
   self.assertEqual(vals[0]["coloc_recomputed"],"False")
 def test_duplicate_with_inconsistent_beta_rejected(self):
  with tempfile.TemporaryDirectory() as td:
   src,genes,r=fixture()
   q=r.copy();q["QTL_beta_SOURCE_ORIENTATION_UNVERIFIED"]="-0.4"
   input=testsource(td,[r,q])
   out=Path(td)/"out";out.mkdir()
   with self.assertRaises(ValueError):m.process(input,src,genes,out)
 def test_false_ancestry_or_coloc_never_assumed(self):
  with tempfile.TemporaryDirectory() as td:
   src,genes,r=fixture()
   input=testsource(td,[r]);out=Path(td)/"out";out.mkdir()
   m.process(input,src,genes,out)
   vals=m.rd(out/"IS_EAS5_OFFICIAL_QTL_ALT_GWAS_HARMONIZED_VARIANTS.tsv")
   self.assertEqual(vals[0]["source_cohort_ancestry_verified"],"False")
   self.assertEqual(vals[0]["matched_QTL_LD_verified"],"False")
if __name__=="__main__":unittest.main()
