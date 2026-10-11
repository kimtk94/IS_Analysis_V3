"""eQTL Catalogue n_variants full cis completeness and non-destructive source recovery."""
import csv,gzip,importlib.util,sys,tempfile,unittest
from pathlib import Path
root=Path(__file__).resolve().parents[1]/"scripts/is"
def load(file,name):
 spec=importlib.util.spec_from_file_location(name,root/file)
 m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
 return m
a=load("audit_is_eas5_qtl_full_cis_permuted_universe.py","fullcis5")
b=load("extend_is_g0015_public_qtl_complete_cis.py","ciswideG0015")
def row(gene,chrom,pos,ref,alt,rs):
 values=[gene,chrom,str(pos),ref,alt,f"chr{chrom}_{pos}_{ref}_{alt}",
   "100","0.20",".0001",".10",".02","SNP","30","100",
   "NA",gene,gene,"3",rs]
 assert len(values)==19
 return "\t".join(values)+"\n"
class TestCIS(unittest.TestCase):
 def test_complete_partial_and_not_detected_are_distinct(self):
  self.assertEqual(a.infer_regional_status({"n_variants":7000},2000),
   "REGIONAL_SLICE_INCOMPLETE_CIS_NOT_SAFE_FOR_FULL_LOCUS_COLOC")
  self.assertEqual(a.infer_regional_status({"n_variants":7000},7000),
   "FULL_CIS_EQUIVALENCE_VS_EBI_NVARIANTS_COUNT_ONLY")
  self.assertEqual(a.infer_regional_status({"n_variants":7000},0),
   "PERMUTED_GENE_TESTED_BUT_NO_NOMINAL_ROWS_IN_REGIONAL_SLICE")
  self.assertEqual(a.infer_regional_status(None,0),
   "NO_EBI_PERMUTED_GENE_RECORD_CANNOT_ASSESS")
 def test_official_permuted_must_have_unique_gene_records(self):
  with tempfile.TemporaryDirectory() as t:
   p=Path(t)/"permuted.tsv.gz"
   with gzip.open(p,"wt") as f:
    f.write("molecular_trait_object_id\tn_traits\tn_variants\tp_beta\tp_perm\n")
    f.write("ENSG123\t1\t1000\t0.003\t0.001\n")
   x,n=a.parse_permuted(p,{"ENSG123"})
   self.assertEqual(n,1)
   self.assertEqual(x["ENSG123"]["n_variants"],1000)
   with gzip.open(p,"at") as f:f.write("ENSG123\t1\t1000\t0.003\t0.001\n")
   with self.assertRaises(ValueError):a.parse_permuted(p,{"ENSG123"})
 def test_duplicate_rsid_not_distinct_cis_variants(self):
  with tempfile.TemporaryDirectory() as t:
   p=Path(t)/"qtl.tsv.gz"
   with gzip.open(p,"wt") as f:
    f.write(row("ENSG123","10",103000001,"A","G","rs1"))
    f.write(row("ENSG123","10",103000001,"A","G","rs2"))
   unq,n,dup,total,_=a.region_variants(p,{"ENSG123"},"10")
   self.assertEqual(n["ENSG123"],2)
   self.assertEqual(len(unq["ENSG123"]),1)
   self.assertEqual(dup["ENSG123"],1)
 def test_extended_cis_multiset_contains_original(self):
  with tempfile.TemporaryDirectory() as t:
   old=Path(t)/"old.tsv.gz";wide=Path(t)/"wide.tsv.gz"
   same=row("ENSG123","10",103000001,"A","G","rs1")
   with gzip.open(old,"wt") as f:f.write(same)
   with gzip.open(wide,"wt") as f:
    f.write(same);f.write(row("ENSG123","10",103000002,"C","T","rs2"))
   v,rec,dup,nn,ne=b.collect(old,wide,{"ENSG123"})
   self.assertEqual((nn,ne,len(v["ENSG123"])),(1,2,2))
   with gzip.open(wide,"wt") as f:f.write(row("ENSG123","10",103000002,"C","T","rs2"))
   with self.assertRaises(ValueError):b.collect(old,wide,{"ENSG123"})
if __name__=="__main__":unittest.main()
