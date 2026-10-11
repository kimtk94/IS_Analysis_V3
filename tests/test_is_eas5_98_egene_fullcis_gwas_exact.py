"""Guard downstream full-cis SNP-GWAS data gates without external downloads."""
import importlib.util,sys,unittest
from collections import Counter
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/is/audit_is_eas5_98_egene_fullcis_gwas_exact.py"
spec=importlib.util.spec_from_file_location("fullcis_bh05_gwas",p)
m=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=m;spec.loader.exec_module(m)
class Fake:
 qtl_g0015=Path("/a");qtl_t1=Path("/b");qtl_t4=Path("/c");qtl_g0018=Path("/d")
class TestFullcisGwas(unittest.TestCase):
 def test_GIGA_other_effect_historical_ID_is_not_REF_ALT(self):
  from pathlib import Path
  import sys
  sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
  import audit_is_eas5_public_qtl_gwas_snp_overlap as mod
  class Ref:
   def get(self,chrom,pos):
    assert (chrom,pos)==("5",29439275)
    return "T"
  source={"dataset":"GCST90104545","phenotype":"AIS","ancestry":"EAS",
    "build":"GRCh37","chr":"5","pos":"29439275",
    "effect_allele":"T","other_allele":"C",
    "variant_id":"5:29439275:C:T","beta":"0.0106","se":"0.0136",
    "p":"0.4351","eaf":"0.7798"}
  qc,r=m.checked_gigastroke_fasta_anchored(source,Ref(),mod)
  self.assertEqual(qc,"FASTA_ANCHORED_HISTORICAL_OTHER_EFFECT_ORDER")
  self.assertEqual(r["variant"],"5:29439275:T:C")
  self.assertAlmostEqual(r["beta_alt"],-.0106)
  self.assertAlmostEqual(r["af_alt"],.2202)
  self.assertFalse(r["palindrome"])
 def test_GIGA_bad_historical_ID_does_not_get_rescued(self):
  import sys
  sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
  import audit_is_eas5_public_qtl_gwas_snp_overlap as mod
  class Ref:
   def get(self,chrom,pos):return "T"
  r={"dataset":"GCST90104545","phenotype":"AIS","ancestry":"EAS",
    "build":"GRCh37","chr":"5","pos":"29439275",
    "effect_allele":"T","other_allele":"C",
    "variant_id":"5:29439275:A:T","beta":"0.0106","se":"0.0136",
    "p":"0.4351","eaf":"0.7798"}
  qc,data=m.checked_gigastroke_fasta_anchored(r,Ref(),mod)
  self.assertEqual(qc,"GIGA_FASTA_ALLELE_PAIR_OR_SOURCE_ID_MISMATCH")
  self.assertIsNone(data)
 def test_all_five_original_regions_to_verified_input_sources(self):
  self.assertEqual(m.file_for("IS_XDATA_G0015","QTD000131",Fake()),Path("/a/IS_XDATA_G0015__QTD000131__GRCh38_10_101403919_105460593.tsv.gz"))
  self.assertEqual(m.file_for("IS_XDATA_G0018","QTD000434",Fake()),Path("/d/IS_XDATA_G0018__QTD000434__GRCh38_11_377434_3387434.tsv.gz"))
  self.assertEqual(m.file_for("IS_XDATA_G0025","QTD000171",Fake()),Path("/c/IS_XDATA_G0025__QTD000171__GRCh38_17_5691879_8725409.tsv.gz"))
 def test_two_GRCh38_variants_colliding_on_same_GRCh37_pair_must_be_quarantined(self):
  from collections import Counter
  assay=("IS_XDATA_G0018","QTD000131","ENSG1")
  key=(assay,"11:1213481:G:A")
  records={};conflicts=set();stats=Counter()
  a={"source_variant_38":"chr11_1183604_G_A","QTL_beta_ALT":0.2}
  b={"source_variant_38":"chr11_1184723_G_A","QTL_beta_ALT":-0.3}
  self.assertEqual(m.insert_uniquely_mapped_variant(records,conflicts,key,a,stats),
     "UNIQUE_GRCH37_MAPPED_SNP")
  self.assertEqual(m.insert_uniquely_mapped_variant(records,conflicts,key,b,stats),
     "QUARANTINED_TWO_DISTINCT_SOURCE_SNP")
  self.assertNotIn(key,records)
  self.assertIn(key,conflicts)
  self.assertEqual(m.insert_uniquely_mapped_variant(records,conflicts,key,a,stats),
     "QUARANTINED_AFTER_COLLISION")
  self.assertEqual(stats["GRCh38_TO_GRCh37_DIFFERENT_SOURCE_SNP_COLLISIONS"],1)
 def test_never_promote_stroke_causal_gene_from_same_SNP_associations(self):
  old_assays=m.ASSAYS;old_egenes=m.EGENES
  try:
   m.ASSAYS=1;m.EGENES=1
   assay=("IS_XDATA_G0015","QTD000131","ENSG1")
   gene={"region_id":assay[0],"QTL_dataset":assay[1],"gene_id":assay[2],
     "gene_symbol":"TEST","archived_646_gene_tissue_ABF_count":"0",
     "official_p_beta":"1e-10","BH_qvalue_across_ALL_STUDY_tested_genes":"1e-8",
     "official_n_cis_tested_genetic_variants":"1"}
   previous={
    assay:{"exact_NONPALINDROMIC_SNP_in_BOTH_BBJ_AIS":"1",
      "exact_GRCh37_SNP_in_BOTH_BBJ_AIS":"1",
      "unique_lifted_valid_QTL_SNP_with_any_GWAS":"1"}
   }
   variant="10:101:T:C"
   q={"source_sig":("0.2","1e-10","0.2","0.02","20","100"),
     "source_variant_38":"chr10_108_T_C",
     "mapped_variant_37":variant,"QTL_p":1e-10,"QTL_beta_ALT":.2,
     "QTL_se":.02,"QTL_AF_ALT":.2,"QTL_maf":.2,
     "strand":"+","palindrome":False}
   gwas={"variant":variant,"beta_alt":.3,"se":.1,"p":1e-10,"af_alt":.22,
    "palindrome":False}
   cnt=Counter({"QTL_SNP_VERIFIED_NUMERIC":1})
   rows,genes,varrows,summary=m.integrate([gene],previous,
     {(assay,variant):q},{variant:gwas},{variant:gwas},{assay:cnt},None)
   self.assertEqual(summary["sameSNP_marginal_threshold_joint_assays"],1)
   self.assertEqual(genes[0]["new_H4_coloc_or_causal_genes"],0)
   self.assertFalse(rows[0]["any_eGene_to_AIS_coloc_computed"])
   self.assertFalse(rows[0]["source_QTL_genotype_matched_LD_verified"])
  finally:
   m.ASSAYS=old_assays;m.EGENES=old_egenes
 def test_narrow_source_gwas_overlap_count_must_be_preserved(self):
  old_assays=m.ASSAYS;old_egenes=m.EGENES
  try:
   m.ASSAYS=1;m.EGENES=1
   assay=("IS_XDATA_G0015","QTD000131","ENSG1")
   f=[{"region_id":assay[0],"QTL_dataset":assay[1],"gene_id":assay[2],
      "gene_symbol":"TEST","archived_646_gene_tissue_ABF_count":"0",
      "official_p_beta":"0.001","BH_qvalue_across_ALL_STUDY_tested_genes":"0.02",
      "official_n_cis_tested_genetic_variants":"1"}]
   past={assay:{"exact_NONPALINDROMIC_SNP_in_BOTH_BBJ_AIS":"1",
     "exact_GRCh37_SNP_in_BOTH_BBJ_AIS":"1",
     "unique_lifted_valid_QTL_SNP_with_any_GWAS":"1"}}
   with self.assertRaises(ValueError):m.integrate(f,past,{}, {}, {},{assay:Counter()},None)
  finally:
   m.ASSAYS=old_assays;m.EGENES=old_egenes
if __name__=="__main__":unittest.main()
