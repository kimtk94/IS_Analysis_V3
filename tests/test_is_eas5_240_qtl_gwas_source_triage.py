"""Regression: never promote study-wide molecular eGene evidence to IS causal genes."""
import importlib.util,sys,unittest
from pathlib import Path
path=Path(__file__).resolve().parents[1]/"scripts/is/build_is_eas5_240_qtl_gwas_source_triage.py"
spec=importlib.util.spec_from_file_location("eas5_genetic_source_triage",path)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
COUNTS={"IS_XDATA_G0015":71,"IS_XDATA_G0023":14,"IS_XDATA_G0007":10,
 "IS_XDATA_G0025":92,"IS_XDATA_G0018":53}
def synthetic():
 genes=[];source=[];region=[]
 i=0
 for rid,n in COUNTS.items():
  for j in range(n):
   pos=i;gid=f"ENSG{pos:011d}";i+=1
   source.append({"region_id":rid,"gene_id":gid,"QTL_dataset":"QTD000131",
                  "exact_NONPALINDROMIC_SNP_in_BOTH_BBJ_AIS":"100"})
   if pos<153:
    source.append({"region_id":rid,"gene_id":gid,"QTL_dataset":"QTD000171",
                   "exact_NONPALINDROMIC_SNP_in_BOTH_BBJ_AIS":"200"})
   bh=1 if pos<66 else 0
   prior=1 if pos<10 else 0
   genes.append({"region_id":rid,"gene_id":gid,
    "gene_symbol":f"TEST{pos}","number_archived_BBJ_GTEx_ABF_assays":str(prior),
    "number_QTL_studies_with_official_fullcis_test":"3",
    "number_QTL_studies_with_genomewide_BH05_eQTL":str(bh),
    "dataset_eGene_FDR05":"QTD000131" if bh else "",
    "lowest_BH_q_across_three_QTL_studies":"0.01" if bh else "0.2"})
 for rid,tier in m.REGION_TIER.items():
  region.append({"region_id":rid,"primary_GWAS_EAS_BBJ_AIS_tier":tier})
 region += [{"region_id":f"OTHER{i:02d}","primary_GWAS_EAS_BBJ_AIS_tier":"T4_OTHER"} for i in range(25)]
 assert len(genes)==240 and len(source)==393 and len(region)==30
 return genes,source,region,{}
class TestPositionalQTLReadiness(unittest.TestCase):
 def test_original_240_gene_set_and_novel_eGene_BH_gate(self):
  genes,alleles,regions,case=synthetic()
  result,summary=m.build(genes,alleles,regions,case)
  self.assertEqual(len(result),240)
  self.assertEqual(len(summary),5)
  self.assertEqual(sum(int(r["source_QTL_eGene_studywide_BH_q_le_0p05"])>0 for r in result),66)
  self.assertEqual(sum(int(r["source_QTL_eGene_studywide_BH_q_le_0p05"])>0 and int(r["original_gene_archived_646_ABF_assays"])==0 for r in result),56)
  self.assertTrue(all(r["gene_causal_status"]=="NOT_ESTABLISHED" for r in result))
  self.assertTrue(all(r["extended_full_cis_SNP_variant_GWAS_matching_rerun"] is False for r in result))
 def test_drifting_region_trait_tier_rejected(self):
  genes,alleles,regions,case=synthetic()
  regions[0]["primary_GWAS_EAS_BBJ_AIS_tier"]="FAKE_SHARED_CAUSAL_LOCUS"
  with self.assertRaises(ValueError):m.build(genes,alleles,regions,case)
 def test_missing_gene_tissue_source_rejected(self):
  genes,alleles,regions,case=synthetic()
  alleles.pop()
  with self.assertRaises(ValueError):m.build(genes,alleles,regions,case)
if __name__=="__main__":unittest.main()
