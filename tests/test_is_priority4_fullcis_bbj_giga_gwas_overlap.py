import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
import audit_is_priority4_fullcis_bbj_giga_gwas_overlap as mod

class FakeCanonical:
 @staticmethod
 def canonical_qtl(chrom,pos,ref,alt,lo,fa,cache):
  if len(ref)==len(alt)==1:
   return "GRCH37_FASTA_REF_VERIFIED",f"{chrom}:{pos}:{ref}:{alt}","+"
  return "QTL_NOT_BIALLELIC_SNP","",None

class TestPriority4GwasCis(unittest.TestCase):
 def test_provenance_and_biallelic_filter(self):
  rows=[
   ["FGF5","ENSG00000138675","BBJ_IS_L001","4","100","A","C","4:100:A:C","True","0.01","0.3","0.1","EBI_REPROCESSED_CIS_TESTED_VARIANT_ONLY"],
   ["NEURL","ENSG00000107954","BBJ_IS_L002","10","200","G","GA","10:200:G:GA","False","0.9","0.1","0.2","EBI_REPROCESSED_CIS_TESTED_VARIANT_ONLY"]]
  columns=["symbol","gene_id","locus","chromosome","position_GRCh38","REF_GRCh38","ALT_GRCh38",
  "variant_id_GRCh38","in_original_archived_GWAS_QTL_intersection","QTD000161_nominal_QTL_p","QTD000161_beta","QTD000161_se","status"]
  with tempfile.TemporaryDirectory() as t,patch.dict(mod.EXPECTED,{"FGF5":1,"NEURL":1},clear=True):
   fp=Path(t)/"input.tsv"
   with fp.open("w",newline="") as f:
    w=csv.writer(f,delimiter="\t");w.writerow(columns);w.writerows(rows)
   inp,mapping,qc,sources=mod.prepare(fp,None,None,FakeCanonical)
   self.assertEqual(len(inp),2)
   self.assertEqual(list(mapping),["4:100:A:C"])
   self.assertEqual(qc["FGF5"]["GRCH37_FASTA_REF_VERIFIED"],1)
   self.assertEqual(qc["NEURL"]["QTL_NOT_BIALLELIC_SNP"],1)
   self.assertEqual(mapping["4:100:A:C"]["records"][0]["in_archived_original_QTL_BBJ_intersection"],"True")
   self.assertFalse(mapping["4:100:A:C"]["records"][0]["same_variant_QTL_GIGA_joint_marginal_threshold"])

 def test_science_manifest_if_locally_available(self):
  p=Path("/srv/is-analysis/results/is/audits/is_priority4_fullcis_bbj_giga_snp_exact_20261011_v1/IS_PRIORITY4_FULLCIS_EBI_BBJ_GIGA_MANIFEST.json")
  if not p.is_file():self.skipTest("Read-only GWAS/EBI science manifest not available")
  data=json.loads(p.read_text())
  self.assertEqual(data["gene_source_records"],20689)
  self.assertEqual(len(data["summaries"]),4)
  self.assertFalse(data["new_coloc_performed"])
  self.assertFalse(data["GTEx_in_study_molecular_LD_verified"])
  self.assertEqual(data["causal_gene_count"],0)
  self.assertTrue(data["original_2225_positional_candidate_genes_unchanged"])
  for r in data["summaries"]:
   self.assertEqual(r["full_ebi_reprocessed_cis_snps"],mod.EXPECTED[r["symbol"]])
   self.assertLessEqual(r["exact_BOTH_BBJ_GIGA_overlap"],r["exact_BBJ_GWAS_overlap"])
   self.assertLessEqual(r["exact_BOTH_BBJ_GIGA_overlap"],r["exact_GIGASTROKE_EAS_AIS_overlap"])
   self.assertLessEqual(r["EBI_only_with_BBJ_GWAS"],r["EBI_only_vs_original_QTL_BBJ_intersection_count"])
if __name__=="__main__":unittest.main()
