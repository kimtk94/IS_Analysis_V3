"""Fail-closed frozen science gates for paired COL4A1/COL4A2 ABF (offline CI).

The 270 posterior rows are diagnostic single-signal results, not causal gene
claims. Prevent accidentally promoting prior-sensitive H4 for either gene.
"""
import csv,json,math,unittest
from collections import defaultdict
from pathlib import Path
BASE=Path(__file__).resolve().parents[1]/"docs/data/is"
def read(name):
 with (BASE/name).open(newline="") as f:
  return list(csv.DictReader(f,delimiter="\t"))
class TestCOL4ABFPriorSensitivity(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.posterior=read("IS_COL4A1_COL4A2_THREE_QTL_P12_AF_N_ABF_SENSITIVITY.tsv")
  cls.summary=read("IS_COL4A1_COL4A2_6_TISSUE_PAIRED_ABF_EXPLORATORY_CONCLUSIONS.tsv")
  cls.tags=read("IS_COL4A1_COL4A2_THREE_TAG_SNP_PUBLIC_SOURCE_COMPARISON.tsv")
  cls.gate=json.loads((BASE/"IS_COL4A1_COL4A2_PAIRED_ABF_SCIENCE_GATE_MANIFEST.json").read_text())
 def test_270_not_causal_and_probabilities_sum_one(self):
  self.assertEqual(len(self.posterior),270)
  self.assertEqual(self.gate["new_shared_causal_genes"],0)
  self.assertEqual(self.gate["three_tissue_matched_6_assays"],6)
  self.assertEqual(len(self.summary),6)
  for row in self.posterior:
   pp=[float(row[f"PP_H{i}"]) for i in range(5)]
   self.assertTrue(all(math.isfinite(x) and 0<=x<=1 for x in pp))
   self.assertAlmostEqual(sum(pp),1,places=7)
   self.assertEqual(row["co_localized_gene_claim"],"NOT_ESTABLISHED")
   self.assertEqual(row["QTL_EAS_ancestry_matched_LD_attested"].upper(),"FALSE")
   self.assertEqual(row["EAS_BBJ_independent_replication_attested"].upper(),"FALSE")
 def test_exact_three_study_two_gene_same_SNP_sample_universe(self):
  expected={"QTD000131":4730,"QTD000171":4648,"QTD000434":3946}
  for did,n in expected.items():
   x=[r for r in self.summary if r["QTL_dataset"]==did]
   self.assertEqual(len(x),2)
   self.assertEqual(set(r["gene_symbol"] for r in x),{"COL4A1","COL4A2"})
   self.assertEqual(set(int(r["matched_nonpalindromic_genomic_SNPs"]) for r in x),{n})
   self.assertTrue(all(r["gene_causal_status"]=="NOT_ESTABLISHED_EXPLORATORY_SINGLE_CAUSAL_ABF_ONLY" for r in x))
 def test_eGene_FDR_and_H4_do_not_select_single_gene(self):
  x={r["gene_symbol"]:r for r in self.summary if r["QTL_dataset"]=="QTD000434"}
  self.assertEqual(x["COL4A1"]["EBI_studywide_BH_eGene"],"False")
  self.assertEqual(x["COL4A2"]["EBI_studywide_BH_eGene"],"True")
  self.assertAlmostEqual(float(x["COL4A1"]["primary_coloc_H4"]),.948198081940169,places=7)
  self.assertAlmostEqual(float(x["COL4A2"]["primary_coloc_H4"]),.95711866720915,places=7)
  self.assertEqual(x["COL4A1"]["main_conditional_H4_SNP"],x["COL4A2"]["main_conditional_H4_SNP"])
  self.assertTrue(self.gate["high_H4_in_two_genes_same_GWAS_SNP_set_is_not_locus_gene_resolution"])
 def test_prior_grid_sensitivity_is_not_hard_validation(self):
  subset=[r for r in self.posterior if
    r["QTL_dataset"]=="QTD000434" and
    r["gene_symbol"]=="COL4A2" and
    r["SNP_selection"]=="ALL_NONPAL_PAIRED" and
    float(r["QTL_N_multiplier"])==1.0]
  self.assertEqual(len(subset),5)
  h4={float(r["p12"]):float(r["PP_H4"]) for r in subset}
  self.assertAlmostEqual(h4[1e-6],.690595693425709,places=7)
  self.assertAlmostEqual(h4[1e-4],.995539730656256,places=7)
 def test_published_lacunar_rs9515201_distinct_from_our_eas_tags(self):
  self.assertEqual(len(self.tags),18)
  genes={r["gene_symbol"] for r in self.tags}
  self.assertEqual(genes,{"COL4A1","COL4A2"})
  x=[r for r in self.tags if r["GRCh37_fasta_REF_ALT_variant"]=="13:111040798:A:C"]
  self.assertEqual(len(x),6)
  self.assertTrue(all(abs(float(r["EAS_AIS_p"])-.03702)<1e-10 for r in x))
  self.assertTrue(all(float(r["QTL_p"])>.05 for r in x))
  self.assertFalse(self.gate["confirmed_sample_ancestry_matched_QTL_LD"])
if __name__=="__main__":unittest.main()
