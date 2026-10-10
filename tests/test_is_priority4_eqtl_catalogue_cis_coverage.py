import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_priority4_eqtl_catalogue_cis_coverage import (
    COLS,GENES,parse_qtd_row,get_permuted,evaluate,
)
class TestQTDReprocessedCis(unittest.TestCase):
 def test_fixed_source_fields_and_allele_key(self):
  row=["ENSG00000138675","4","79759933","G","C","chr4_79759933_G_C",
       "131","0.494286","0.0632102","-0.170516","0.09115","SNP",
       "177","350","NA","ENSG00000138675","ENSG00000138675","13.616","rs12511948"]
  self.assertEqual(len(row),len(COLS))
  d=parse_qtd_row("\t".join(row)+"\n")
  self.assertEqual(d["gene_id"],"ENSG00000138675")
  self.assertEqual(d["pvalue"],"0.0632102")
  bad=row.copy();bad[3]="A"
  with self.assertRaises(ValueError):parse_qtd_row("\t".join(bad)+"\n")
  with self.assertRaises(ValueError):parse_qtd_row("\t".join(row[:-1])+"\n")
 def test_permuted_full_cis_denominators_if_mounted(self):
  base=Path("/srv/is-analysis/results/is/audits/is_priority4_qtd000161_full_cis_universe_20261011_v1")
  perm=base/"QTD000161.permuted.tsv.gz"
  if not perm.exists():self.skipTest("Open eQTL Catalogue file unavailable")
  s=get_permuted(perm)
  self.assertEqual({k:int(v["n_variants"]) for k,v in s.items()},
     {"ENSG00000138675":5511,"ENSG00000197826":5490,
      "ENSG00000138172":4868,"ENSG00000107954":4820})
 def test_one_gene_original_archive_vs_reprocessed_not_equivalent_if_mounted(self):
  base=Path("/srv/is-analysis/results/is/audits")
  q=base/"is_priority4_qtd000161_full_cis_universe_20261011_v1"
  files=[q/"BBJ_IS_L001.QTD000161.region.priority4.tsv",
         q/"BBJ_IS_L002.QTD000161.region.priority4.tsv"]
  perm=q/"QTD000161.permuted.tsv.gz"
  src=base/"legacy_coloc_drive_recovery_20261011_v1/input"
  if not all(p.exists() for p in files+[perm]):self.skipTest("Region query not locally available")
  with tempfile.TemporaryDirectory() as t:
   x=evaluate(perm,src,files,Path(t))
   self.assertTrue(x["reprocessed_summary_not_original_gtEx_nominal"])
   self.assertFalse(x["exact_original_archive_coloc_uses_all_tested_variants_verified"])
   self.assertTrue(all(0<=v["shared_SNP_count_div_full_catalogue_cis_count"]<=1 for v in x["rows"]))
   self.assertEqual(x["original_candidate_genes"],2225)
if __name__=="__main__":unittest.main()
