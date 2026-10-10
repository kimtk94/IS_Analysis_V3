import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_priority4_qtd_full_cis_reconstruction import construct

class TestFullCisReconstruction(unittest.TestCase):
    def test_each_gene_unique_cis_equals_independent_permuted_denominator(self):
        root=Path("/srv/is-analysis/results/is/audits")
        qdir=root/"is_priority4_qtd000161_full_cis_wide_region_20261011_v2"
        perm=root/"is_priority4_qtd000161_full_cis_universe_20261011_v1/QTD000161.permuted.tsv.gz"
        src=root/"legacy_coloc_drive_recovery_20261011_v1/input"
        if not all(x.exists() for x in [qdir/"BBJ_IS_L001.QTD000161.wide.priority4.tsv",qdir/"BBJ_IS_L002.QTD000161.wide.priority4.tsv",perm]):
            self.skipTest("Reprocessed wide-cis public inputs unavailable")
        with tempfile.TemporaryDirectory() as tmp:
            j=construct(qdir,perm,src,Path(tmp))
            self.assertEqual(j["genes"],4)
            self.assertEqual(j["complete_ebi_reprocessed_universes"],4)
            self.assertEqual(j["EBI_reprocessed_gene_variant_rows"],20689)
            self.assertFalse(j["GTEx_original_FastQTL_test_universe_proven"])
            self.assertFalse(j["archived_GTEx_V8_label_source_pipeline_authenticated"])
            self.assertFalse(j["can_promote_coloc_susie"])
            self.assertTrue(j["all_2225_positional_genes_preserved"])
            data={r["symbol"]:r for r in j["results"]}
            counts={"FGF5":5511,"C4orf22":5490,"CALHM2":4868,"NEURL":4820}
            self.assertEqual({k:x["tabix_wide_unique_variant_count"] for k,x in data.items()},counts)
            for row in data.values():
                self.assertGreater(row["exact_shared_variants"],0)
                self.assertGreater(row["EBI_reprocessed_only_variant_keys"],0)
                self.assertEqual(row["permutation_minus_wide_unique"],0)
if __name__=="__main__":unittest.main()
