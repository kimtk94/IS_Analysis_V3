import gzip
import importlib.util
from pathlib import Path
import unittest

FILE=Path(__file__).resolve().parents[1]/"scripts/is/audit_qtd000161_historic_permuted_sha.py"
spec=importlib.util.spec_from_file_location("qtd_perm",FILE)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

HEADER="molecular_trait_id\tn_variants\tp_beta\n"
def fixture(change=None):
    lines=[]
    for gene,(name,count,p) in m.HISTORICAL_GENES.items():
        n,px=(count,p)
        if change and gene in change:n,px=change[gene]
        lines.append(f"{gene}\t{n}\t{px}")
    return gzip.compress((HEADER+"\n".join(lines)+"\n").encode(),mtime=0)

class Tests(unittest.TestCase):
    def test_actual_archive_sha_is_pinned(self):
        self.assertEqual(len(m.HISTORICAL_SHA256),64)
        self.assertEqual(sum(v[1] for v in m.HISTORICAL_GENES.values()),20689)
    def test_all_genes_correct_but_unknown_bytes_still_hold(self):
        r=m.inspect(fixture())
        self.assertTrue(r["all_four_gene_permutation_pbeta_n_variants_match_history"])
        self.assertFalse(r["whole_permuted_gzip_same_bytes_as_20261011_archive"])
        self.assertIn("HOLD",r["permutation_source_gate"])
    def test_wrong_beta_p_holds(self):
        change={"ENSG00000138675":(5511,0.4)}
        r=m.inspect(fixture(change))
        self.assertFalse(r["all_four_gene_permutation_pbeta_n_variants_match_history"])
        self.assertFalse(r["four_gene_permutation_comparison"]["FGF5"]["p_beta_close"])
    def test_wrong_n_holds(self):
        r=m.inspect(fixture({"ENSG00000107954":(4819,0.40605)}))
        self.assertFalse(r["all_four_gene_permutation_pbeta_n_variants_match_history"])
    def test_missing_gene_holds(self):
        raw=gzip.compress((HEADER+"ENSG00000138675\t5511\t0.0278074\n").encode())
        r=m.inspect(raw)
        self.assertFalse(r["all_four_gene_permutation_pbeta_n_variants_match_history"])
    def test_source_invalid_gzip_holds(self):
        r=m.run(lambda _:b"wrong-gzip-header")
        self.assertEqual(r["permutation_source_gate"],"HOLD_OFFICIAL_PERMUTATION_SOURCE_UNAVAILABLE")
        self.assertFalse(r["coloc_performed"])
    def test_unapproved_large_file_holds(self):
        old=m.MAX_COMPRESSED
        try:
            m.MAX_COMPRESSED=5
            result=m.run(lambda _:fixture())
            self.assertIn("HOLD",result["permutation_source_gate"])
        finally:m.MAX_COMPRESSED=old
    def test_duplicate_gene_is_rejected(self):
        v="ENSG00000138675\t5511\t0.0278074\n"
        data=gzip.compress((HEADER+v+v).encode())
        with self.assertRaises(ValueError):m.inspect(data)
    def test_no_coloc_promotion_if_hash_or_stats_fail(self):
        r=m.inspect(fixture())
        self.assertFalse(r["QTD_all_full_nominal_original_2023_gzip_sha256_authenticated"])
        self.assertFalse(r["is_GTEx_v8_original_FastQTL_statistics"])
        self.assertEqual(r["scientific_status"],"HOLD")
if __name__=="__main__":
    unittest.main()
