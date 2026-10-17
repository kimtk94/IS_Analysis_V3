import csv
import tempfile
import unittest
from pathlib import Path
from scripts.is.build_is_priority4_evidence_scorecard import integrate,GENES,EXPECTED,LOCUS

class TestEvidenceScorecard(unittest.TestCase):
    def make_comparison(self,p):
        columns=["gene","locus","QTD_nominal_p_lt_1e4",
          "both_GWAS_allele_exact_nonpal_unique","BBJ_and_GIGA_GWS_same_SNP",
          "triple_QTD_BBJ_GIGA_threshold"]
        with p.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=columns,delimiter="\t");w.writeheader()
            for gene in GENES:
                w.writerow(dict(gene=gene,locus=LOCUS[gene],
                  QTD_nominal_p_lt_1e4=EXPECTED[gene],
                  both_GWAS_allele_exact_nonpal_unique=2000,
                  BBJ_and_GIGA_GWS_same_SNP=0,
                  triple_QTD_BBJ_GIGA_threshold=0))
    def test_missing_stages_are_not_promoted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=root/"compare.tsv";self.make_comparison(p)
            rows=integrate(p,root/"missing_p01.tsv",root/"missing_p02.tsv",root/"result")
            self.assertEqual(len(rows),4)
            cal=next(x for x in rows if x["gene"]=="CALHM2")
            self.assertEqual(cal["P01_status"],"NOT_EXECUTED_OR_MISSING")
            self.assertEqual(cal["P02_status"],"NOT_EXECUTED_OR_MISSING")
            self.assertEqual(cal["causal_gene_status"],"NOT_ESTABLISHED")
            with self.assertRaises(FileExistsError):
                integrate(p,root/"missing_p01.tsv",root/"missing_p02.tsv",root/"result")
    def test_missing_comparison_is_blocked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            with self.assertRaises(FileNotFoundError):
                integrate(root/"none",root/"missing1",root/"missing2",root/"out")
    def test_unexpected_qtl_count_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=root/"compare.tsv";self.make_comparison(p)
            text=p.read_text().replace("CALHM2\tBBJ_IS_L002\t3","CALHM2\tBBJ_IS_L002\t4")
            p.write_text(text)
            with self.assertRaises(ValueError):
                integrate(p,root/"missing1",root/"missing2",root/"out")
