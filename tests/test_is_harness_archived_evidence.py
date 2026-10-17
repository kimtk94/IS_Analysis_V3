"""Offline independent checks of archived EAS5 FASTA/QTL evidence flags."""
import importlib.util
import pathlib
import tempfile
import unittest
import csv
import gzip

SCRIPT=pathlib.Path(__file__).resolve().parents[1]/"scripts/is/harness_archived_evidence.py"
spec=importlib.util.spec_from_file_location("is_archived_evidence_test",SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def write_tsv(path,columns,items):
    path.parent.mkdir(parents=True,exist_ok=True)
    op=gzip.open if path.suffix==".gz" else open
    with op(path,"wt",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=columns,delimiter="\t")
        writer.writeheader()
        writer.writerows(items)

class ArchivedEvidenceTests(unittest.TestCase):
    def test_source_flags_are_not_new_fasta_validation_or_coloc(self):
        with tempfile.TemporaryDirectory() as directory:
            root=pathlib.Path(directory)
            archived=root/mod.ARCHIVED_RELATIVE
            ready=root/mod.READY_RELATIVE
            overlap=root/mod.OVERLAP_RELATIVE
            fields=sorted(mod.REQUIRED_ARCHIVE)
            rows=[
                {"region_id":"R1","QTL_dataset":"Q1","gene_id":"ENSG1",
                 "mapped_genomic_variant_GRCh37":"1:100:A:C",
                 "BBJ_allele_exact_Fasta37":"True","EAS_AIS_allele_exact_Fasta37":"True",
                 "QTL_p":"0.001","BBJ_p":"0.0001","EAS_AIS_p":"0.02",
                 "QTL_beta_GRCh37_ALT":"0.3","BBJ_beta_ALT":"0.2","EAS_AIS_beta_ALT":"0.1"},
                {"region_id":"R1","QTL_dataset":"Q1","gene_id":"ENSG1",
                 "mapped_genomic_variant_GRCh37":"1:101:A:G",
                 "BBJ_allele_exact_Fasta37":"False","EAS_AIS_allele_exact_Fasta37":"True"},
                {"region_id":"R2","QTL_dataset":"Q2","gene_id":"ENSG2",
                 "mapped_genomic_variant_GRCh37":"1:100:A:C",
                 "BBJ_allele_exact_Fasta37":"True","EAS_AIS_allele_exact_Fasta37":"True"},
            ]
            write_tsv(archived,fields,rows)
            write_tsv(overlap,sorted(mod.REQUIRED_OVERLAP),[
                {"region_id":"R1","QTL_dataset":"Q1","gene_id":"ENSG1",
                 "GRCh37_FASTA_VERIFIED_SNP_rows":"2","FASTA_REF_CONFLICT_rows":"0",
                 "both_GWAS_exact_variant_rows":"2","both_GWAS_exact_nonpalindromic_rows":"1",
                 "eligible_for_coloc":"False","blocked_on":"cohort-matched LD"},
            ])
            write_tsv(ready,sorted(mod.REQUIRED_READY),[
                {"region_id":"R1","QTL_dataset":"Q1","gene_id":"ENSG1",
                 "possible_single_signal_ABF_followup":"True",
                 "can_now_assert_source_genome_GRCh37_REF_after_liftover":"False",
                 "can_now_assert_OFFICIAL_ALT_EFFECT_DIRECTION":"False",
                 "QTL_source_ancestry_verified":"False",
                 "molecular_QTL_matched_LD_verified":"False"},
            ])
            r=mod.audit(root)
            self.assertEqual(r["archived"]["archived_rows"],3)
            self.assertEqual(r["archived"]["both_archived_fasta_flag_rows"],2)
            self.assertEqual(r["archived"]["both_fasta_flag_unique_variant_strings"],1)
            self.assertEqual(r["overlap_report"]["reported_fasta_verified_rows_sum"],2)
            self.assertEqual(r["readiness"]["abf_candidate_rows"],1)
            self.assertIs(r["abf_flagged_groups"][0]["matched_qtl_ld_verified"],False)
            self.assertFalse(any("coloc_executed" in name for name in r))

    def test_false_string_and_unrecognized_flags_are_not_truthy(self):
        self.assertIs(mod.flag("False"),False)
        self.assertIs(mod.flag("True"),True)
        self.assertIsNone(mod.flag("mystery"))
        self.assertIsNone(mod.nonnegative("-1"))
        self.assertIsNone(mod.nonnegative("not-a-count"))

if __name__=="__main__":
    unittest.main()
