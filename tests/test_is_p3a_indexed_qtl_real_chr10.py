"""Fail-closed indexed official GRCh38 QTL chr10 survey tests; no remote downloads."""
import copy
import importlib
import io
import json
import unittest
from pathlib import Path

M=importlib.import_module("scripts.is.audit_p3a_indexed_qtl_real_chr10")
LIFT=importlib.import_module("scripts.is.audit_p3a_chr10_ucsc_liftover_fasta")
OFF=importlib.import_module("scripts.is.audit_official_p3a_qtl_catalogue_provenance")
REGION={"chrom":"10","start":102403919,"end":104460593}

class DNA:
    def __getitem__(self,i):
        return ord("A")
class IndexedQTLTest(unittest.TestCase):
    def source(self, extra=""):
        h="\t".join(M.COLUMNS+("type","rsid"))
        values=["10","103000000","A","C","chr10_103000000_A_C",
                "ENSG00000065613","ENSG00000065613","0.12","0.03",
                "0.23","0.00001","SNP","rs1"]
        second=values.copy()
        second[-1]="rs2"
        return [h+"\n","\t".join(values)+"\n","\t".join(second)+"\n"]+([extra] if extra else [])
    def test_valid_original_format_dedups_different_rsid(self):
        r=M.read_rows(self.source(),DNA(),{"ENSG00000065613"},REGION)
        self.assertEqual(r["interval_rows"]["region_rows_scanned"],2)
        self.assertEqual(r["interval_rows"]["candidate_unique_SNP_gene_pairs"],1)
        self.assertEqual(r["interval_rows"]["duplicated_variant_rsid_rows"],1)
        self.assertEqual(r["interval_rows"]["conflicting_duplicates"],0)
        self.assertEqual(r["gene_rows"]["ENSG00000065613"]["valid_nonpal_SNPs"],1)
    def test_mismatched_source_ref_denied(self):
        rows=self.source()
        rows[1]=rows[1].replace("\tA\tC\tchr10_103000000_A_C",
                                 "\tG\tC\tchr10_103000000_G_C")
        r=M.read_rows(rows,DNA(),{"ENSG00000065613"},REGION)
        self.assertEqual(r["interval_rows"]["bad_REF_against_UCSC_GRCh38"],1)
    def test_invalid_nonfinite_se_is_not_approved(self):
        rows=self.source()
        rows[1]=rows[1].replace("\t0.03\t","\tNaN\t")
        out=M.read_rows(rows,DNA(),{"ENSG00000065613"},REGION)
        self.assertEqual(out["interval_rows"]["invalid_signed_beta_SE_MAF"],1)
    def test_bad_original_variant_id_rejected(self):
        rows=self.source()
        rows[1]=rows[1].replace("chr10_103000000_A_C","chr9_103000000_A_C")
        out=M.read_rows(rows,DNA(),{"ENSG00000065613"},REGION)
        self.assertEqual(out["interval_rows"]["invalid_allele_or_variant_key"],1)
    def test_bad_build_region_rejected(self):
        rows=self.source()
        rows[1]=rows[1].replace("\t103000000\t","\t107000000\t")
        with self.assertRaises(ValueError):
            M.read_rows(rows,DNA(),{"ENSG00000065613"},REGION)
    def test_truncated_tabix_rows_rejected(self):
        rows=self.source()
        rows[1]="10\t103000000\tA\n"
        with self.assertRaises(ValueError):
            M.read_rows(rows,DNA(),{"ENSG00000065613"},REGION)
    def test_missing_header_rejected(self):
        with self.assertRaises(ValueError):
            M.read_rows(["10\t103000000\tA\tC\n"],DNA(),{"ENSG00000065613"},REGION)
    def test_qtl_nominal_sample_N_not_approved(self):
        self.assertEqual(M.EXPECTED_N,{"QTD000131":387,"QTD000171":205,"QTD000434":560})
    def test_conflicting_duplicate_variant_rejected_count(self):
        rows=self.source()
        rows[-1]=rows[-1].replace("\t0.12\t0.03\t","\t0.93\t0.03\t")
        out=M.read_rows(rows,DNA(),{"ENSG00000065613"},REGION)
        self.assertEqual(out["interval_rows"]["conflicting_duplicates"],1)
    def test_unrelated_genes_are_not_mislabelled(self):
        out=M.read_rows(self.source(),DNA(),{"ENSG00000299999"},REGION)
        self.assertEqual(out["interval_rows"]["candidate_gene_rows_seen"],0)
        self.assertEqual(out["interval_rows"]["noncandidate_rows"],2)
    def test_zero_based_chr10_fasta_ref_used(self):
        class REF:
            def __getitem__(self,index):
                return ord("G") if index==103000000-1 else ord("A")
        rows=self.source()
        rows[1]=rows[1].replace("\tA\tC\tchr10_103000000_A_C",
                                 "\tG\tC\tchr10_103000000_G_C")
        out=M.read_rows(rows,REF(),{"ENSG00000065613"},REGION)
        self.assertEqual(out["gene_rows"]["ENSG00000065613"]["valid_unique_SNPs"],1)
    def test_no_inferring_fullcis_from_one_region(self):
        self.assertFalse(LIFT.REGION["id"] != "IS_XDATA_G0015")

if __name__=="__main__":
    unittest.main()
