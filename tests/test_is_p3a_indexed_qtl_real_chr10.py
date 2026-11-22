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
    def test_actual_documented_column_order_without_commented_header(self):
        field={
            "molecular_trait_id":"ENSG00000065613","chromosome":"10",
            "position":"103000000","ref":"A","alt":"C",
            "variant":"chr10_103000000_A_C",
            "ma_samples":"12","maf":"0.23","pvalue":"0.00001",
            "beta":"0.12","se":"0.03","type":"SNP",
            "ac":"22","an":"300","r2":"0.99",
            "molecular_trait_object_id":"ENSG00000065613",
            "gene_id":"ENSG00000065613","median_tpm":"1.3","rsid":"rs1"}
        first="\t".join(M.OFFICIAL_COLUMN_ORDER)
        second="\t".join(field[k] for k in M.OFFICIAL_COLUMN_ORDER)
        z=M.read_rows([first+"\n",second+"\n"],DNA(),
                      {"ENSG00000065613"},REGION)
        self.assertEqual(z["gene_rows"]["ENSG00000065613"]["valid_unique_SNPs"],1)
    def test_no_inferring_fullcis_from_one_region(self):
        self.assertFalse(LIFT.REGION["id"] != "IS_XDATA_G0015")

    def test_indel_is_classified_not_mislabelled_bad_SNP(self):
        rows=self.source()
        cols=rows[1].strip().split("\t")
        h=rows[0].strip().split("\t")
        vals=dict(zip(h,cols))
        vals.update(ref="AA",alt="A",
                    variant="chr10_103000000_AA_A",type="INDEL")
        rows[1]="\t".join(vals[k] for k in h)+"\n"
        x=M.read_rows(rows,DNA(),{"ENSG00000065613"},REGION)
        self.assertEqual(x["interval_rows"]["excluded_variant_reason_classes"][
            "MULTIBASE_ALLELE_INDEL_OR_COMPLEX"],1)
        self.assertEqual(x["interval_rows"]["excluded_variant_source_type_labels"][
            "INDEL"],1)

    def test_variant_key_mismatch_has_own_reason(self):
        rows=self.source()
        rows[1]=rows[1].replace("chr10_103000000_A_C","chr9_103000000_A_C")
        x=M.read_rows(rows,DNA(),{"ENSG00000065613"},REGION)
        self.assertEqual(x["interval_rows"]["excluded_variant_reason_classes"][
            "NONCANONICAL_VARIANT_KEY"],1)
        self.assertEqual(x["interval_rows"]["excluded_variant_reason_classes"][
            "MULTIBASE_ALLELE_INDEL_OR_COMPLEX"],0)

    def test_non_acgt_allele_classification(self):
        rows=self.source()
        rows[1]=rows[1].replace("\tA\tC\tchr10_103000000_A_C",
                                 "\tN\tC\tchr10_103000000_N_C")
        x=M.read_rows(rows,DNA(),{"ENSG00000065613"},REGION)
        self.assertEqual(x["interval_rows"]["excluded_variant_reason_classes"][
            "NON_ACGT_ALLELE"],1)

    def test_equal_reference_and_alternate_classification(self):
        rows=self.source()
        rows[1]=rows[1].replace("\tA\tC\tchr10_103000000_A_C",
                                 "\tA\tA\tchr10_103000000_A_A")
        x=M.read_rows(rows,DNA(),{"ENSG00000065613"},REGION)
        self.assertEqual(x["interval_rows"]["excluded_variant_reason_classes"][
            "SAME_REF_ALT"],1)

    def test_original_excluded_classes_are_mutually_exclusive(self):
        rows=self.source()
        rows[1]=rows[1].replace("\tA\tC\tchr10_103000000_A_C",
                                 "\tAA\tA\tchr10_103000000_AA_A")
        x=M.read_rows(rows,DNA(),{"ENSG00000065613"},REGION)
        v=x["interval_rows"]
        self.assertEqual(sum(v["excluded_variant_reason_classes"].values()),
                         v["invalid_allele_or_variant_key"])

    def test_bounds_are_observed_not_theoretical_fullcis(self):
        x=M.read_rows(self.source(),DNA(),{"ENSG00000065613"},REGION)
        p=x["gene_rows"]["ENSG00000065613"]
        self.assertEqual(p["observed_GRCh38_min_pos"],103000000)
        self.assertEqual(p["observed_GRCh38_max_pos"],103000000)

if __name__=="__main__":
    unittest.main()
