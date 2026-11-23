"""Synthetic scientific and source-integrity tests for original QTL AN-based coloc pathway."""
from __future__ import annotations

import gzip
import importlib
import json
from pathlib import Path
import unittest

M=importlib.import_module("scripts.is.audit_p3a_official_coloc_ana_qtl_source_probe")


def official_texts():
    columns=("# **an** - total alleles. For autosomal variants, sample size is equal to an/2.\n"
             "# **alt** - also the effect allele\n**maf**\n**beta**\n**se**")
    code=('eQTL_dataset=list(N = (eqtl_sumstats$an)[1]/2, '
          'varbeta = eqtl_sumstats$se^2, MAF = eqtl_sumstats$maf, '
          'beta = eqtl_sumstats$beta, type = "quant")')
    return columns, code, code


def source_gzip(rows=3, an=1120, maf=.21, nonfinite=False, missing=None):
    fields=["molecular_trait_id","gene_id","variant","chromosome","position",
            "ref","alt","maf","an","beta","se","pvalue","rsid","type"]
    if missing:
        fields.remove(missing)
    data=["\t".join(fields)]
    for i in range(rows):
        v={
            "molecular_trait_id":"ENSG00000148795","gene_id":"ENSG00000148795",
            "variant":"chr10_102830531_C_G","chromosome":"10","position":str(102830531+i),
            "ref":"C","alt":"G","maf":str(maf),"an":str(an),"beta":"0.04",
            "se":"nan" if nonfinite else "0.01", "pvalue":"0.005", "rsid":"rs123", "type":"SNP"
        }
        data.append("\t".join(v[f] for f in fields))
    return gzip.compress(("\n".join(data)+"\n").encode())


def manifest():
    return json.loads(Path(
        "docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json"
    ).read_text(encoding="utf-8"))


def preview():
    return M.parse_qtl_header_rows(source_gzip())


def ftp():
    return {"manifest_dataset_id":M.DATASET_ID,
            "official_manifest_contains_full_cis_sumstats_path":True}


class FakeHTTPResponse:
    def __init__(self,status=206,content_range="bytes 0-131071/50000000",
                 url=M.REPO):
        self.status=status
        self.headers={"Content-Range":content_range}
        self._url=url
    def geturl(self):
        return self._url
    def read(self,n):
        return b"\x1f\x8b"+b"q"*(M.MAX_GZIP_RANGE-2)
    def __enter__(self): return self
    def __exit__(self,*args):return False


class OfficialQTLOriginalSourceAuditTests(unittest.TestCase):
    def test_original_coloc_explicit_an_div2_maf_beta_varbeta(self):
        c,t,i=official_texts()
        a=M.official_method_checks(c,t,i)
        self.assertTrue(a["QTL_beta_SE_and_MAF_proven_official_fields"])
        self.assertTrue(a["binary_stroke_GWAS_requires_type_cc_and_verified_case_fraction_s"])
        self.assertFalse(a["AN_div2_validated_exact_regression_sample_N_for_QTD000434"])

    def test_no_false_individual_phenotype_sample_size_approval(self):
        a=M.assemble(manifest(),M.official_method_checks(*official_texts()),
                     ftp(),preview(),{})
        self.assertEqual(a["n_original_gene_tissue_groups"],36)
        self.assertEqual(a["n_original_ROSMAP_assays"],15)
        self.assertEqual(a["n_approved_QTL_N_via_AN_div2"],0)
        self.assertEqual(a["n_formal_coloc_abf_executed"],0)
        self.assertTrue(all(not g["coloc_abf_eligible"] for g in a["groups"]))

    def test_source_first_qtl_header_and_preview_are_not_genewise(self):
        a=preview()
        self.assertEqual(a["preview_rows_checked"],3)
        self.assertEqual(a["preview_sample_count_min_AN_over_2"],560)
        self.assertTrue(a["preview_is_NOT_chr10_15_gene_per_locus_verification"])
        self.assertEqual(a["actual_15_ROSMAP_gene_QTL_AN_MAF_beta_SE_verified"],0)

    def test_partial_gzip_prefix_works(self):
        data=source_gzip(rows=3)
        a=M.parse_qtl_header_rows(data[:len(data)-3])
        self.assertGreaterEqual(a["preview_rows_checked"],2)

    def test_variants_with_invalid_AN_rejected(self):
        for an in [0,-2,1119]:
            with self.subTest(an=an),self.assertRaises(ValueError):
                M.parse_qtl_header_rows(source_gzip(an=an))

    def test_nonfinite_beta_se_rejected(self):
        with self.assertRaises(ValueError):
            M.parse_qtl_header_rows(source_gzip(nonfinite=True))

    def test_maf_outside_binomial_boundary_rejected(self):
        for maf in [0,-.1,.51]:
            with self.subTest(maf=maf),self.assertRaises(ValueError):
                M.parse_qtl_header_rows(source_gzip(maf=maf))

    def test_missing_AN_or_maf_field_rejected(self):
        for field in ["an","maf","beta","se","molecular_trait_id","alt"]:
            with self.subTest(field=field),self.assertRaises(ValueError):
                M.parse_qtl_header_rows(source_gzip(missing=field))

    def test_corrupt_gzip_rejected(self):
        with self.assertRaises(ValueError):
            M.parse_qtl_header_rows(b"plain text content")

    def test_wrong_original_coloc_sample_size_formula_rejected(self):
        c,t,i=official_texts()
        with self.assertRaises(ValueError):
            M.official_method_checks(c,t.replace("$an","$maf"),i)

    def test_wrong_original_coloc_beta_varbeta_rejected(self):
        c,t,i=official_texts()
        with self.assertRaises(ValueError):
            M.official_method_checks(c,t,i.replace("eqtl_sumstats$se^2","eqtl_sumstats$beta^2"))

    def test_wrong_effect_allele_source_claim_rejected(self):
        c,t,i=official_texts()
        with self.assertRaises(ValueError):
            M.official_method_checks(c.replace("also the effect allele","different allele"),t,i)

    def test_QTL_only_does_not_approve_binary_GWAS(self):
        a=M.assemble(manifest(),M.official_method_checks(*official_texts()),
                     ftp(),preview(),{})
        self.assertEqual(a["n_original_GIGASTROKE_GWAS_beta_logOR_and_case_fraction_verified"],0)
        self.assertIn("STROKE_GWAS",a["scientific_status"])

    def test_bounded_remote_only_known_url_with_HTTP_206(self):
        url=("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
             "QTS000025/QTD000434/QTD000434.all.tsv.gz")
        def good_opener(req,timeout):
            self.assertEqual(req.headers["Range"],"bytes=0-131071")
            self.assertEqual(timeout,40)
            return FakeHTTPResponse(url=url)
        raw,metadata=M.bounded_source_get(url,good_opener)
        self.assertEqual(len(raw),M.MAX_GZIP_RANGE)
        self.assertFalse(metadata["whole_source_file_verified_or_downloaded"])
        with self.assertRaises(ValueError):
            M.bounded_source_get(url,lambda req,timeout:FakeHTTPResponse(status=200,url=url))

    def test_remote_redirect_other_domain_rejected(self):
        url=("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
             "QTS000025/QTD000434/QTD000434.all.tsv.gz")
        with self.assertRaises(ValueError):
            M.bounded_source_get(url,lambda req,timeout:FakeHTTPResponse(
                url="https://example.com/elsewhere.gz"))

    def test_non_standard_cross_study_QTL_path_rejected(self):
        url=("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
             "QTS000025/QTD000999/QTD000999.all.tsv.gz")
        with self.assertRaises(ValueError):
            M.bounded_source_get(url)

    def test_invalid_preview_does_not_promote_full_cis(self):
        m=manifest();m["groups"]=m["groups"][:-1]
        with self.assertRaises(ValueError):
            M.assemble(m,M.official_method_checks(*official_texts()),ftp(),preview(),{})

    def test_no_individual_snp_rows_saved(self):
        a=M.assemble(manifest(),M.official_method_checks(*official_texts()),
                     ftp(),preview(),{})
        s=json.dumps(a)
        self.assertNotIn("chr10_102830531_C_G",s)
        self.assertNotIn("rs123",s)
        self.assertFalse(a["donor_IDS_or_SNP_effect_rows_saved"])


if __name__=="__main__":
    unittest.main()
