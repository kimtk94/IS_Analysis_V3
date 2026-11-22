"""Original 36-group QTL sample model and SE/MAF/AN variance-proxy negative tests."""
import importlib
import json
import math
import unittest
from pathlib import Path

M=importlib.import_module("scripts.is.audit_p3a_36_qtl_sdY_proxy_input_gate")
QTL=M.QTL
SOURCE_HEADERS=QTL.OFFICIAL_COLUMN_ORDER
GENE="ENSG00000065613"


def tsv(pos=103000000,an="774",ac="154",maf="0.198966",
        beta="0.14",se="0.055",pvalue="0.0003",gene=GENE):
    d={k:"." for k in SOURCE_HEADERS}
    d.update({"molecular_trait_id":gene,"chromosome":"10",
        "position":str(pos),"ref":"A","alt":"C",
        "variant":f"chr10_{pos}_A_C","gene_id":gene,
        "beta":beta,"se":se,"pvalue":pvalue,"maf":maf,
        "type":"SNP","an":an,"ac":ac,
        "ma_samples":"101","rsid":"rs123"})
    return "\t".join(d[k] for k in SOURCE_HEADERS)


class FakeTabix:
    contigs=("10",)
    def __init__(self,lines):self.lines=lines
    def __enter__(self):return self
    def __exit__(self,*args):return False
    def fetch(self,chrom,start,end):
        assert chrom=="10"
        assert start==M.UNION.EXPECTED_WINDOW["start"]-1
        assert end==M.UNION.EXPECTED_WINDOW["end"]
        return iter(self.lines)


def artifact_dict():
    return json.loads(Path("docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json").read_text())


def source_inputs():
    manifest=artifact_dict()
    keys=manifest["groups"]
    union={"schema":"is-p3a-36-expanded-original-cis-reference-union-v1",
           "technical_status":"PASS_THREE_REAL_INDEXED_REFERENCE_UNION_SURVEYS",
           "scientific_status":"HOLD","n_approved_coloc_abf":0,
           "cis_union_reference_Ensembl105_GRCh38_1based":M.UNION.EXPECTED_WINDOW,
           "groups":[{"QTL_dataset":g["QTL_dataset"],"gene_id":g["gene_id"]}
                     for g in keys]}
    allele={"schema":"is-p3a-primary-gwas-yaml-QTL-allele-N-36-v1",
            "technical_status":"PASS_PRIMARY_PUBLISHER_METADATA_AND_36_QTL_AN_SOURCE_READ",
            "formal_coloc_abf_eligible":0,
            "groups":[{"QTL_dataset":g["QTL_dataset"],"gene_id":g["gene_id"],
                       "original_variant_allele_N_summary":{
                           "count_valid_nonpal_original_QTL":1,
                           "count_an_valid_even_autosomal":1,
                           "count_an_missing_or_invalid":0,
                           "min_variant_an_over_2":M.STUDY_N[g["QTL_dataset"]],
                           "max_variant_an_over_2":M.STUDY_N[g["QTL_dataset"]],
                           "n_distinct_variant_an_over_2":1,
                           "top_variant_an_over_2_count":[
                               {"variant_n":M.STUDY_N[g["QTL_dataset"]],"records":1}],
                           "count_QTL_maf_vs_minor_ac_an_delta_gt_0_01":0,
                       }} for g in keys]}
    match={"schema":"is-p3a-36-real-source-reciprocal-SNP-GWAS-harmonization-v1",
           "technical_status":"PASS_QTL_GWAS_SOURCE_SNP_HARMONIZATION_AGGREGATES_ONLY",
           "n_formal_coloc_abf_eligible":0,
           "rows":[{"QTL_dataset":g["QTL_dataset"],"gene_id":g["gene_id"],
                    "QTL_source_indexed_nonpal_valid":1,
                    "QTL_uniquely_reciprocally_lifted_REF_valid":1,
                    "per_ancestry":{
                         a:{"dataset":("GCST90104545" if a=="EAS" else "GCST90104540"),
                            "n_nonpal_QTL_GWAS_recipient_allele_matches":1}
                         for a in ("EAS","EUR")}}
                   for g in keys]}
    cis={"schema":"is-p3a-36-original-eQTL-permuted-n_variants-nominal-completeness-v1",
         "technical_status":"PASS_3_ORIGINAL_PERMUTATION_AND_NOMINAL_REGION_SOURCES_READ",
         "n_gene_tissue_groups_matching_full_cis_tested_variant_COUNT":36,
         "n_formal_coloc_abf_eligible":0,
         "groups":[{"QTL_dataset":g["QTL_dataset"],"gene_id":g["gene_id"],
                    "published_permutation_vs_recovered_nominal_unique_agree":True,
                    "original_permutation_rows_for_trait":1} for g in keys]}
    official={"technical_status":"PASS_3_SOURCE_IDENTITIES_IN_IMMUTABLE_CATALOGUE",
              "QTL_datasets":{q:{"nominal_catalogue_sample_size":n}
                              for q,n in M.STUDY_N.items()}}
    stats={q:{g["gene_id"]:{
           "original_strict_nonpal_SNPs":1,"positive_even_genotyped_AN_count":1,
           "invalid_AN_AC_count":0,"min_observed_SNP_N":M.STUDY_N[q],
           "max_observed_SNP_N":M.STUDY_N[q],
           "n_distinct_observed_SNP_N":1,
           "most_common_N_observed":[{"variant_n":M.STUDY_N[q],"records":1}],
           "maf_vs_source_ac_an_delta_gt_0p01":0}
           for g in keys if g["QTL_dataset"]==q} for q in M.STUDY_N}
    return manifest,union,allele,match,cis,official,stats


class SourceVarianceTests(unittest.TestCase):
    def test_true_residual_se_proxy_positive_and_formula(self):
        x=M.proxy_value(.05,.2,387)
        self.assertAlmostEqual(x,.05*math.sqrt(2*.2*.8*387))
        self.assertGreater(x,0)

    def test_cannot_compute_proxy_with_zero_standard_error(self):
        with self.assertRaises(ValueError):M.proxy_value(0,.2,387)

    def test_cannot_compute_proxy_with_nan_maf(self):
        with self.assertRaises(ValueError):M.proxy_value(.05,float("nan"),387)

    def test_invalid_sample_n_is_not_silently_accepted(self):
        with self.assertRaises(ValueError):M.proxy_value(.05,.2,0)

    def test_quantile_interpolates_numeric_data(self):
        self.assertEqual(M.quantile([1.,2.,3.,4.,5.],.5),3)
        self.assertEqual(M.quantile([1.,2.,3.,4.],.5),2.5)

    def test_summarize_keeps_values_descriptive_only(self):
        x=M.summarize([2.0,4.0,6.0])
        self.assertEqual(x["median"],4)
        self.assertEqual(x["n"],3)

    def test_100pct_genotyped_N_source_and_original_SE(self):
        line=tsv()
        result,provenance=M.scan_one("QTD000131","fake",{GENE},
                                      open_tabix=lambda _:FakeTabix([line]))
        row=result[GENE]
        self.assertEqual(row["original_strict_nonpal_SNPs"],1)
        self.assertEqual(row["min_observed_SNP_N"],387)
        self.assertEqual(row["n_distinct_observed_SNP_N"],1)
        self.assertEqual(row["per_SNP_regression_SE_MAF_AN_residual_proxy_all"]["n"],1)
        self.assertFalse(row["quantitative_source_phenotype_sdY_independently_observed"])
        self.assertEqual(provenance["original_full_region_rows_read"],1)

    def test_source_variant_N_varies_but_study_cohort_N_does_not(self):
        a=tsv(pos=103000000,an="774")
        b=tsv(pos=103000001,an="770")
        result,_=M.scan_one("QTD000131","fake",{GENE},
                              open_tabix=lambda _:FakeTabix([a,b]))
        x=result[GENE]
        self.assertEqual(x["min_observed_SNP_N"],385)
        self.assertEqual(x["max_observed_SNP_N"],387)
        self.assertEqual(x["n_distinct_observed_SNP_N"],2)
        self.assertEqual(x["n_records_at_dataset_nominal_N"],1)

    def test_exact_duplicate_rsid_does_not_inflate_proxy_N(self):
        result,_=M.scan_one("QTD000131","fake",{GENE},
                     open_tabix=lambda _:FakeTabix([tsv(),tsv()]))
        self.assertEqual(result[GENE]["original_strict_nonpal_SNPs"],1)

    def test_bad_QTL_AC_AN_flagged_instead_of_proxy_imputed(self):
        result,_=M.scan_one("QTD000131","fake",{GENE},
                     open_tabix=lambda _:FakeTabix([tsv(an="774",ac="900")]))
        self.assertEqual(result[GENE]["invalid_AN_AC_count"],1)
        self.assertEqual(result[GENE]["per_SNP_regression_SE_MAF_AN_residual_proxy_all"]["n"],0)

    def test_molecular_trait_source_mismatch_is_fatal(self):
        d=tsv(gene="ENSG00000000001")
        with self.assertRaises(ValueError):
            M.scan_one("QTD000131","fake",{"ENSG00000000001"},
                 open_tabix=lambda _:FakeTabix([d.replace("ENSG00000000001","ENSG00000065613",1)]))

    def test_original_36_group_sci_and_source_gates(self):
        m,u,a,w,c,o,s=source_inputs()
        keys=M.pin_evidence(m,u,a,w,c,o)
        self.assertEqual(len(keys),36)
        self.assertEqual(len({g for q,g in keys}),20)

    def test_final_72_outcome_paths_not_promoted_without_sdY(self):
        m,u,a,w,c,o,s=source_inputs()
        result=M.evaluate(m,u,a,w,c,o,s,{})
        self.assertEqual(result["n_real_gene_tissue_assays"],36)
        self.assertEqual(result["n_groups_with_variant_N_heterogeneity"],0)
        self.assertEqual(result["n_formal_coloc_abf_eligible"],0)
        self.assertEqual(sum(len(g["GWAS_ancestry_effect_size_route"]) for g in result["groups"]),72)
        self.assertTrue(all(not z["may_set_sdY_equal_1_without_per_gene_verification"]
                            for z in result["groups"]))

    def test_not_allow_synthetic_manuscript_count_if_real_source_swapped(self):
        m,u,a,w,c,o,s=source_inputs()
        s["QTD000131"][m["groups"][0]["gene_id"]]["original_strict_nonpal_SNPs"]=2
        with self.assertRaises(ValueError):M.evaluate(m,u,a,w,c,o,s,{})

    def test_case_control_route_never_fabricates_GWAS_N_eff(self):
        m,u,a,w,c,o,s=source_inputs()
        v=M.evaluate(m,u,a,w,c,o,s,{})
        self.assertFalse(v["GWAS_original_per_SNP_N_eff_fabricated"])
        self.assertTrue(v["proxy_not_direct_sdY_and_not_approved_input"])
        self.assertEqual(v["n_formal_coloc_susie_executed"],0)


if __name__=="__main__":
    unittest.main()
