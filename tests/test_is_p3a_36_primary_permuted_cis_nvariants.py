"""P3A 36-source primary n_variants vs full cis nominal source reconciliation traps."""
import gzip
import importlib
import io
import json
import unittest
from pathlib import Path
from unittest.mock import patch

M=importlib.import_module("scripts.is.audit_p3a_36_primary_permuted_cis_nvariants")
MANIFEST=Path("docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json")


def fake_input():
    m=json.loads(MANIFEST.read_text(encoding="utf-8"))
    keys=m["groups"]
    genes={g["gene_id"] for g in keys}
    refs={x:{"phenotype_pos_GRCh38":103100000,
             "predicted_cis_1Mb_left":102100000,
             "predicted_cis_1Mb_right":104100000} for x in genes}
    union={"schema":"is-p3a-36-expanded-original-cis-reference-union-v1",
           "technical_status":"PASS_THREE_REAL_INDEXED_REFERENCE_UNION_SURVEYS",
           "scientific_status":"HOLD","n_genuine_assay_groups":36,
           "n_approved_coloc_abf":0,
           "cis_union_reference_Ensembl105_GRCh38_1based":M.UNION.EXPECTED_WINDOW,
           "groups":[{"QTL_dataset":x["QTL_dataset"],"gene_id":x["gene_id"],
                      "extended_context_gene_association_rows":10}
                     for x in keys]}
    ref={"schema":"is-p3a-official-105-phenotype-cis-window-context-audit-v1",
         "technical_status":"PASS_REF_POS_20_GENES_AND_36_KEYS",
         "scientific_status":"HOLD","genes":refs,
         "groups":[{"QTL_dataset":x["QTL_dataset"],"gene_id":x["gene_id"]}
                   for x in keys]}
    off={"technical_status":"PASS_3_SOURCE_IDENTITIES_IN_IMMUTABLE_CATALOGUE",
         "QTL_datasets":{
             q:{"nominal_catalogue_sample_size":n,"source_https_probe":
                "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/"
                  +q+"/"+q+".all.tsv.gz"} for q,n in M.QTL_N.items()}}
    allele={"schema":"is-p3a-primary-gwas-yaml-QTL-allele-N-36-v1",
            "technical_status":"PASS_PRIMARY_PUBLISHER_METADATA_AND_36_QTL_AN_SOURCE_READ",
            "formal_coloc_abf_eligible":0,
            "groups":[{"QTL_dataset":x["QTL_dataset"],"gene_id":x["gene_id"],
                       "original_variant_allele_N_summary":
                        {"count_valid_nonpal_original_QTL":6}}
                       for x in keys]}
    perms={q:{x["gene_id"]:[{"n_variants":8,
                              "molecular_trait_id_original":x["gene_id"]}]
              for x in keys if x["QTL_dataset"]==q} for q in M.QTL_N}
    nominal={q:{x["gene_id"]:{
                   "nominal_gene_rows_physical":10,
                   "nominal_all_variant_count_unique":8,
                   "observed_molecular_trait_ids":[x["gene_id"]],
                   "observed_outside_reference_1Mb":0,
                   "source_variant_class_counts":{"SNP":7,"INDEL":1},
                   "duplicate_rsid_source_rows":2}
                for x in keys if x["QTL_dataset"]==q} for q in M.QTL_N}
    return m,union,ref,off,allele,perms,nominal


class Response:
    status=200
    def __init__(self,b,url):self.b=b;self.url=url
    def __enter__(self):return self
    def __exit__(self,*a):return False
    def read(self,maxbytes):return self.b[:maxbytes]
    def geturl(self):return self.url


def make_perm(gene="ENSG00000000001",count="9",header="molecular_trait_id\tn_variants"):
    b=header+"\n"+gene+"\t"+count+"\n"
    return gzip.compress(b.encode("utf-8"))


class OriginalCisPermutationTests(unittest.TestCase):
    def test_source_permutation_url_must_match_nominal_stem(self):
        m,u,r,o,a,p,n=fake_input()
        url=M.permutation_url("QTD000131",o)
        self.assertTrue(url.endswith("/QTD000131.permuted.tsv.gz"))

    def test_36_real_groups_20_genes_pinned(self):
        m,u,r,o,a,p,n=fake_input()
        keyed=M.pin_inputs(m,u,r,o,a)
        self.assertEqual(len(keyed),36)
        self.assertEqual(len({x[1] for x in keyed}),20)

    def test_all_permutation_n_variants_equal_nominal_source(self):
        m,u,r,o,a,p,n=fake_input()
        result=M.evaluate(m,u,r,o,a,p,n,{})
        self.assertEqual(result["n_gene_tissue_groups_matching_full_cis_tested_variant_COUNT"],36)
        self.assertEqual(result["n_groups_without_matching_source_n_variants"],0)
        self.assertEqual(result["n_formal_coloc_abf_eligible"],0)
        self.assertTrue(all(not g["formal_coloc_abf_eligible"] for g in result["groups"]))

    def test_nominal_n_variants_shortfall_not_approved(self):
        m,u,r,o,a,p,n=fake_input()
        x=next(iter(n["QTD000131"].values()))
        x["nominal_all_variant_count_unique"]=7
        result=M.evaluate(m,u,r,o,a,p,n,{})
        self.assertEqual(result["n_groups_without_matching_source_n_variants"],1)
        self.assertTrue(any(g["status"]=="HOLD_CIS_NVAR_OR_PHENOTYPE_MISMATCH"
                            for g in result["groups"]))

    def test_permutation_multi_trait_rows_stay_unresolved(self):
        m,u,r,o,a,p,n=fake_input()
        x=next(iter(p["QTD000131"].values()))
        x.append({"n_variants":8})
        result=M.evaluate(m,u,r,o,a,p,n,{})
        self.assertEqual(result["n_groups_without_matching_source_n_variants"],1)

    def test_nominal_molecular_trait_mismatch_stays_unresolved(self):
        m,u,r,o,a,p,n=fake_input()
        x=next(iter(n["QTD000131"].values()))
        x["observed_molecular_trait_ids"]=["ENSGWRONG"]
        result=M.evaluate(m,u,r,o,a,p,n,{})
        self.assertEqual(result["n_groups_without_matching_source_n_variants"],1)

    def test_source_row_change_is_blocking_not_dropped(self):
        m,u,r,o,a,p,n=fake_input()
        x=next(iter(n["QTD000131"].values()))
        x["nominal_gene_rows_physical"]-=1
        with self.assertRaises(ValueError):M.evaluate(m,u,r,o,a,p,n,{})

    def test_off_cis_source_row_is_blocking(self):
        m,u,r,o,a,p,n=fake_input()
        x=next(iter(n["QTD000131"].values()))
        x["observed_outside_reference_1Mb"]=1
        with self.assertRaises(ValueError):M.evaluate(m,u,r,o,a,p,n,{})

    def test_wrong_enrichment_of_snp_count_rejected(self):
        m,u,r,o,a,p,n=fake_input()
        x=next(iter(n["QTD000131"].values()))
        x["nominal_all_variant_count_unique"]=5
        with self.assertRaises(ValueError):M.evaluate(m,u,r,o,a,p,n,{})

    def test_permutation_source_gzip_and_total_tested_variants(self):
        m,u,r,o,a,p,n=fake_input()
        q="QTD000131";gene=next(iter(p[q]))
        perm=gzip.compress(("molecular_trait_id\tn_variants\tp_beta\n"+
                              gene+"\t8\t0.02\n").encode("utf-8"))
        url=M.permutation_url(q,o)
        records,prov=M.fetch_permuted(q,o,
                      fetch=lambda req,timeout:Response(perm,url))
        self.assertEqual(records[gene][0]["n_variants"],8)
        self.assertEqual(prov["number_original_permutation_rows"],1)
        self.assertEqual(len(prov["gzip_SHA256"]),64)

    def test_permutation_missing_n_variants_rejected(self):
        m,u,r,o,a,p,n=fake_input()
        q="QTD000131";gene=next(iter(p[q]))
        b=gzip.compress(("molecular_trait_id\tp_beta\n"+gene+"\t0.02\n").encode())
        with self.assertRaises(ValueError):
            M.fetch_permuted(q,o,fetch=lambda req,timeout:Response(b,M.permutation_url(q,o)))

    def test_non_numeric_permutation_tested_variant_count_rejected(self):
        m,u,r,o,a,p,n=fake_input()
        q="QTD000131";gene=next(iter(p[q]))
        b=gzip.compress(("molecular_trait_id\tn_variants\n"+gene+"\tNaN\n").encode())
        with self.assertRaises(ValueError):
            M.fetch_permuted(q,o,fetch=lambda req,timeout:Response(b,M.permutation_url(q,o)))

    def test_incorrect_publisher_host_redirect_rejected(self):
        m,u,r,o,a,p,n=fake_input()
        q="QTD000131";gene=next(iter(p[q]))
        b=gzip.compress(("molecular_trait_id\tn_variants\n"+gene+"\t8\n").encode())
        with self.assertRaises(ValueError):
            M.fetch_permuted(q,o,fetch=lambda req,timeout:Response(b,"https://fake.example/malware"))

    def test_source_variant_n_not_statistical_effective_n(self):
        self.assertEqual(M.PERMUTED_HEADER,{"molecular_trait_id","n_variants"})
        m,u,r,o,a,p,n=fake_input()
        result=M.evaluate(m,u,r,o,a,p,n,{})
        self.assertFalse(result["all_gene_tissue_model_effective_N_approved"])
        self.assertEqual(result["n_formal_coloc_susie_executed"],0)


if __name__=="__main__":
    unittest.main()
