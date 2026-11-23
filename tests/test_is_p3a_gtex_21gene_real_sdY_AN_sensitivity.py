"""Offline adversarial tests for GTEx 21 original source QTL sdY and N QC."""
import copy
import hashlib
import importlib
import json
import math
import unittest
from pathlib import Path
from unittest.mock import patch

M=importlib.import_module("scripts.is.audit_p3a_gtex_21gene_real_sdY_AN_sensitivity")
C=M.QTL.OFFICIAL_COLUMN_ORDER


def source_manifest():
    return json.loads(Path(
        "docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json"
    ).read_text(encoding="utf-8"))


def make_lines(dataset="QTD000131",n_rows_per_gene=4,
               invalid_an=False,duplicate=False,wrong_beta_dup=False,
               mismatched_window=False):
    names=list(M.DATASETS[dataset]["gene_variant_counts"])
    n=M.DATASETS[dataset]["catalogue_n_hypothesis"]
    rows=[]
    for k,gene in enumerate(names):
        for i in range(n_rows_per_gene):
            pos=(102000000 + k*100 + i*2 if not mismatched_window
                 else 98000000 + k*100 + i*2)
            maf=[.03,.07,.22,.34][i % 4]
            AN=2*n if i!=1 else 2*(n-5)
            if invalid_an and k==0 and i==0:
                AN=-1
            sd=0.9
            se=math.sqrt(sd*sd/(2*n*maf*(1-maf)))
            d={
                "molecular_trait_id":gene,"chromosome":"10","position":str(pos),
                "ref":"A","alt":"C","variant":f"chr10_{pos}_A_C",
                "ma_samples":"4","maf":str(maf),"pvalue":"0.5",
                "beta":"0.03","se":str(se),"type":"SNP","ac":"50",
                "an":str(AN),"r2":"0.97",
                "molecular_trait_object_id":gene,
                "gene_id":gene,"median_tpm":"2","rsid":f"rs{i}",
            }
            line="\t".join(d[h] for h in C)+"\n"
            rows.append(line)
            if duplicate and k==0 and i==0:
                d["beta"]="0.04" if wrong_beta_dup else "0.03"
                rows.append("\t".join(d[h] for h in C)+"\n")
    return rows


def mocked_scan(dataset,lines):
    counts={gene:4 for gene in M.DATASETS[dataset]["gene_variant_counts"]}
    local=copy.deepcopy(M.DATASETS)
    local[dataset]["gene_variant_counts"]=counts
    local[dataset]["expected_nonpal_source_assay_records"]=sum(counts.values())
    local[dataset]["region_source_rows"]=len(lines)
    local[dataset]["region_stream_SHA256"]=hashlib.sha256(
        ("\t".join(C)+"\n"+"".join(lines)).encode()).hexdigest()
    with patch.object(M,"DATASETS",local), \
         patch.object(M,"MIN_GENE_VARIANTS",3), \
         patch.object(M,"MIN_FULL_AN_VARIANTS",2), \
         patch.object(M,"MIN_MAF_STRATUM",2):
        return M.scan_source(lines,dataset,source_manifest())


class GenuineGTExSDYSourceTests(unittest.TestCase):
    def test_original_21_manifest_and_variant_counts(self):
        g=source_manifest()
        self.assertEqual(len(g["groups"]),36)
        self.assertEqual(sum(x["QTL_dataset"]=="QTD000131" for x in g["groups"]),10)
        self.assertEqual(sum(x["QTL_dataset"]=="QTD000171" for x in g["groups"]),11)
        self.assertEqual(sum(M.DATASETS["QTD000131"]["gene_variant_counts"].values()),44693)
        self.assertEqual(sum(M.DATASETS["QTD000171"]["gene_variant_counts"].values()),45570)

    def test_gtex_aorta_all_gene_sdY_OLS_estimates(self):
        a=mocked_scan("QTD000131",make_lines())
        self.assertEqual(a["source_genuine_candidate_group_count"],10)
        self.assertEqual(a["source_gene_assay_nonpal_SNP_total"],40)
        self.assertEqual(a["source_variant_AN_heterogeneous_gene_assays"],10)
        for r in a["groups"]:
            self.assertAlmostEqual(
                r["original_coloc_sdY_est_at_nominal_study_N_UNAPPROVED"][
                    "sdY_est_coloc_zero_intercept_OLS"],.9,places=10)
            self.assertFalse(r["published_QTL_model_gene_N_directly_verified"])
            self.assertFalse(r["formal_coloc_abf_eligible"])
            self.assertEqual(r["source_an_over_2_modal_study_N_records"],3)

    def test_gtex_cortex_all_gene_sdY_OLS_estimates(self):
        a=mocked_scan("QTD000171",make_lines("QTD000171"))
        self.assertEqual(a["source_genuine_candidate_group_count"],11)
        self.assertEqual(a["source_gene_assay_nonpal_SNP_total"],44)
        self.assertEqual(a["n_QTL_estimates_only_UNAPPROVED"],11)
        self.assertEqual(a["source_variant_AN_heterogeneous_gene_assays"],11)

    def test_nominal_N_from_AN_is_only_hypothesis(self):
        a=mocked_scan("QTD000131",make_lines())
        for r in a["groups"]:
            self.assertEqual(r["source_an_over_2_hypothesis_only"],387)
            self.assertLess(r["source_an_over_2_min"],387)
            self.assertTrue(r["source_variant_heterogeneous_N_not_accepted_as_model_N"])
        self.assertEqual(a["n_QTL_N_regression_approved"],0)

    def test_no_raw_original_variant_rows_persisted(self):
        a=mocked_scan("QTD000131",make_lines())
        s=json.dumps(a)
        self.assertNotIn("rs0",s)
        self.assertNotIn("chr10_102000000_A_C",s)
        self.assertNotIn('"beta":',s)
        self.assertFalse(a["raw_source_QTL_effects_or_variant_IDs_exported"])

    def test_duplicate_same_source_variant_is_deduped(self):
        a=mocked_scan("QTD000131",make_lines(duplicate=True))
        self.assertEqual(a["source_gene_assay_nonpal_SNP_total"],40)
        self.assertEqual(a["source_all_region_original_rows"],41)

    def test_duplicate_with_conflicting_beta_fails(self):
        with self.assertRaises(ValueError):
            mocked_scan("QTD000131",make_lines(
                duplicate=True,wrong_beta_dup=True))

    def test_negative_AN_source_row_fails(self):
        with self.assertRaises(ValueError):
            mocked_scan("QTD000131",make_lines(invalid_an=True))

    def test_original_approved_region_exact(self):
        r=M.original_region()
        self.assertEqual(r["region"],"10:101837472-105353833")

    def test_original_source_unrelated_chr10_window_drift(self):
        with self.assertRaises(ValueError):
            mocked_scan("QTD000131",make_lines(mismatched_window=True))

    def test_original_source_stream_sha_mismatch_blocks(self):
        dataset="QTD000131"
        lines=make_lines(dataset)
        local=copy.deepcopy(M.DATASETS)
        local[dataset]["gene_variant_counts"]={k:4 for k in local[dataset]["gene_variant_counts"]}
        local[dataset]["expected_nonpal_source_assay_records"]=40
        local[dataset]["region_source_rows"]=40
        local[dataset]["region_stream_SHA256"]="0"*64
        with patch.object(M,"DATASETS",local), \
             patch.object(M,"MIN_GENE_VARIANTS",3), \
             patch.object(M,"MIN_FULL_AN_VARIANTS",2):
            with self.assertRaises(ValueError):
                M.scan_source(lines,dataset,source_manifest())

    def test_noncanonical_QTL_study_source_url_blocked(self):
        import urllib.request
        class NoNetwork:
            pass
        with self.assertRaises(ValueError):
            M.ORIGINAL.fetch_original_index(
                "https://evil.example/source/QTD000131.all.tsv.gz",opener=NoNetwork())
        with self.assertRaises(ValueError):
            M.ORIGINAL.fetch_original_index(
                "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/"
                "QTD999999/QTD999999.all.tsv.gz",opener=NoNetwork())

    def test_scalar_N_sensitivity_sqrt_ratio(self):
        x=[(.03,.15),(.22,.10),(.32,.09),(.08,.14)]
        f=M.OLS.estimate_coloc_sdY(x,387)["sdY_est_coloc_zero_intercept_OLS"]
        lo=M.OLS.estimate_coloc_sdY(x,366)["sdY_est_coloc_zero_intercept_OLS"]
        self.assertAlmostEqual(lo/f,math.sqrt(366/387),places=12)

    def test_odd_even_MAF_stratum_are_not_LD_independent(self):
        a=mocked_scan("QTD000131",make_lines())
        self.assertTrue(all(v["sorted_odd_even_NOT_independent_due_to_LD"]
                            for v in a["groups"]))
        self.assertTrue(all(v["MAF_low_high_ratio"] is not None for v in a["groups"]))

    def test_study_alias_not_mistaken_for_N_phenotype_count(self):
        self.assertEqual(M.DATASETS["QTD000131"]["catalogue_n_hypothesis"],387)
        self.assertEqual(M.DATASETS["QTD000171"]["catalogue_n_hypothesis"],205)
        self.assertTrue(M.DATASETS["QTD000131"]["region_stream_SHA256"])
        self.assertTrue(M.DATASETS["QTD000171"]["region_stream_SHA256"])

    def test_empty_original_source_fails_closed(self):
        with self.assertRaises(ValueError):
            mocked_scan("QTD000131",[])


if __name__=="__main__":
    unittest.main()
