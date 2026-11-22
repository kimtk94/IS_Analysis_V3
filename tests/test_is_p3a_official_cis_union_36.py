"""Synthetic failure traps for original QTL reference union scan; no network."""
import importlib
import json
import unittest
from pathlib import Path
from unittest.mock import patch

M=importlib.import_module("scripts.is.audit_p3a_official_cis_union_36")
MANIFEST=Path("docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json")


def fake_reference():
    manifest=json.loads(MANIFEST.read_text(encoding="utf8"))
    ids=sorted({g["gene_id"] for g in manifest["groups"]})
    points={gene:103452402 for gene in ids}
    points[ids[0]]=102837472
    points[ids[-1]]=104353833
    genes={}
    for key,pos in points.items():
        genes[key]={
          "gene_id":key,"chromosome":"10",
          "official_gene_reference_found":True,
          "phenotype_pos_GRCh38":pos,
          "predicted_cis_1Mb_left":pos-1000000,
          "predicted_cis_1Mb_right":pos+1000000,
          "reference_window_is_study_specific_QTL_nominal_cis_proof":False,
        }
    ref={"schema":"is-p3a-official-105-phenotype-cis-window-context-audit-v1",
       "technical_status":"PASS_REF_POS_20_GENES_AND_36_KEYS",
       "scientific_status":"HOLD","n_unique_ref_position_found":20,
       "n_real_groups":36,"n_group_coloc_abf_eligible":0,
       "n_reference_cis_entirely_in_prior_context":1,
       "full_gene_specific_nominal_tested_cis_proven":False,
       "reference_105_positions_cannot_approve_dataset_specific_TSS":True,
       "official_Ensembl_105_phenotype_reference":{
          "publisher_MD5":M.CIS.EXPECTED_MD5},
       "genes":genes,
       "groups":[{"region_id":x["region_id"],"gene_id":x["gene_id"],
         "QTL_dataset":x["QTL_dataset"]} for x in manifest["groups"]]}
    return manifest,ref


def source_fixture(gene):
    fields={x:"." for x in M.SOURCE_FIELDS}
    fields.update({"molecular_trait_id":gene,"chromosome":"10",
      "position":"103452402","ref":"A","alt":"C",
      "variant":"chr10_103452402_A_C",
      "ma_samples":"42","maf":"0.22","pvalue":"0.0001",
      "beta":"0.1","se":"0.04","type":"SNP",
      "ac":"50","an":"300","r2":"1",
      "molecular_trait_object_id":gene,"gene_id":gene,
      "median_tpm":"5","rsid":"rs123"})
    return "\t".join(fields[k] for k in M.SOURCE_FIELDS)


class FakeDNA:
    def __getitem__(self,i):return ord("A")


class FakeTabix:
    contigs=["10"]
    def __init__(self,records):self.records=records
    def __enter__(self):return self
    def __exit__(self,*a):return None
    def fetch(self,chrom,start,end):
        assert chrom=="10"
        assert start==M.EXPECTED_WINDOW["start"]-1
        assert end==M.EXPECTED_WINDOW["end"]
        return iter(self.records)


class SourceUnionTests(unittest.TestCase):
    def test_reference_union_exact_20_gene_bounds(self):
        m,r=fake_reference()
        self.assertEqual(M.bounds(r["genes"]),M.EXPECTED_WINDOW)
        self.assertEqual(len(r["genes"]),20)

    def test_extremity_lost_is_rejected(self):
        m,r=fake_reference()
        g=next(iter(r["genes"]))
        r["genes"].pop(g)
        with self.assertRaises(ValueError):M.bounds(r["genes"])

    def test_phenotype_start_mismatch_rejected(self):
        m,r=fake_reference()
        g=next(iter(r["genes"]))
        r["genes"][g]["predicted_cis_1Mb_left"]+=1
        with self.assertRaises(ValueError):M.bounds(r["genes"])

    def test_fake_source_cis_approval_rejected(self):
        m,r=fake_reference()
        g=next(iter(r["genes"]))
        r["genes"][g]["reference_window_is_study_specific_QTL_nominal_cis_proof"]=True
        with self.assertRaises(ValueError):M.bounds(r["genes"])

    def test_positional_source_observations_accumulate(self):
        m,r=fake_reference()
        gene=list(r["genes"])[1]
        pos=r["genes"][gene]["phenotype_pos_GRCh38"]
        summary={gene:{"source_rows":0,"observed_min_GRCh38_pos":None,
          "observed_max_GRCh38_pos":None,
          "observed_rows_outside_Ensembl105_reference_cis":0}}
        M.update_position_stats(pos,gene,r["genes"],summary)
        M.update_position_stats(pos-1000001,gene,r["genes"],summary)
        self.assertEqual(summary[gene]["source_rows"],2)
        self.assertEqual(summary[gene]["observed_rows_outside_Ensembl105_reference_cis"],1)
        self.assertEqual(summary[gene]["observed_min_GRCh38_pos"],pos-1000001)

    def test_indexed_read_tracks_original_beta_se_and_gene(self):
        m,r=fake_reference()
        gene=list(r["genes"])[1]
        row=source_fixture(gene)
        record=M.scan_one("https://example.invalid/test.all.tsv.gz",
           M.EXPECTED_WINDOW,FakeDNA(),{gene:r["genes"][gene]},
           tabix_factory=lambda url:FakeTabix([row]))
        self.assertEqual(record["interval_rows"]["region_rows_scanned"],1)
        self.assertEqual(record["gene_rows"][gene]["valid_nonpal_SNPs"],1)
        self.assertEqual(record["gene_study_reference_cis_checks"][gene]["source_rows"],1)
        self.assertEqual(record["source_nominal_header_from_pinned_release"],list(M.SOURCE_FIELDS))

    def test_truncated_original_tabix_row_rejected(self):
        m,r=fake_reference()
        gene=list(r["genes"])[1]
        with self.assertRaises(ValueError):
            M.scan_one("https://example.invalid/test",
              M.EXPECTED_WINDOW,FakeDNA(),{gene:r["genes"][gene]},
              tabix_factory=lambda url:FakeTabix(["10\t103452402"]))

    def test_no_chr10_remote_index_rejected(self):
        m,r=fake_reference()
        gene=list(r["genes"])[1]
        class NoIndex(FakeTabix):
            contigs=["8"]
        with self.assertRaises(ValueError):
            M.scan_one("https://example.invalid/test",M.EXPECTED_WINDOW,
              FakeDNA(),{gene:r["genes"][gene]},
              tabix_factory=lambda url:NoIndex([]))

    def test_join_36_holds_even_if_all_real_sources_pass(self):
        m,r=fake_reference()
        groups=m["groups"]
        baseline={"groups":[{"region_id":g["region_id"],
          "QTL_dataset":g["QTL_dataset"],"gene_id":g["gene_id"],
          "GRCh38_context_gene_associations_seen":5,
          "GRCh38_context_valid_nonpal_unique":4} for g in groups]}
        datasets={}
        for dataset in M.EXPECTED_GROUPS:
            genes={x["gene_id"] for x in groups if x["QTL_dataset"]==dataset}
            datasets[dataset]={
                "technical_status":"PASS_REAL_INDEXED_SOURCE",
                "gene_rows":{g:{"candidate_rows":9,"valid_nonpal_SNPs":8}
                             for g in genes},
                "gene_study_reference_cis_checks":{
                    g:{"source_rows":9,"observed_min_GRCh38_pos":103452402,
                       "observed_max_GRCh38_pos":103452402,
                       "observed_rows_outside_Ensembl105_reference_cis":0}
                    for g in genes}}
        result=M.join_results(m,baseline,r,M.EXPECTED_WINDOW,datasets)
        self.assertEqual(result["technical_status"],
             "PASS_THREE_REAL_INDEXED_REFERENCE_UNION_SURVEYS")
        self.assertEqual(result["n_genuine_assay_groups"],36)
        self.assertEqual(result["total_new_candidate_rows"],36*4)
        self.assertEqual(result["n_approved_coloc_abf"],0)
        self.assertEqual(result["n_executed_coloc_abf"],0)
        self.assertFalse(result["study_original_gene_specific_cis_tested_variants_proven"])
        self.assertTrue(all(not x["formal_coloc_abf_eligible"] for x in result["groups"]))

    def test_superset_loss_rejected(self):
        m,r=fake_reference()
        baseline={"groups":[{"region_id":g["region_id"],
           "QTL_dataset":g["QTL_dataset"],"gene_id":g["gene_id"],
           "GRCh38_context_gene_associations_seen":5,
           "GRCh38_context_valid_nonpal_unique":4} for g in m["groups"]]}
        datasets={}
        for dataset in M.EXPECTED_GROUPS:
            genes={x["gene_id"] for x in m["groups"] if x["QTL_dataset"]==dataset}
            datasets[dataset]={
              "technical_status":"PASS_REAL_INDEXED_SOURCE",
              "gene_rows":{g:{"candidate_rows":3,"valid_nonpal_SNPs":2} for g in genes},
              "gene_study_reference_cis_checks":{g:{"source_rows":3}
                                                   for g in genes}}
        with self.assertRaises(ValueError):
            M.join_results(m,baseline,r,M.EXPECTED_WINDOW,datasets)

    def test_reference_source_checksum_drift_rejected(self):
        m,r=fake_reference()
        r["official_Ensembl_105_phenotype_reference"]["publisher_MD5"]="fake"
        old={"chrom":"10","start":102403919,"end":104460593}
        with patch.object(M.BASE,"verify_inputs",return_value=old):
            with self.assertRaises(ValueError):
                M.check_inputs(m,{"groups":[]},r,{},{})

    def test_original_key_not_overwritten_by_other_QTL_tissue(self):
        m,r=fake_reference()
        keys={(g["region_id"],g["QTL_dataset"],g["gene_id"]) for g in m["groups"]}
        self.assertEqual(len(keys),36)
        self.assertEqual(len({k[2] for k in keys}),20)

    def test_union_contains_original_and_two_sided_extensions(self):
        reg=M.EXPECTED_WINDOW
        self.assertEqual(M.OLD_WINDOW["start"]-reg["start"],566447)
        self.assertEqual(reg["end"]-M.OLD_WINDOW["end"],893240)
        self.assertEqual(reg["end"]-reg["start"]+1,3516362)


if __name__=="__main__":
    unittest.main()
