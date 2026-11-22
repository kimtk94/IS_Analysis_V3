"""Offline negative-case tests, no original GWAS download or causal analyses."""
import copy
import importlib
import json
from pathlib import Path
import unittest

g=importlib.import_module("scripts.is.p3a_publisher_md5_gate")
j=importlib.import_module("scripts.is.p3a_join_original_source_evidence")

BASE=Path("docs/is/evidence/IS_CANDIDATE_WIDE_COLOC_BASELINE_20261122.json")

def samples():
    base=json.loads(BASE.read_text(encoding="utf-8"))
    today={}
    before={}
    for ancestry,sp in g.SOURCES.items():
        today[ancestry]={
            "schema":"is-p3a-publisher-original-byte-integrity-v1",
            "technical_status":"PASS_PUBLISHER_MD5_SHA256_SIZE",
            "ancestry":ancestry,"accession":sp["accession"],
            "file_size_bytes":sp["bytes"],
            "actual_file_md5":sp["publisher_md5_pinned"],
            "actual_file_sha256":sp["full_sha256_pinned"],
            "publisher_yaml_md5":sp["publisher_md5_pinned"],
            "full_gwas_file_verified_with_publisher_md5":True,
            "full_gwas_sample_size_approved":False,
            "archived_P3A_GWAS_allele_effect_verified":False,
            "scientific_status":"HOLD",
        }
        before[ancestry]={"research_ready":False,"results":[{
            "ancestry":ancestry,"accession":sp["accession"],
            "full_source_status":"TRANSFER_AND_STRUCTURE_QC_PASS",
            "total_variant_rows":sp["original_rows_in_prior_source_qc"],
            "bytes_downloaded":sp["bytes"],
            "sha256_calculated_local_not_publisher_verified":sp["full_sha256_pinned"],
            "per_variant_sample_size_column_present":False,
            "source_independent_sha256_verified":False
        }]}
    eqtl={
      "scientific_decision":"SOURCE_ROW_QC_PASS_COLOC_HOLD",
      "coloc_ABF_beta_varbeta_route":"BLOCKED_MISSING_VARBETA",
      "coloc_ABF_pvalue_MAF_N_route":"BLOCKED_MISSING_MAF_AND_CELLTYPE_N",
      "actual_per_celltype_N":"NOT_ESTABLISHED_BY_416_COHORT_ROWS",
      "source_nominal_build":"GRCh38"
    }
    supplement={
      "schema":"is-ci-public-source-gate-v1",
      "gigastroke_status":"PASS",
      "bbj_2021":{"status":"DISCOVERED","file_downloaded":False,"cohort_independent":False},
      "independent_replication":"NOT_ESTABLISHED"
    }
    return [base,today["EUR"],today["EAS"],before["EUR"],before["EAS"],eqtl,supplement]

class OriginalSourceAuditTests(unittest.TestCase):
    def test_official_GWAS_metadata_roundtrip(self):
        sp=g.SOURCES["EAS"]
        url,meta=g.build_urls(sp)
        self.assertTrue(url.endswith("GCST90104545_buildGRCh37.tsv.gz"))
        self.assertTrue(meta.endswith("-meta.yaml"))
        fields=g.publisher_metadata("\n".join([
            "gwas_id: GCST90104545",
            "genome_assembly: GRCh37",
            "data_file_name: GCST90104545_buildGRCh37.tsv.gz",
            "data_file_md5sum: 904e12e33ec3c2eeefe4dc93bec4995c",
        ]),sp)
        self.assertEqual(fields["data_file_md5sum"],sp["publisher_md5_pinned"])

    def test_wrong_publisher_md5_rejected(self):
        spec=g.SOURCES["EUR"]
        text="\n".join([
            "gwas_id: GCST90104540","genome_assembly: GRCh37",
            "data_file_name: GCST90104540_buildGRCh37.tsv.gz",
            "data_file_md5sum: "+"0"*32
        ])
        with self.assertRaisesRegex(ValueError,"MD5 drift"):
            g.publisher_metadata(text,spec)

    def test_build_mismatch_rejected(self):
        spec=g.SOURCES["EUR"]
        text="\n".join([
            "gwas_id: GCST90104540","genome_assembly: GRCh38",
            "data_file_name: GCST90104540_buildGRCh37.tsv.gz",
            "data_file_md5sum: "+spec["publisher_md5_pinned"]
        ])
        with self.assertRaisesRegex(ValueError,"accession/assembly"):
            g.publisher_metadata(text,spec)

    def test_md5_source_results_technically_pass_science_hold(self):
        values=samples()
        status=j.summarize(*values)
        self.assertEqual(status["technical_status"],
                         "PASS_ORIGINAL_GWAS_PUBLISHER_DIGEST_AND_SOURCE_EVIDENCE")
        self.assertEqual(status["scientific_status"],"HOLD")
        self.assertEqual(status["per_group_audit_rows"],0)
        self.assertEqual(status["unresolved_original_group_ID_count"],36)
        self.assertFalse(status["independent_BBJ_excluded_stroke_replications"])
        self.assertFalse(status["source_effect_allele_harmonization_approved"])
        self.assertEqual(status["new_coloc_abf_computations"],0)
        self.assertEqual(status["publisher_GWAS"]["EAS"]["original_rows_from_prior_full_QC"],6725216)

    def test_prior_SNP_count_drift_rejected(self):
        vals=samples()
        vals[4]["results"][0]["total_variant_rows"]=123
        with self.assertRaisesRegex(ValueError,"source identity"):
            j.summarize(*vals)

    def test_publisher_file_checksum_mismatch_rejected(self):
        vals=samples()
        vals[2]["actual_file_md5"]="0"*32
        with self.assertRaisesRegex(ValueError,"not a validated"):
            j.summarize(*vals)

    def test_false_sample_N_approval_rejected(self):
        vals=samples()
        vals[1]["full_gwas_sample_size_approved"]=True
        with self.assertRaisesRegex(ValueError,"not a validated"):
            j.summarize(*vals)

    def test_original_eQTL_missings_not_overridden(self):
        vals=samples()
        vals[5]["coloc_ABF_beta_varbeta_route"]="PASS"
        with self.assertRaisesRegex(ValueError,"eQTL"):
            j.summarize(*vals)

    def test_noncandidate_BBJ_release_not_independent_replication(self):
        vals=samples()
        vals[6]["bbj_2021"]["cohort_independent"]=True
        with self.assertRaisesRegex(ValueError,"independent replication"):
            j.summarize(*vals)

    def test_36_assay_group_ids_cannot_be_fabricated(self):
        vals=samples()
        vals[0]["assay_group_universe"]["abf_followup_source_flagged_groups"]=35
        with self.assertRaisesRegex(ValueError,"wrong legacy P3A"):
            j.summarize(*vals)

    def test_real_36_assays_recovered_with_20_unique_gene_ids(self):
        src=json.loads(Path("docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json").read_text())
        result=j.summarize(*samples(),restored_keys=src)
        self.assertEqual(result["per_group_audit_rows"],36)
        self.assertEqual(result["unresolved_original_group_ID_count"],0)
        self.assertEqual(result["unique_recovered_gene_ids"],20)
        self.assertTrue(result["restored_36_real_gene_QTL_keys"])
        self.assertFalse(result["server_original_positive_stat_group_assignment_verified"])
        self.assertEqual(result["server_positive_stat_group_membership_unresolved_count"],36)
        self.assertEqual(result["new_coloc_abf_computations"],0)
        self.assertTrue(all(not x["eligible"] for x in result["per_group_blockers_if_available"]))

    def test_recovered_manifest_false_gene_is_rejected(self):
        src=json.loads(Path("docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json").read_text())
        src["groups"][0]["gene_id"]="ENSG12345678901"
        src["groups"][1]["gene_id"]="ENSG12345678901"
        with self.assertRaisesRegex(ValueError,"restored actual"):
            j.summarize(*samples(),restored_keys=src)

    def test_source_positive_group_assignment_cannot_be_invented(self):
        src=json.loads(Path("docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json").read_text())
        src["groups"][0]["original_P3A_source_stat_positive_proven"]=True
        with self.assertRaisesRegex(ValueError,"restored actual"):
            j.summarize(*samples(),restored_keys=src)

if __name__=="__main__":
    unittest.main()
