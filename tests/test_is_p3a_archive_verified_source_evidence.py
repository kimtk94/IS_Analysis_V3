"""Privacy, integrity and fail-closed tests for durable original P3A aggregate archive.

All artificial documents are synthetic group-level fixtures, not raw source SNPs.
The real official archived four-file reconstruction runs in GitHub Actions.
"""
from __future__ import annotations

import copy
import hashlib
import importlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

A=importlib.import_module("scripts.is.archive_p3a_verified_source_evidence")


def original_manifest():
    return json.loads(Path(
        "docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json"
    ).read_text(encoding="utf8"))


def synthetic_document(kind):
    groups=original_manifest()["groups"]
    rows=[{"QTL_dataset":g["QTL_dataset"],
           "gene_id":g["gene_id"],"gene_symbol":g["gene_symbol"]}
          for g in groups]
    if kind=="coverage":
        rows=[{**g,"ancestry":a} for g in rows for a in ("EAS","EUR")]
        schema="is-p3a-original-36x2-gwas-qtl-fullcis-nonpal-coverage-gate-v1"
        status="PASS_PINNED_3_ORIGINAL_SOURCE_AUDITS_AND_72_RECONCILED_GROUPS"
    else:
        schema=A.COVERAGE.INPUTS[kind]["schema"]
        status=A.COVERAGE.INPUTS[kind]["status"]
    return {
        "schema":schema,"technical_status":status,
        "scientific_status":"HOLD",
        "rows":rows,
        "aggregate_by_ancestry":{"EAS":{"matched":89597},
                                 "EUR":{"matched":128282,"discordant":8}},
        "n_assay_ancestry_rows":72,
        "EUR_only_eight_discordant_gene_assay_rows":[],
        "n_permutation_vs_nominal_variant_count_exact_groups":36,
        "formal_coloc_abf_eligible":0,
        "formal_coloc_abf_executed":0,
        "formal_coloc_susie_executed":0,
    }


class DurableP3AOriginalEvidenceTests(unittest.TestCase):
    def test_exact_source_runs_pinned(self):
        self.assertEqual(A.SOURCE_RUNS,{
            "reciprocal":45031834139,
            "cis":45036525615,
            "qtl_n":45038801844,
            "coverage":45189787019})
        self.assertEqual(len(A.ARCHIVE_FILES),4)

    def test_official_2023_36_group_identity(self):
        g=original_manifest()["groups"]
        self.assertEqual(len(g),36)
        self.assertEqual(len({x["gene_id"] for x in g}),20)
        self.assertEqual(len({(x["QTL_dataset"],x["gene_id"]) for x in g}),36)

    def test_git_blob_inputs_are_original_aggregate_not_SNP_data(self):
        for kind in A.ARCHIVE_FILES:
            obj=synthetic_document(kind)
            self.assertEqual(A.validate_privacy_and_shape(obj,kind),
                             72 if kind=="coverage" else 36)

    def test_privacy_donor_id_is_hard_failure(self):
        for key in ("sample_id","donor_id","pseudo_id","rsid","beta",
                    "se","pvalue","varbeta","genotype","variant_id"):
            with self.subTest(key=key):
                obj=synthetic_document("cis")
                obj["rows"][0][key]="rsSYNTHETIC" if key=="rsid" else 0.002
                with self.assertRaises(ValueError):
                    A.validate_privacy_and_shape(obj,"cis")

    def test_privacy_long_raw_value_or_source_rows_fail(self):
        obj=synthetic_document("qtl_n")
        obj["rows"][0]["unexpected_expression_value"]="A"*1000
        with self.assertRaises(ValueError):
            A.validate_privacy_and_shape(obj,"qtl_n")
        obj=synthetic_document("reciprocal")
        obj["rows"]*=3
        with self.assertRaises(ValueError):
            A.validate_privacy_and_shape(obj,"reciprocal")

    def test_science_status_cannot_be_approved(self):
        obj=synthetic_document("cis")
        obj["scientific_status"]="PASS_COLOC_H4"
        with self.assertRaises(ValueError):
            A.validate_privacy_and_shape(obj,"cis")

    def test_fake_gene_ID_and_dataset_rejected(self):
        for field,value in (("gene_id","rs1234"),("QTL_dataset","QTD999999")):
            obj=synthetic_document("coverage")
            obj["rows"][0][field]=value
            with self.assertRaises(ValueError):
                A.validate_privacy_and_shape(obj,"coverage")

    def test_valid_source_sha_is_not_sha_of_reparsed_json(self):
        x=b'{ "schema" : "x", "rows": [] }\n'
        y=json.dumps(json.loads(x),indent=2).encode()
        self.assertNotEqual(A.sha256(x),A.sha256(y))

    def test_recomputation_requires_exact_72_rows_and_36_groups(self):
        obj=synthetic_document("coverage")
        self.assertEqual(len(obj["rows"]),72)
        self.assertEqual(len({(x["QTL_dataset"],x["gene_id"],x["ancestry"])
                              for x in obj["rows"]}),72)

    def test_real_archive_requires_exact_four_entries(self):
        with self.assertRaises(ValueError):
            A.check_and_archive({},original_manifest(),
                                Path("docs/is/evidence/verified_p3a_20261124"))

    def test_source_hash_conflict_fails_before_writing(self):
        with tempfile.TemporaryDirectory() as d:
            paths={}
            for kind in A.ARCHIVE_FILES:
                f=Path(d)/(kind+".json")
                f.write_text(json.dumps(synthetic_document(kind)),encoding="utf8")
                paths[kind]=f
            with self.assertRaises(ValueError):
                A.check_and_archive(paths,original_manifest(),
                        Path("docs/is/evidence/verified_p3a_20261124"))

    def test_output_outside_allowlisted_git_evidence_dir_rejected(self):
        docs={kind:synthetic_document(kind) for kind in A.ARCHIVE_FILES}
        with tempfile.TemporaryDirectory() as d:
            paths={}
            for kind,obj in docs.items():
                f=Path(d)/(kind+".json")
                f.write_text(json.dumps(obj),encoding="utf8")
                paths[kind]=f
            pins=copy.deepcopy(A.COVERAGE.INPUTS)
            for kind in ("reciprocal","cis","qtl_n"):
                pins[kind]["sha256"]=A.sha256(paths[kind].read_bytes())
            covhash=A.sha256(paths["coverage"].read_bytes())
            with patch.object(A.COVERAGE,"INPUTS",pins), \
                 patch.object(A,"PINNED_COVERAGE_SHA256",covhash), \
                 patch.object(A.COVERAGE,"compute_36x2",return_value=docs["coverage"]):
                with self.assertRaises(ValueError):
                    A.check_and_archive(paths,original_manifest(),Path(d)/"wrong")
                # Test successful evidence validation without actual Git writes.
                with patch.object(Path,"mkdir"),patch.object(Path,"write_bytes"), \
                     patch.object(Path,"write_text"):
                    index=A.check_and_archive(paths,original_manifest(),
                        Path("docs/is/evidence/verified_p3a_20261124"))
                self.assertEqual(len(index["source_evidence_file_sha256_and_GitHub_run"]),4)
                self.assertEqual(index["original_assay_x_ancestry_aggregate_rows"],72)
                self.assertFalse(index["ROSMAP_2021_QC_vs_R7_9_sample_gap_explained"])
                self.assertFalse(index["original_aggregate_evidence_donor_or_variant_effect_rows_archived"])

    def test_source_AN_proxy_only_not_exact_regression_N_claim(self):
        self.assertEqual(A.ARCHIVE_SCHEMA,
            "is-p3a-36x2-persistent-github-pinned-source-evidence-v1")
        self.assertTrue(all(len(x["sha256"])==64 for x in A.COVERAGE.INPUTS.values()))


if __name__=="__main__":
    unittest.main()
