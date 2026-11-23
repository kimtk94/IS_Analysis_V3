"""Fail-closed tests for coloc original R beta/SE cc without case-fraction gate."""
import hashlib
import importlib
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

M=importlib.import_module("scripts.is.audit_p3a_coloc_original_cc_beta_no_s_gate")


def groups():
    return json.loads(Path(
        "docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json"
    ).read_text(encoding="utf-8"))


def passing_R_stdout():
    return "".join(f"ORIGINAL_COLOC_TEST_{name}=PASS\n" for name in M.EXPECTED_KEYS)


class OriginalSourceCaseControlInputGateTests(unittest.TestCase):
    def test_genuine_36_group_manifest_is_not_reduced(self):
        g=groups()
        self.assertEqual(len(g["groups"]),36)
        self.assertEqual(len({r["gene_id"] for r in g["groups"]}),20)
        self.assertEqual(sum(r["QTL_dataset"]=="QTD000434" for r in g["groups"]),15)

    def test_original_source_commit_is_pinned_not_main(self):
        self.assertEqual(len(M.UPSTREAM_COMMIT),40)
        self.assertRegex(M.UPSTREAM_COMMIT,r"^[a-f0-9]{40}$")
        self.assertNotIn("master",M.UPSTREAM_RAW_BASE)
        self.assertNotIn("/main/",M.UPSTREAM_RAW_BASE)

    def test_git_blob_sha1_has_real_git_object_prefix(self):
        b=b"hello\n"
        expected=hashlib.sha1(b"blob 6\x00hello\n",usedforsecurity=False).hexdigest()
        self.assertEqual(M.git_blob_sha1(b),expected)
        self.assertNotEqual(M.git_blob_sha1(b),hashlib.sha1(
            b,usedforsecurity=False).hexdigest())

    def test_original_source_payload_modified_rejected(self):
        for path in M.SOURCE_BLOB_SHA1:
            with self.subTest(path=path),self.assertRaises(ValueError):
                M.verify_source_code("check_dataset <- function(d) {TRUE}",path)

    def test_unknown_source_path_rejected(self):
        with self.assertRaises(ValueError):
            M.verify_source_code("anything", "R/random.R")

    def test_runtime_synthetic_success_has_all_six_keys(self):
        d=M.parse_runtime_check(passing_R_stdout())
        self.assertEqual(len(d),6)
        self.assertTrue(all(d.values()))

    def test_runtime_missing_test_fails_closed(self):
        fewer="".join(
            f"ORIGINAL_COLOC_TEST_{name}=PASS\n"
            for name in M.EXPECTED_KEYS[:-1])
        with self.assertRaises(ValueError):
            M.parse_runtime_check(fewer)

    def test_runtime_fail_test_fails_closed(self):
        with self.assertRaises(ValueError):
            M.parse_runtime_check(passing_R_stdout().replace(
                "=PASS","=FAIL",1))

    def test_runtime_unknown_test_rejected(self):
        with self.assertRaises(ValueError):
            M.parse_runtime_check(passing_R_stdout()+
                                  "ORIGINAL_COLOC_TEST_backdoor=PASS\n")

    def test_runtime_duplicate_rejected(self):
        with self.assertRaises(ValueError):
            M.parse_runtime_check(passing_R_stdout()+
                                  f"ORIGINAL_COLOC_TEST_{M.EXPECTED_KEYS[0]}=PASS\n")

    def test_runtime_empty_result_rejected(self):
        with self.assertRaises(ValueError):
            M.parse_runtime_check("")

    def test_wrong_original_group_count_blocks(self):
        with tempfile.TemporaryDirectory() as d:
            g=groups()
            g["groups"]=g["groups"][:-1]
            with self.assertRaises(ValueError):
                M.audit_original_contract(g,Path(d))

    def test_positive_audit_never_approves_formal_coloc_or_fabricates_s(self):
        with tempfile.TemporaryDirectory() as d:
            fake_files=lambda name: (f"fake original {name}\n".encode(),
                {"upstream_file":name,"git_blob_sha1":M.SOURCE_BLOB_SHA1[name]})
            with patch.object(M,"fetch_pinned_r_source",side_effect=fake_files), \
                 patch.object(M,"TEST_R",Path(
                    "tests/is/test_coloc_cc_beta_no_s_original_source.R")), \
                 patch.object(M.subprocess,"run",return_value=SimpleNamespace(
                    returncode=0,stdout=passing_R_stdout(),stderr="")) as pr:
                result=M.audit_original_contract(groups(),Path(d))
            self.assertEqual(result["n_original_coloc_R_checks_passed"],6)
            self.assertFalse(result[
                "direct_beta_plus_varbeta_case_fraction_s_required_by_original_coloc"])
            self.assertTrue(result[
                "pvalue_MAF_N_only_case_fraction_s_required_by_original_coloc"])
            self.assertFalse(result["original_source_QTL_gene_regression_N_validated"])
            self.assertFalse(result["QTL_estimated_sdY_approved_for_formal_coloc"])
            self.assertEqual(result["n_current_formal_coloc_executed"],0)
            self.assertTrue(pr.called)
            self.assertEqual(pr.call_args.kwargs["timeout"],120)

    def test_R_source_runtime_nonzero_exit_blocks(self):
        with tempfile.TemporaryDirectory() as d:
            with patch.object(M,"fetch_pinned_r_source",side_effect=lambda name:
                (b"fake\n",{"upstream_file":name})), \
                 patch.object(M.subprocess,"run",return_value=SimpleNamespace(
                    returncode=2,stdout="",stderr="error")):
                with self.assertRaises(ValueError):
                    M.audit_original_contract(groups(),Path(d))

    def test_original_check_no_s_in_GWAS_beta_path_not_same_as_quant_N(self):
        # Enforces both branches remain explicitly distinguished in result.
        keys=M.EXPECTED_KEYS
        self.assertIn("cc_beta_varbeta_without_s",keys)
        self.assertIn("cc_pvalue_only_missing_s_rejected",keys)
        self.assertIn("quant_beta_varbeta_missing_sdY_N_MAF_rejected",keys)

    def test_no_raw_GWAS_or_QTL_effects_logged_in_report(self):
        with tempfile.TemporaryDirectory() as d:
            with patch.object(M,"fetch_pinned_r_source",side_effect=lambda name:
                (b"fake\n",{"upstream_file":name})), \
                 patch.object(M.subprocess,"run",return_value=SimpleNamespace(
                    returncode=0,stdout=passing_R_stdout(),stderr="")):
                r=M.audit_original_contract(groups(),Path(d))
        out=json.dumps(r)
        self.assertNotIn("rs123",out)
        self.assertNotIn("variant_beta",out)
        self.assertFalse(r["cross_build_GRCh37_GRCh38_complete_QTL_GWAS_cis_harmonisation_approved"])


if __name__=="__main__":
    unittest.main()
