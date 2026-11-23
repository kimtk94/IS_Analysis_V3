"""Original ROSMAP QC-report stages are not proof of Release7 sample identity."""
import importlib
import json
import unittest
from pathlib import Path

M = importlib.import_module("scripts.is.audit_p3a_rosmap_official_qc_report_stages")


def report_doc(before=581, after=576, matching=576, nonmatch=0, duplicate=False):
    row = lambda k, n: f"<pre>## {k}: {n}</pre>"
    s = "<h1>ROSMAP QC report</h1><h2>Concordance between expression and genotype data</h2>"
    for k, n in [
        ("Number of samples in study before quality control", before),
        ("Number of samples in study after quality control", after),
        ("Number of genotyped individuals in the current study ROSMAP", 581),
        ("Number of non-matching genotypes", nonmatch),
        ("Number of matching genotypes", matching),
    ]:
        s += row(k,n)
    if duplicate:
        s += row("Number of samples in study after quality control",after)
    s += "<script>Number of matching genotypes: 9999</script>"
    return s.encode()


def real_manifest():
    return json.loads(Path(
        "docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json"
    ).read_text(encoding="utf-8"))


def original2021():
    return {"2021_original_Zenodo_RNA_genotype_QC_pass_N": 569,
            "metadata_publisher_MD5": "038f2fc144129d401144f3924109e1df"}


def release7():
    return {"dataset_id":"QTD000434", "study":"ROSMAP",
            "sample_group":"brain_naive", "Release7_declared_tissue_N":560}


class ROSMAP_QC_Stage_EvidenceTests(unittest.TestCase):
    def test_real_public_report_expected_counts_and_no_script_contamination(self):
        r = M.parse_original_qc_report(report_doc())
        self.assertEqual(r["counts_only"], {
            "qc_before":581, "qc_after":576, "genotype_sample_count":581,
            "matched_genotypes":576, "nonmatching_genotypes":0})
        self.assertFalse(r["report_stage_source_release_certified_same_as_Release7"])

    def test_36_assay_science_status_and_counts_are_not_sample_exclusions(self):
        r = M.assemble(real_manifest(), M.parse_original_qc_report(report_doc()),
                       original2021(), release7(), "qcsha", "r7sha")
        self.assertEqual(len(r["groups"]), 36)
        self.assertEqual(r["n_ROSMAP_QTD000434_groups"], 15)
        self.assertEqual(r["cross_source_descriptive_only"][
            "official_report_after_QC_minus_Release7_nominal_N"], 16)
        self.assertEqual(r["cross_source_descriptive_only"][
            "2021_metadata_both_QC_minus_Release7_nominal_N"],9)
        self.assertEqual(r["cross_source_descriptive_only"][
            "official_report_after_QC_minus_2021_metadata_both_QC"],7)
        self.assertTrue(r["cross_source_descriptive_only"][
            "these_differences_are_not_verified_sequential_sample_exclusions"])
        self.assertEqual(r["n_coloc_abf_eligible"], 0)
        self.assertTrue(all(not x["scientifically_eligible_for_coloc"] for x in r["groups"]))

    def test_failed_publisher_identity_rejected(self):
        with self.assertRaises(ValueError):
            M.parse_original_qc_report(report_doc().replace(b"ROSMAP QC report",b"OTHER report"))

    def test_missing_publisher_field_rejected(self):
        raw = report_doc().replace(b"Number of genotyped individuals in the current study ROSMAP",
                                   b"Unknown field")
        with self.assertRaises(ValueError):
            M.parse_original_qc_report(raw)

    def test_changed_publisher_n_rejected(self):
        with self.assertRaises(ValueError):
            M.parse_original_qc_report(report_doc(after=575))

    def test_duplicate_public_field_rejected(self):
        with self.assertRaises(ValueError):
            M.parse_original_qc_report(report_doc(duplicate=True))

    def test_matched_concordance_inconsistent_rejected(self):
        with self.assertRaises(ValueError):
            M.parse_original_qc_report(report_doc(matching=575, nonmatch=1))

    def test_wrong_RELEASE7_n_rejected(self):
        bad=release7();bad["Release7_declared_tissue_N"]=569
        with self.assertRaises(ValueError):
            M.assemble(real_manifest(), M.parse_original_qc_report(report_doc()),
                       original2021(), bad, "qcsha", "r7sha")

    def test_wrong_original_2021_metadata_rejected(self):
        bad=original2021();bad["2021_original_Zenodo_RNA_genotype_QC_pass_N"]=576
        with self.assertRaises(ValueError):
            M.assemble(real_manifest(), M.parse_original_qc_report(report_doc()),
                       bad, release7(), "qcsha", "r7sha")

    def test_json_has_no_individual_ids_or_made_up_sdy(self):
        r=M.assemble(real_manifest(),M.parse_original_qc_report(report_doc()),
                     original2021(),release7(),"sha","sha")
        public=json.dumps(r)
        self.assertNotIn("PRIVATE_SAMPLE_",public)
        self.assertNotIn('"sdY": 1',public)
        self.assertFalse(r["original_sample_IDs_in_results"])
        self.assertFalse(r["cross_source_descriptive_only"][
            "ROSMAP_15_gene_expression_sdY_measured"])

    def test_missing_or_duplicate_groups_rejected(self):
        m=real_manifest()
        m["groups"]=m["groups"][:-1]
        with self.assertRaises(ValueError):
            M.assemble(m, M.parse_original_qc_report(report_doc()),
                       original2021(), release7(), "sha", "sha")

    def test_only_public_source_counts_not_sequential_lineage(self):
        r = M.parse_original_qc_report(report_doc())
        self.assertNotIn("donor_id",json.dumps(r))
        self.assertFalse(r["source_individual_identifiers_exported"])


if __name__=="__main__":
    unittest.main()
