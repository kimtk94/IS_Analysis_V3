"""Synthetic source fixtures for ROSMAP Release6/7/8 metadata/input inventory."""
import importlib
import json
import unittest
from pathlib import Path
from unittest.mock import patch

M = importlib.import_module("scripts.is.audit_p3a_rosmap_release_input_inventory")


def manifest():
    return json.loads(Path(
        "docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json"
    ).read_text(encoding="utf-8"))


def metadata(methods=("ge", "exon", "tx", "txrev", "leafcutter"), n=560,
             duplicate=False, study="QTS000025", group="brain_naive"):
    header = "study_id\tdataset_id\tstudy_label\tsample_group\tsample_size\tquant_method\n"
    rows = [f"{study}\tQTD{434+i:06}\tROSMAP\t{group}\t{n}\t{method}\n"
            for i, method in enumerate(methods)]
    if duplicate:
        rows.append(rows[0])
    return (header+"".join(rows)).encode()


def ftp(inputs=False, sample_n=560, wrong_path=False):
    fields = ["study_id", "dataset_id", "study_label", "sample_group",
              "sample_size", "quant_method", "ftp_path", "ftp_cs_path", "ftp_lbf_path"]
    if inputs:
        fields.append("sample_meta")
    base = "ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/"
    path = ("QTS000000/QTD000434/" if wrong_path else "QTS000025/QTD000434/")
    values = ["QTS000025","QTD000434","ROSMAP","brain_naive",str(sample_n),"ge",
              base+"sumstats/"+path+"QTD000434.all.tsv.gz",
              base+"susie/"+path+"QTD000434.credible_sets.tsv.gz",
              base+"susie/"+path+"QTD000434.lbf_variable.txt.gz"]
    if inputs:
        values.append("ftp://example.invalid/private/sample.tsv")
    return ("\t".join(fields)+"\n"+"\t".join(values)+"\n").encode()


def qtl():
    return {"input_requirements_proven_by_official_workflow":{"sample_meta":True},
            "original_QTD000434_normalized_matrix_recovered":False,
            "original_QTD000434_sample_meta_recovered":False}


def valid_release():
    return {v:M.release_audit(v,metadata()) for v in M.RELEASE_PATHS}


class OriginalPublisherReleaseInputAuditTests(unittest.TestCase):
    def test_R6_R7_R8_exact_560_and_study_id(self):
        for r in M.RELEASE_PATHS:
            got=M.release_audit(r,metadata())
            self.assertEqual(got["gene_expression_nominal_tissue_N"],560)
            self.assertEqual(got["gene_expression_dataset_id"],"QTD000434")
            self.assertFalse(got["qtlmap_input_sample_intersection_measured"])

    def test_r8_new_majiq_still_560(self):
        raw=metadata(("ge","exon","tx","txrev","leafcutter","majiq"))
        got=M.release_audit("R8",raw)
        self.assertEqual(got["ROSMAP_study_quantification_count"],6)

    def test_wrong_N_each_release_fails(self):
        for r in M.RELEASE_PATHS:
            with self.assertRaises(ValueError):
                M.release_audit(r,metadata(n=569))

    def test_missing_ge_QTD000434_fails(self):
        with self.assertRaises(ValueError):
            M.release_audit("R7",metadata(("exon","tx","txrev","leafcutter")))

    def test_duplicate_dataset_fails_closed(self):
        with self.assertRaises(ValueError):
            M.release_audit("R6",metadata(duplicate=True))

    def test_wrong_study_identity_fails_closed(self):
        with self.assertRaises(ValueError):
            M.release_audit("R7",metadata(study="QTS999999"))

    def test_wrong_qtl_group_fails_closed(self):
        with self.assertRaises(ValueError):
            M.release_audit("R6",metadata(group="wrong_tissue"))

    def test_official_manifest_is_output_inventory_not_R7_input_access_assertion(self):
        got=M.ftp_audit(ftp())
        self.assertTrue(got["official_manifest_contains_full_cis_sumstats_path"])
        self.assertFalse(got["has_direct_qtlmap_count_matrix_column"])
        self.assertFalse(got["has_direct_qtlmap_sample_meta_column"])
        self.assertFalse(got["Release7_normalized_expression_public_availability_determined"])
        self.assertFalse(got["file_absence_outside_manifest_claimed"])

    def test_unexpected_sample_meta_source_column_is_not_silently_ignored(self):
        with self.assertRaises(ValueError):
            M.ftp_audit(ftp(inputs=True))

    def test_ftp_wrong_release_study_N_rejected(self):
        with self.assertRaises(ValueError):
            M.ftp_audit(ftp(sample_n=569))

    def test_ftp_path_unsupported_study_rejected(self):
        with self.assertRaises(ValueError):
            M.ftp_audit(ftp(wrong_path=True))

    def test_parser_checks_actual_qtlmap_phenotype_sample_intersection_contract(self):
        fake={"source_status":"PASS_OFFICIAL_QTLMAP_NORMALIZED_MATRIX_AND_SAMPLE_INTERSECTION"}
        with patch.object(M.SMALL.PRIOR,"parse_qtlmap",return_value=fake) as p:
            parsed=M.qtlmap_model_audit(b"real qtlmap docs mocked by checksum source test")
        p.assert_called_once()
        self.assertFalse(parsed["original_QTD000434_normalized_matrix_recovered"])

    def test_parser_rejects_mutated_qtlmap_semantics(self):
        fake={"source_status":"HOLD_OFFICIAL_QTLMAP_CHANGED"}
        with patch.object(M.SMALL.PRIOR,"parse_qtlmap",return_value=fake):
            with self.assertRaises(ValueError):
                M.qtlmap_model_audit(b"wrong")

    def test_all_36_hold_and_n_560_not_new_release7_filter(self):
        a=M.assemble(manifest(),valid_release(),M.ftp_audit(ftp()),qtl(),{})
        self.assertTrue(a["560_sample_N_already_in_R6_not_newly_R7"])
        self.assertFalse(a["2021_Zenodo_QC_pass_569_to_560_exclusion_reason_proven"])
        self.assertFalse(a["R6_R7_R8_sample_overlap_independently_verified"])
        self.assertEqual(len(a["groups"]),36)
        self.assertEqual(a["n_original_ROSMAP_groups"],15)
        self.assertTrue(all(not x["coloc_eligible"] for x in a["groups"]))
        self.assertEqual(a["n_gene_specific_original_Release7_sdY_measured"],0)
        self.assertEqual(a["n_coloc_abf_executed"],0)

    def test_no_real_donor_identifiers_or_phenotype_values_emitted(self):
        s=json.dumps(M.assemble(manifest(),valid_release(),M.ftp_audit(ftp()),qtl(),{}))
        self.assertNotIn("PRIVATE_DONOR_",s)
        self.assertNotIn('"sdY": 1',s)
        self.assertNotIn('"estimated_N": 560',s)

    def test_manifest_mismatch_halts(self):
        m=manifest(); m["groups"]=m["groups"][:-1]
        with self.assertRaises(ValueError):
            M.assemble(m,valid_release(),M.ftp_audit(ftp()),qtl(),{})

    def test_absent_publisher_schema_halts(self):
        with self.assertRaises(ValueError):
            M.ftp_audit(b"dataset_id\tstudy_id\nQTD000434\tQTS000025\n")


if __name__ == "__main__":
    unittest.main()
