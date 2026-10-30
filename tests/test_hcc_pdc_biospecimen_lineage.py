"""Unit tests: PDC case/sample/aliquot lineage, no patient ID leakage."""
import copy
import pathlib
import sys
import unittest
from unittest.mock import patch
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
import probe_pdc_biospecimen_lineage as module

SOURCE=[
 {"case_id":"C1","case_submitter_id":"PAT-1","samples":[
     {"sample_id":"S1","sample_submitter_id":"PAT-1-T","sample_type":"Primary Tumor",
      "aliquots":[{"aliquot_id":"A1","aliquot_submitter_id":"PAT-1-T-A",
                   "aliquot_run_metadata":[{"label":"tmt_126","alias":"batch1"}]}]},
     {"sample_id":"S2","sample_submitter_id":"PAT-1-N","sample_type":"Solid Tissue Normal",
      "aliquots":[{"aliquot_id":"A2","aliquot_submitter_id":"PAT-1-N-A"}]}]},
 {"case_id":"C2","case_submitter_id":"PAT-2","samples":[
     {"sample_id":"S3","sample_submitter_id":"PAT-2-T","sample_type":"Primary Tumor",
      "aliquots":[{"aliquot_id":"A3","aliquot_submitter_id":"PAT-2-T-A"}]}]},
]

class LineageTests(unittest.TestCase):
    def test_lines_from_tumor_and_normal(self):
        links=module.flatten(copy.deepcopy(SOURCE))
        r=module.count_summary(SOURCE,links)
        self.assertEqual(r["distinct_cases"],2)
        self.assertEqual(r["distinct_samples"],3)
        self.assertEqual(r["distinct_aliquots"],3)
        self.assertEqual(r["sample_types"]["Primary Tumor"]["unique_cases"],2)
        self.assertEqual(r["sample_types"]["Solid Tissue Normal"]["unique_cases"],1)
        self.assertNotIn("PAT-1",str(r))
        self.assertEqual(r["cases_with_tumor_and_normal_sample"],1)
        self.assertEqual(r["cases_without_aliquot"],0)
        self.assertEqual(r["run_metadata_records"],1)
        self.assertEqual(r["distinct_experiment_label_pairs"],0)


    def test_paginated_structure_and_counts(self):
        with patch.object(module,"request_graphql",
            return_value={"paginatedCasesSamplesAliquots":{"total":2,"casesSamplesAliquots":SOURCE}}):
            cases,meta=module.get_cases("PDC000198")
        self.assertEqual(len(cases),2)
        self.assertEqual(meta["repository_reported_cases"],2)
        self.assertNotIn("PAT-1",str(meta))

    def test_unauditable_incomplete_pagination(self):
        with patch.object(module,"request_graphql",
            return_value={"paginatedCasesSamplesAliquots":{"total":999,"casesSamplesAliquots":[]}}):
            with self.assertRaises(ValueError):
                module.get_cases("PDC000198")

    def test_header_prefix_joins_are_not_rna(self):
        links=module.flatten(SOURCE)
        fake_payload=b"Gene\tpat1ta\tpat2ta\nAPOA1\t1\t0.3\n"
        with patch.object(module,"get_file_link",return_value=({}, "https://pdc.cancer.gov/f")):
            with patch.object(module,"bounded_prefix",return_value=(fake_payload,206,True)):
                result=module.header_linkage("PDC000198",links)
        self.assertEqual(result["matrix_header_column_count"],3)
        self.assertFalse(result["rna_protein_join_verified"])
        self.assertFalse(result["matrix_gene_expression_values_verified"])

if __name__=="__main__":
    unittest.main()
