"""Evidence-depth report must reject file absence and never promote biology."""
import json
import pathlib
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
from evaluate_parallel_p0_evidence import EVIDENCE_ITEMS, summarize, markdown

def write_fixture(root,name,doc):
    item=next(v for v in EVIDENCE_ITEMS if v[0]==name)
    dest=pathlib.Path(root)/f"hcc-p0-{item[1]}"/item[2]
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(doc),encoding="utf-8")

def complete(root):
    write_fixture(root,"ng_supplement",{"ng":{"attempts":[
        {"status":"XLSX_SOURCE_FETCHED","xlsx_schema":{"sheets":4}}]}})
    write_fixture(root,"ng_quant_workbook_shape",{"ng_supplements":{
        "Supplementary_Data_2":{"status":"PUBLIC_XLSX_SHAPE_VERIFIED"}}})
    write_fixture(root,"ng_all_public_supplement_layouts",{"audit":{
        "supplements_targeted":14,
        "supplements_downloaded_and_layout_screened":14,
        "wide_matrix_geometry_hints":[]}})
    write_fixture(root,"pdc_protein_sample_link",{"records":{"PDC000198":{
        "status":"PUBLIC_MATRIX_METADATA_CORRESPONDENCE_SCREENED",
        "schema_comparison":{"any_exact_normalized_overlap":True}}}})
    write_fixture(root,"node_rna_access",{"NODE":{
        "status":"PUBLIC_LANDING_HTML_ONLY",
        "downloaded_expression_matrix":False}})
    write_fixture(root,"hccdb_rna_access",{"landing_checks":{
        "main":{"status":"LANDING_HTML_ONLY"}},
        "matrix_download_status":"NOT_VERIFIED"})
    write_fixture(root,"geo_metadata",{"records":[{
        "download_status":"FETCHED_PUBLIC_SUPPLEMENT"}]})
    write_fixture(root,"xing_clinical",{"audit":{
        "status":"MATCHED_COMPLETE_CASE_PRESCREEN_AUDITED",
        "n_matched_structurally":152,
        "source_data_sha256_verified":True}})
    write_fixture(root,"xing_etiology",{"audit":{
        "status":"ETIOLOGY_CODEBOOK_PREFLIGHT_AUDITED",
        "matched_protein_to_clinical_rows":152,
        "source_checksums_verified":True,
        "etiology_patient_group_assignment_confirmed":False}})

class Tests(unittest.TestCase):
    def test_all_eight_reads_are_not_proof_of_omics_or_event(self):
        with tempfile.TemporaryDirectory() as tmp:
            complete(tmp)
            r=summarize(tmp)
            self.assertTrue(r["all_source_reports_read"])
            self.assertEqual(r["source_reports_read"],9)
            self.assertEqual(r["evidence"]["xing_protein_clinical"]["structural_matched_pairs"],152)
            self.assertEqual(r["evidence"]["ng_quant_workbook_index"]["workbooks_with_shape_verified"],1)
            self.assertEqual(r["evidence"]["ng_all_public_supplement_layouts"]["supplements_screened"],14)
            self.assertFalse(r["evidence"]["node_rna"]["patient_sample_join_established"])
            self.assertEqual(r["research_claim_gates"]["V02_patient_matched_quantitative_rna_protein"],"BLOCKED")
            self.assertIn("source-evidence depth",markdown(r))
            self.assertFalse(r["patient_identifiers_or_raw_values_in_summary"])

    def test_missing_file_never_becomes_source_verified(self):
        with tempfile.TemporaryDirectory() as tmp:
            complete(tmp)
            item=next(v for v in EVIDENCE_ITEMS if v[0]=="xing_clinical")
            (pathlib.Path(tmp)/f"hcc-p0-{item[1]}"/item[2]).unlink()
            r=summarize(tmp)
            self.assertFalse(r["all_source_reports_read"])
            self.assertEqual(r["source_reports_read"],8)
            self.assertEqual(r["evidence"]["xing_protein_clinical"]["state"],"UNVERIFIED")
            self.assertEqual(r["research_claim_gates"]["V12_event_coding_time_origin_censoring"],"PARTIAL")

    def test_sources_with_unverified_status_are_not_scientific_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            complete(tmp)
            write_fixture(tmp,"ng_supplement",{"ng":{"attempts":[
                {"status":"UNVERIFIED"}]}})
            write_fixture(tmp,"node_rna_access",{"NODE":{"status":"HTTP_ERROR"}})
            r=summarize(tmp)
            self.assertTrue(r["all_source_reports_read"])
            self.assertEqual(r["evidence"]["ng_supplement_structure"]["state"],"UNVERIFIED")
            self.assertEqual(r["evidence"]["node_rna"]["state"],"UNVERIFIED")
            self.assertEqual(r["research_claim_gates"]["V02_patient_matched_quantitative_rna_protein"],"BLOCKED")

    def test_corrupt_json_or_large_file_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            complete(tmp)
            item=next(v for v in EVIDENCE_ITEMS if v[0]=="xing_etiology")
            p=pathlib.Path(tmp)/f"hcc-p0-{item[1]}"/item[2]
            p.write_text("{invalid",encoding="utf-8")
            self.assertEqual(summarize(tmp)["source_reports_read"],8)

if __name__=="__main__":unittest.main()
