"""Stage36 conservative PDC study GCT source and 330-aliquot proof contract."""
import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage36_original_pdc_gct_source import (
    page_inventory,is_safe_page,classify_anchor,parse_gct13_header,
    alias_map,compare_quant_column_names,review
)

class Stage36Tests(unittest.TestCase):
    def gct(self):
        raw=(
            "#1.3\n"
            "2\t2\t1\t5\n"
            "id\tDescription\tSAMPLE_A\tSAMPLE_B\n"
            "aliquot_id\tna\tU1\tU2\n"
            "aliquot_submitter_id\tna\tX1\tX2\n"
            "sample_type\tna\tPrimary Tumor\tSolid Tissue Normal\n"
            "case_id\tna\tC1\tC2\n"
            "aliquot_alias\tna\tQUANT_A\tQUANT_B\n"
            "P1\tprotein\t2.0\t3.0\n"
        )
        return parse_gct13_header(raw.encode())

    def official(self):
        return {
            "aliquot_id":{"U1":[{"aliquot":"U1"}],"U2":[{"aliquot":"U2"}]},
            "aliquot_submitter_id":{"X1":[{"aliquot":"U1"}],"X2":[{"aliquot":"U2"}]},
        }

    def test_gct13_metadata_exact_aliquot_candidates(self):
        a=self.gct()
        self.assertEqual(a["declared_samples_n"],2)
        self.assertEqual(a["column_metadata_count_n"],5)
        x=alias_map(a,self.official())
        self.assertEqual(x["columns_with_one_unambiguous_official_aliquot_n"],2)
        self.assertFalse(x["GCT_to_official_aliquot_330_candidate"])
        self.assertFalse(x["GCT_to_quant_matrix_330_bridge_proven"])

    def test_alias_can_link_two_columns_but_cannot_prove_numeric(self):
        x=compare_quant_column_names(self.gct(),["QUANT_A","QUANT_B"])
        self.assertEqual(x["quant_to_GCT_unique_column_candidates_n"],2)
        self.assertFalse(x["same_original_numeric_assay_certified"])
        self.assertNotIn("QUANT_A",str(x))

    def test_duplicate_official_aliquot_cannot_pass(self):
        g=self.gct()
        g["metadata"]["aliquot_id"]=["U1","U1"]
        g["metadata"]["aliquot_submitter_id"]=["X1","X1"]
        x=alias_map(g,self.official())
        self.assertEqual(x["duplicate_aliquot_columns_n"],1)
        self.assertEqual(x["distinct_official_aliquot_uuids_n"],1)

    def test_disagreeing_uuid_vs_submitter_rejected(self):
        g=self.gct()
        g["metadata"]["aliquot_id"][1]="U1"
        x=alias_map(g,self.official())
        self.assertEqual(x["columns_with_ambiguous_official_aliquot_n"],1)

    def test_missing_column_metadata_never_full_match(self):
        g=self.gct()
        del g["metadata"]["aliquot_submitter_id"]
        x=alias_map(g,self.official())
        self.assertGreater(x["mandatory_source_annotation_fields_missing_n"],0)

    def test_gct_format_rejects_inconsistent_width(self):
        x=b"#1.3\n2\t2\t1\t1\nid\tDescription\tA\naliquot_id\tna\tX1\n"
        with self.assertRaises(ValueError):parse_gct13_header(x)
        with self.assertRaises(ValueError):parse_gct13_header(b"#1.2\n2\t2\nName\tDescription\tA\tB\n")

    def test_public_page_link_not_inferred_from_neighbor_study(self):
        h=b'<a href="/a.gct">PDC999999 heatmap</a><a href="/pdc/study/PDC000198">PDC000198</a>'
        o=page_inventory(h,"https://pdc.cancer.gov/pdc/explore-quantitation-data","PDC000198")
        self.assertEqual(o["target_study_explicit_GCT_links_n"],0)
        self.assertEqual(o["all_page_explicit_GCT_links_unassigned_n"],1)
        self.assertFalse(o["full_original_GCT_retrieved"])

    def test_unapproved_gct_url_never_probed(self):
        self.assertTrue(is_safe_page("https://pdc.cancer.gov/pdc/study/PDC000198"))
        self.assertFalse(is_safe_page("http://pdc.cancer.gov/pdc/study/PDC000198"))
        self.assertFalse(is_safe_page("https://evil.example/pdc/study/PDC000198"))
        self.assertIsNone(classify_anchor("javascript:alert(1)","heatmap","https://pdc.cancer.gov"))

    def test_postflight_does_not_accept_false_biological_go(self):
        d={"schema":"hcc-stage36-original-pdc-GCT-discovery-v1",
           "studies":[{
               "pdc_study_id":p,"status":"SOURCE_GCT_ACCESS_HOLD",
               "GCT_original_file_verified":False,
               "quant_header_to_biological_aliquot_verified":False,
               "any_urls_or_sample_identifiers_exported":False,
               "pages":[{"status":"ACCESS_HOLD","official_page_checked":True}]*3,
               "pages_successful_n":0
           } for p in ("PDC000198","PDC000199")]}
        self.assertEqual(review(d)["G30_D"],"STOP_NO_INDEPENDENT_MATCHED_RNA_PROTEIN")
        d["studies"][0]["GCT_original_file_verified"]=True
        with self.assertRaises(ValueError):review(d)

if __name__=="__main__":unittest.main()
