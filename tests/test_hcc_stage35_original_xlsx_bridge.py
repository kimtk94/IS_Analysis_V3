"""Source-derived XLSX exact overlapping tags are only preflight, not identity."""
import io
import pathlib
import sys
import unittest
from openpyxl import Workbook
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage35_original_experimental_xlsx_bridge import (
    text_tokens,compare_worksheet,in_scope_tokens,inspect_xlsx,review
)

class Stage35Tests(unittest.TestCase):
    @staticmethod
    def source():
        return {
            "official_aliquot_uuid":{"UUID_A","UUID_B"},
            "official_aliquot_submitter":{"AL_1","AL_2"},
            "official_plex":{"plex-1"},
            "official_run":{"run-1"},
            "cdap_analytical_sample":{"cdap-1"},
            "cdap_biological_label":{"CDAP_A","CDAP_B"},
        }
    def test_source_exact_row_candidates_not_certified(self):
        sheet=[("CDAP label","Aliquot ID"),("CDAP_A","UUID_A"),("CDAP_B","UUID_B")]
        r=compare_worksheet(sheet,self.source())
        self.assertEqual(r["single_source_row_pair_candidates_n"],2)
        self.assertEqual(r["distinct_source_labels_with_single_candidate_n"],2)
        self.assertFalse(r["row_cooccurrence_is_mapping_proof"])

    def test_ambiguous_row_never_unique(self):
        sheet=[("CDAP_A","CDAP_B","UUID_A","UUID_B")]
        r=compare_worksheet(sheet,self.source())
        self.assertEqual(r["rows_with_ambiguous_multiple_link_tokens_n"],1)
        self.assertEqual(r["distinct_source_labels_with_single_candidate_n"],0)

    def test_numeric_codes_not_stringified(self):
        self.assertEqual(text_tokens([127, "0127", 127.0]),{"0127"})

    def test_workbook_in_memory_no_values_exported(self):
        b=Workbook()
        s=b.active
        s.append(["CDAP label","Aliquot ID"])
        s.append(["CDAP_A","UUID_A"])
        stream=io.BytesIO()
        b.save(stream)
        result=inspect_xlsx(stream.getvalue(),self.source())
        self.assertEqual(result["worksheet_count_n"],1)
        self.assertEqual(result["worksheet_anonymous_statistics"][0]["distinct_source_labels_with_single_candidate_n"],1)
        self.assertNotIn("CDAP_A",str(result))
        self.assertNotIn("UUID_A",str(result))

    def test_ratio_source_column_selection(self):
        rows=[{"AnalyticalSample":"cdap-1","126":"POOL",
               "127N":"CDAP_A","Ratios":"127N/126"}]
        a,b,n=in_scope_tokens(rows,["126","127N"])
        self.assertEqual((len(a),len(b),n),(1,1,1))
        self.assertEqual(b,{"CDAP_A"})

    def test_postflight_requires_holds(self):
        base={"schema":"hcc-stage35-original-PDC-XLSX-bridge-v1",
             "sources":[{"pdc_study_id":p,"status":"ACCESS_HOLD",
                "file_names_urls_or_original_values_exported":False,
                "source_CDP_to_PDC_mapping_certified":False,
                "quant_header_to_aliquot_certified":False,
                "RNA_protein_patient_join_certified":False} for p in ("PDC000198","PDC000199")]}
        self.assertEqual(review(base)["G30_D"],"STOP_NO_MATCHED_EXTERNAL_PATIENT_MULTIMODAL_DATA")
        base["sources"][0]["quant_header_to_aliquot_certified"]=True
        with self.assertRaises(ValueError):
            review(base)

if __name__=="__main__":
    unittest.main()
