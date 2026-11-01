"""CDAP *.sample.txt design companion / patient aliquot exact ID bridge guards."""
import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from probe_stage30_gao_sample_design_companion import compare_rows,verify

class SampleTests(unittest.TestCase):
    def test_row_with_both_quant_and_aliquot_exact_names(self):
        table=[["plex","column","aliquot"],
               ["PRIVATE_PLEX","PRIVATE_QUANT_X","PRIVATE_ALIQUOT_T"],
               ["PRIVATE_PLEX","PRIVATE_QUANT_Y","PRIVATE_ALIQUOT_N"]]
        refs={"aliquot_id":{"PRIVATE_ALIQUOT_T":1,"PRIVATE_ALIQUOT_N":2},
              "aliquot_submitter_id":{}}
        x=compare_rows(table,["PRIVATE_QUANT_X","PRIVATE_QUANT_Y"],refs)
        self.assertEqual(x["companion_quant_header_exact_common_tokens_n"],2)
        self.assertEqual(x["companion_official_aliquot_UUID_or_submitter_exact_common_tokens_n"],2)
        self.assertEqual(x["companion_rows_containing_both_quant_alias_and_aliquot_alias_n"],2)
        self.assertNotIn("PRIVATE_",str(x))
        self.assertFalse(x["mapping_complete_and_unambiguous_verified"])
    def test_truncated_access_not_biological_failure(self):
        x={"schema":"hcc-stage30-PDC-sample-design-companion-v1",
           "studies":[{"pdc_study_id":p,
                "status":"NO_SAMPLE_DESIGN_COMPANION_IN_ASSEMBLY_CATEGORY",
                "source_file_urls_rows_and_aliquot_names_exported":False,
                "validated_full_assay_col_to_patient_join":False,
                "RNA_patient_join_verified":False}
                for p in ("PDC000198","PDC000199")]}
        self.assertEqual(verify(x)["G30_B"],"HOLD_NO_VALIDATED_ASSAY_COLUMN_TO_ALIQUOT_JOIN")
        x["studies"][0]["validated_full_assay_col_to_patient_join"]=True
        with self.assertRaises(ValueError):verify(x)
if __name__=="__main__":unittest.main()
