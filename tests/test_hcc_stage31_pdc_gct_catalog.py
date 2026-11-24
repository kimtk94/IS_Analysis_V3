"""Gao PDC GCT / sample-mapping catalog aggregate evidence gates."""
import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage31_pdc_nonassembly_mapping_catalog import review
class TestCatalog(unittest.TestCase):
    def test_catalog_technical_access_is_not_patient_match(self):
        data={"schema":"hcc-stage31-PDC-other-file-category-v1",
              "studies":[{"pdc_study_id":pid,"status":"FILE_CATEGORY_PREVIEW_COMPLETED",
                          "file_names_exported":False,"sample_ids_exported":False}
                         for pid in ("PDC000198","PDC000199")]}
        r=review(data)
        self.assertEqual(r["G30_D"],"STOP_NO_VALIDATED_EXTERNAL_RNA_PROTEIN_CASES")
    def test_no_identifier_export(self):
        data={"schema":"hcc-stage31-PDC-other-file-category-v1",
              "studies":[{"pdc_study_id":pid,"status":"SOURCE_ACCESS_HOLD",
                          "file_names_exported":False,"sample_ids_exported":False}
                         for pid in ("PDC000198","PDC000199")]}
        data["studies"][0]["sample_ids_exported"]=True
        with self.assertRaises(ValueError):review(data)
if __name__=="__main__":unittest.main()
