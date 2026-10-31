"""Only GraphQL public API field names, never guessed GCT endpoints."""
import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage36b_pdc_public_gct_api_schema import introspect,verify

class Stage36bTests(unittest.TestCase):
    def test_schema_fields_are_public_only(self):
        def fake(q):
            if "__schema" in q:
                return {"__schema":{"queryType":{"name":"Query"}}},"sha1"
            return {"__type":{"fields":[
                {"name":"quantDataMatrix","args":[{"name":"pdc_study_id"},{"name":"data_type"}]},
                {"name":"biospecimenPerStudy","args":[]},
                {"name":"heatmapExport","args":[{"name":"study_id"}]}]}},"sha2"
        r=introspect(fake)
        self.assertEqual(r["status"],"PUBLIC_GRAPHQL_SCHEMA_INTROSPECTED")
        self.assertEqual(r["potential_GCT_or_quant_routes_n"],2)
        self.assertEqual(r["GCT_specific_source_fields_n"],1)
        self.assertFalse(r["source_schema_is_original_study_GCT_file"])
        self.assertNotIn("biospecimenPerStudy",str(r["public_api_relevant_fields"]))

    def test_empty_schema_is_hold(self):
        x=introspect(lambda q:({"__schema":{"queryType":{}}},"xx"))
        self.assertEqual(x["status"],"ACCESS_HOLD")
        self.assertEqual(x["failure_class"],"ValueError")

    def test_graphql_error_class_only(self):
        def fail(_):
            raise RuntimeError("private sample payload which must not be logged")
        x=introspect(fail)
        self.assertEqual(x["failure_class"],"RuntimeError")
        self.assertNotIn("private",str(x))

    def test_biology_never_claimed_from_schema(self):
        d={"schema":"hcc-stage36b-PDC-public-quant-GCT-GraphQL-route-v1",
           "evidence":{"status":"ACCESS_HOLD",
                       "source_schema_is_original_study_GCT_file":False,
                       "GCT_330_biological_mapping_verified":False,
                       "patient_ids_or_data_exported":False}}
        self.assertEqual(verify(d)["G30_D"],"STOP_NO_MATCHED_EXTERNAL_RNA_PROTEIN")
        d["evidence"]["GCT_330_biological_mapping_verified"]=True
        with self.assertRaises(ValueError):verify(d)

if __name__=="__main__":unittest.main()
