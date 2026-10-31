"""Stage30 official PDC direct quant header↔aliquot test: source joins, ID privacy."""
import io,pathlib,sys,unittest,json
from unittest.mock import patch
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from probe_stage30_gao_quant_api_header_to_biospecimen import (
    official_biospecimen,score,get_quant_header,verify)

def source_rows():
    return [
      {"case_id":"CASE_PRIVATE_1","case_submitter_id":"CASE_SUB_PRIVATE_1",
       "sample_id":"SAMPLE_T_PRIVATE","sample_submitter_id":"S_T_PRIVATE",
       "aliquot_id":"ALIQUOT_T_PRIVATE","aliquot_submitter_id":"AT_PRIVATE",
       "sample_type":"Primary Tumor"},
      {"case_id":"CASE_PRIVATE_1","case_submitter_id":"CASE_SUB_PRIVATE_1",
       "sample_id":"SAMPLE_N_PRIVATE","sample_submitter_id":"S_N_PRIVATE",
       "aliquot_id":"ALIQUOT_N_PRIVATE","aliquot_submitter_id":"AN_PRIVATE",
       "sample_type":"Normal"},
      {"case_id":"CASE_PRIVATE_2","case_submitter_id":"CASE_SUB_PRIVATE_2",
       "sample_id":"SAMPLE_B_PRIVATE","sample_submitter_id":"S_B_PRIVATE",
       "aliquot_id":"ALIQUOT_B_PRIVATE","aliquot_submitter_id":"AB_PRIVATE",
       "sample_type":"Other"},
    ]

class QuantAliasTests(unittest.TestCase):
    def getindex(self):
        with patch("probe_stage30_gao_quant_api_header_to_biospecimen.pdc_query",
                return_value=({"biospecimenPerStudy":source_rows()},"a"*64)):
            # fake small result allowed by patching 100-row input requirement
            return official_biospecimen("PDC000198")
    def test_multi_sample_case_index_is_not_ambiguity(self):
        rows=source_rows()*35
        with patch("probe_stage30_gao_quant_api_header_to_biospecimen.pdc_query",
                return_value=({"biospecimenPerStudy":rows},"a"*64)):
            a=official_biospecimen("PDC000198")
        self.assertEqual(len(a["case_id"]["CASE_PRIVATE_1"]),2)
        x=score(["CASE_PRIVATE_1"],a["case_id"])
        self.assertEqual(x["exact_source_header_alias_n"],1)
        self.assertEqual(x["distinct_exact_matched_samples_n"],2)
        self.assertEqual(x["alias_keys_with_multiple_biospecimen_records_n"],1)
    def test_exact_aliquot_names_bridge(self):
        idx={"AT_PRIVATE":[{"case":"C1","sample":"S1","aliquot":"A1","type":"Tumor"}],
             "AN_PRIVATE":[{"case":"C1","sample":"S2","aliquot":"A2","type":"Normal"}]}
        x=score(["AT_PRIVATE","AN_PRIVATE"],idx)
        self.assertEqual(x["exact_source_header_alias_n"],2)
        self.assertEqual(x["distinct_exact_matched_cases_n"],1)
        self.assertEqual(x["exact_matched_Tumor_alias_n"],1)
        self.assertEqual(x["exact_matched_Normal_alias_n"],1)
        self.assertNotIn("PRIVATE",str(x))
    def test_no_fuzzy_join(self):
        idx={"AA_PRIVATE":[{"case":"C1","sample":"S1","aliquot":"A1","type":"Tumor"}]}
        self.assertEqual(score(["AA_PRIVATE-extra"],idx)["exact_source_header_alias_n"],0)
    def test_header_first_rows_are_private(self):
        raw=b'{"data":{"quantDataMatrix":[["Gene/Aliquot","AT_PRIVATE","AN_PRIVATE"],["GENE_PRIVATE",0.2,-1.1],'
        class Fake:
            status=200
            def __enter__(self):return self
            def __exit__(self,*args):return False
            def read(self,n):return raw[:n]
        with patch("probe_stage30_gao_quant_api_header_to_biospecimen.urllib.request.urlopen",
                   return_value=Fake()):
            summary,ids=get_quant_header("PDC000198")
        self.assertEqual(len(ids),2)
        self.assertEqual(summary["bounded_postheader_numeric_cells_n"],2)
        self.assertTrue(summary["first_header_gene_label_present"])
        self.assertNotIn("PRIVATE",str(summary))
    def test_summary_no_biology_upgrade(self):
        d={"schema":"hcc-stage30-direct-pdc-header-biospecimen-match-v1",
           "studies":[{"pdc_study_id":k,"status":"ACCESS_HOLD",
                "RNA_sample_identity_verified":False,
                "full_quantitative_matrix_verified":False,
                "individual_case_sample_aliquot_or_gene_IDs_exported":False}
                for k in ("PDC000198","PDC000199")]}
        self.assertEqual(verify(d)["G30_D"],"STOP_NO_INDEPENDENT_SAME_PATIENT_RNA_PROTEIN_JOIN")
        d["studies"][0]["RNA_sample_identity_verified"]=True
        with self.assertRaises(ValueError):verify(d)
if __name__=="__main__":unittest.main()
