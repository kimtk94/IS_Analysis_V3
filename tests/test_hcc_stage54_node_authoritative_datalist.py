#!/usr/bin/env python3
import json, pathlib, sys, unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage54_node_authoritative_datalist import (
    PROJECT,SEARCH,DATA_LIST,OED,SourceError,search_payload,source_request,
    search_candidates,project_general_url,datalist_payload,parse_file_page,
    summarize_records,inspect_node,inspect_niaid,niaid_query_url,niaid_summary,review)

class MockResponse:
    def __init__(self,url,content,status=200):
        self.url=url;self.status=status;self.content=json.dumps(content).encode()
    def geturl(self):return self.url
    def read(self,n):return self.content[:n]
    def __enter__(self):return self
    def __exit__(self,*args):return False

def fake_open(routes):
    def f(request,timeout=None):
        payload=json.loads(request.data) if request.data is not None else None
        k=(request.full_url,json.dumps(payload,sort_keys=True) if payload is not None else "")
        data=routes[k]
        return MockResponse(request.full_url,data)
    return f

def key(url,body=None):
    return url, json.dumps(body,sort_keys=True) if body is not None else ""

class Checks(unittest.TestCase):
    def test_only_curated_source_http_routes(self):
        with self.assertRaises(SourceError):
            source_request("https://attacker.invalid/a",opener=lambda *a,**k:None)
        with self.assertRaises(SourceError):
            source_request("https://www.biosino.org/node/api/private/delete",opener=lambda *a,**k:None)
    def test_browse_json_carries_expected_project_filter(self):
        self.assertEqual(search_payload()["queryWord"],PROJECT)
        self.assertEqual(search_payload()["leftStatQueries"][0]["name"],"project")
    def test_candidate_needs_general_info_crosscheck(self):
        page={"data":{"pageInfo":{"totalPages":1,"content":[
            {"id":"55","projectNo":"OEP000321"},{"id":"56","projectNo":"OEP123123"}]}}}
        res=search_candidates(page)
        self.assertEqual(res["candidates"],["55","56"])
        self.assertEqual(res["explicit_accession_matches_in_search_n"],1)
    def test_source_http_json_full_pass_on_verified_project(self):
        browse={"data":{"pageInfo":{"totalPages":1,"content":[{"id":55,"projectNo":PROJECT}]}}}
        general={"data":{"projectNo":PROJECT}}
        listing={"data":{"total":2,"list":[
            {"datNo":"OED001111","security":"Public","name":"HCCDB25_mRNA_level3.txt"},
            {"datNo":"OED001112","security":"public","name":"HCCDB25.sample.txt"}]}}
        routes={key(SEARCH,search_payload()):browse,
                key(project_general_url("55")):general,
                key(DATA_LIST,datalist_payload("55",1)):listing}
        result=inspect_node(fake_open(routes))
        self.assertEqual(result["status"],"NODE_SOURCE_AUTHORIZED_DATALIST_COMPLETE")
        self.assertTrue(result["source_authorized_public_OED_list"])
        self.assertEqual(result["official_datalist"]["public_OED_ids_verified_n"],2)
        self.assertEqual(result["official_datalist"]["public_exact_triplet_filename_count"]["HCCDB25.sample.txt"],1)
        self.assertTrue(review(result)["G30_D"].startswith("STOP"))
    def test_wrong_general_project_never_issues_datalist(self):
        browse={"data":{"pageInfo":{"totalPages":1,"content":[{"id":55}]}}}
        routes={key(SEARCH,search_payload()):browse,
                key(project_general_url("55")):{"data":{"projectNo":"OEP000322"}}}
        x=inspect_node(fake_open(routes))
        self.assertEqual(x["status"],"NODE_PROJECT_IDENTITY_UNCONFIRMED")
        self.assertFalse(x["source_authorized_project_identity"])
    def test_missing_data_body_is_source_hold(self):
        routes={key(SEARCH,search_payload()):{"code":500,"msg":"transient"}}
        x=inspect_node(fake_open(routes))
        self.assertEqual(x["status"],"NODE_API_SOURCE_HOLD")
        self.assertEqual(x["error_code"],"NODE_RESPONSE_DATA_MISSING")
    def test_nonpublic_file_never_exposed(self):
        x=summarize_records([[{"datNo":"OED004444","security":"restricted","name":"HCCDB25.patient.txt"},
                             {"datNo":"OED004445","security":"public","name":"other"}]],2)
        self.assertEqual(x["public_OED_ids_verified_n"],1)
        self.assertNotIn("OED004444",json.dumps(x))
        self.assertEqual(x["public_exact_triplet_filename_count"]["HCCDB25.patient.txt"],0)
    def test_duplicates_fail_closed(self):
        with self.assertRaises(SourceError):
            summarize_records([[{"datNo":"OED004444","security":"public"},
                                {"datNo":"OED004444","security":"public"}]],2)
    def test_partial_list_remains_partial(self):
        x=summarize_records([[{"datNo":"OED000111","security":"public"}]],3)
        self.assertFalse(x["complete_file_listing"])
    def test_niaid_not_authoritative_by_itself(self):
        record={"_id":"node_oep000321","identifier":"OEP000321","distribution":[
                {"identifier":"OED000888","contentUrl":"https://www.biosino.org/node/download/node/data/public/OED000888"}]}
        x=niaid_summary({"hits":[record]})
        self.assertEqual(x["exact_NODE_dataset_identity_records_n"],1)
        self.assertEqual(x["indexed_public_OED_like_identifiers_n"],1)
        self.assertFalse(x["source_confirmed_project_membership"])
        self.assertIn("q=OEP000321",niaid_query_url())
    def test_niaid_text_mentions_no_exact_match(self):
        x=niaid_summary({"hits":[{"_id":"node_oep000123","name":"OEP000321",
                                 "identifier":"OEP000123"}]})
        self.assertEqual(x["exact_NODE_dataset_identity_records_n"],0)
    def test_biological_promotion_blocked(self):
        with self.assertRaises(SourceError):
            review({"source_files_downloaded":True,
                    "PDC_patient_biopsy_join_verified":False})

if __name__=="__main__":unittest.main()
