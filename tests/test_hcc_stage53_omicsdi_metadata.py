#!/usr/bin/env python3
import json
import pathlib
import sys
import unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage53_omicsdi_metadata import (
    documented_search_url,exact_dataset_url,identifier_match,search_hit_summary,
    file_metadata_summary,inspect_query,review,bounded_json)

class MockResponse:
    status=200
    def __init__(self,url,data,mime="application/json"):
        self.url=url
        self.data=json.dumps(data).encode("utf-8")
        self.headers={"Content-Type":mime}
    def geturl(self):return self.url
    def read(self,n):return self.data[:n]
    def __enter__(self):return self
    def __exit__(self,*args):return False

def mock_opener(routes):
    def opener(req,timeout=None):
        value=routes[req.full_url]
        return MockResponse(req.full_url,value)
    return opener

class Stage53Tests(unittest.TestCase):
    def test_only_documented_queries(self):
        self.assertIn("query=OEP000321",documented_search_url("OEP000321"))
        with self.assertRaises(ValueError):documented_search_url("OED000001")
    def test_only_safe_detail_source_ids(self):
        self.assertEqual(exact_dataset_url("node","OEP000321"),
                         "https://www.omicsdi.org/ws/dataset/node/OEP000321")
        with self.assertRaises(ValueError):exact_dataset_url("https://evil.example","OEP000321")
        with self.assertRaises(ValueError):exact_dataset_url("node","../../secrets")
    def test_text_mentions_not_authorized_detail(self):
        data={"count":1,"datasets":[{"id":"PRJ999999","source":"node",
                                      "description":"may contain OEP000321 transcriptome"}]}
        x=search_hit_summary(data,"OEP000321")
        self.assertEqual(x["exact_primary_accession_hit_n"],0)
        self.assertEqual(x["text_mentions_without_exact_accession_n"],1)
        self.assertEqual(x["exact_source_identifiers"],[])
    def test_exact_match_requires_identifier_field(self):
        good={"id":"OEP000321","source":"node"}
        self.assertTrue(identifier_match(good,"OEP000321"))
        self.assertFalse(identifier_match({"title":"OEP000321"},"OEP000321"))
    def test_proven_index_match_still_hold_without_rna(self):
        q="OEP000321"
        urls={
            documented_search_url(q):{"count":1,"datasets":[{"id":q,"source":"node"}]},
            exact_dataset_url("node",q):{"id":q,"title":"HCC transcriptome"},
            exact_dataset_url("node",q,True):{"files":[{"accession":"OED063528","name":"HCCDB25_mRNA_level3.txt"}]},
        }
        out=inspect_query(q,opener=mock_opener(urls))
        self.assertEqual(out["status"],"INDEX_EXACT_ACCESSION_MATCHED")
        detail=out["source_detail_checks"][0]
        self.assertEqual(detail["status"],"OMDI_FILE_LIST_METADATA_READ")
        self.assertEqual(detail["file_metadata"]["unique_OED_accession_fields_n"],1)
        self.assertEqual(detail["file_metadata"]["HCCDB25_exact_file_names"]["HCCDB25_mRNA_level3.txt"],1)
        self.assertFalse(detail["file_metadata"]["NODE_project_membership_source_certified"])
        self.assertFalse(review(out)["official_NODE_project_OED_data_membership_verified"])
    def test_empty_index_matches_are_not_zero_proof_of_dataset(self):
        q="HCCDB25"
        out=inspect_query(q,opener=mock_opener({
            documented_search_url(q):{"count":0,"datasets":[]}}))
        self.assertEqual(out["status"],"INDEX_QUERY_COMPLETE_NO_EXACT_MATCH")
        self.assertFalse(out["NODE_official_OED_file_id_verified"])
    def test_untrusted_or_truncated_payload_rejected(self):
        with self.assertRaises(ValueError):bounded_json("http://evil.invalid/data",mock_opener({}))
    def test_sensitive_identifiers_not_emitted(self):
        m=file_metadata_summary({"files":[
            {"id":"OED001234","name":"HCCDB25.sample.txt","url":"https://example.org/SAMPLE-P01"}]})
        self.assertEqual(m["unique_OED_accession_fields_n"],1)
        self.assertNotIn("OED001234",json.dumps(m))
        self.assertNotIn("SAMPLE-P01",json.dumps(m))
    def test_duplicate_ids_count_separately(self):
        m=file_metadata_summary([{"id":"OED001234"},{"accession":"OED001234"}])
        self.assertEqual(m["OED_accession_fields_n"],2)
        self.assertEqual(m["unique_OED_accession_fields_n"],1)
    def test_no_clinical_go(self):
        with self.assertRaises(ValueError):
            review({"queried_accession":"OEP000321","NODE_official_OED_file_id_verified":True,
                    "actual_HCCDB25_processed_RNA_file_downloaded":False})

if __name__=="__main__":unittest.main()
