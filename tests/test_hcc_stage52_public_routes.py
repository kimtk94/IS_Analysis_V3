#!/usr/bin/env python3
import pathlib
import sys
import unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage52_public_routes import (SITES,FILES,audit,decision,read_public,
                                         scan_js,public_route_only,static_js_url)

class Response:
    def __init__(self,url,raw,mime,partial=True):
        self.status=206 if partial else 200
        self.url=url
        self.raw=raw
        self.headers={"Content-Type":mime}
        if partial:self.headers["Content-Range"]=f"bytes 0-{len(raw)-1}/{len(raw)}"
    def geturl(self):return self.url
    def read(self,n):return self.raw[:n]
    def __enter__(self):return self
    def __exit__(self,*args):return False

class Tests(unittest.TestCase):
    def test_source_bounded_anonymous_literal_extraction(self):
        raw=(b'const x="/node/api/project/detail"; const y="/node/api/sample/info";'
             b'const z="/node/api/OES012345";')
        x=scan_js(raw)
        self.assertIn("/node/api/project/detail",x["candidate_public_route_literals"])
        self.assertNotIn("/node/api/OES012345",x["candidate_public_route_literals"])
        self.assertFalse(x["file_payload_or_authenticated_api_verified"])
    def test_only_source_literal_paths(self):
        self.assertEqual(public_route_only("/api/project/detail"),"/api/project/detail")
        self.assertIsNone(public_route_only("/api/project?secret=1"))
        self.assertIsNone(public_route_only("/api/OED063528"))
        self.assertIsNone(public_route_only("https://evil.example/api/download"))
        self.assertIsNone(public_route_only("/api/project/../files"))
    def test_origin_port_and_credentials(self):
        self.assertEqual(static_js_url("hccdb25","/static/js/app.123.js"),
                         "http://lifeome.net:809/static/js/app.123.js")
        self.assertIsNone(static_js_url("hccdb25","http://lifeome.net:808/static/js/app.js"))
        self.assertIsNone(static_js_url("hccdb25","//evil.example/a.js"))
        self.assertIsNone(static_js_url("node","https://evil.example/a.js"))
        self.assertIsNone(static_js_url("node","/node/app.js?token=abc"))
    def test_bundle_chunk_and_triplet_references_are_not_files(self):
        js=('var x="chunk-2feaa";'
            'const f="HCCDB25_mRNA_level3.txt";'
            'const s="HCCDB25.sample.txt";'
            'const p="HCCDB25.patient.txt"').encode()
        x=scan_js(js)
        self.assertEqual(x["chunk_module_name_hints_n"],1)
        self.assertTrue(all(x["exact_triplet_filename_references"][f]==1 for f in FILES))
        self.assertFalse(x["file_payload_or_authenticated_api_verified"])
    def test_mock_complete_public_shell_is_still_hold(self):
        root=SITES["hccdb25"][0]
        js="http://lifeome.net:809/static/js/app.js"
        pages={root:(b'<html><script src="/static/js/app.js"></script></html>',"text/html"),
               js:(b'const p="/download"; const r="/repository";',"application/javascript")}
        def opener(req,timeout=None):
            if req.full_url not in pages:raise ValueError("unexpected URL")
            body,mime=pages[req.full_url]
            return Response(req.full_url,body,mime)
        x=audit("hccdb25",opener=opener)
        self.assertEqual(x["status"],"PUBLIC_STATIC_ROUTE_EVIDENCE_COLLECTED")
        self.assertEqual(x["identified_public_routes_n"],2)
        self.assertFalse(x["rna_expression_matrix_downloaded"])
        self.assertEqual(decision(x)["G30_C"],"HOLD_HCCDB25_NODE_RNA_ORIGINAL_FILES")
    def test_modified_go_claim_is_rejected(self):
        with self.assertRaises(ValueError):
            decision({"site":"node","rna_expression_matrix_downloaded":True,
                      "official_pdc_aliquot_rna_sample_join":False})
    def test_malformed_range_denied(self):
        url=SITES["node"][0]
        def opener(req,timeout=None):
            x=Response(url,b"hello","text/html")
            x.headers["Content-Range"]="bytes 10-14/20"
            return x
        with self.assertRaises(ValueError):
            read_public(url,100,opener=opener)
    def test_crosshost_redirect_denied(self):
        url=SITES["node"][0]
        def opener(req,timeout=None):
            return Response("https://evil.example/",b"hello","text/html")
        with self.assertRaises(ValueError):
            read_public(url,100,opener=opener)

if __name__=="__main__":unittest.main()
