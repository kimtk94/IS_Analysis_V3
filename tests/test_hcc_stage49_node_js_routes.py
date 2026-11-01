#!/usr/bin/env python3
import pathlib
import sys
import unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage49_node_js_routes import (PROJECT,ScriptParser,asset_url,inspect_js,review,
                                          audit)

class Response:
    def __init__(self,url,raw,mime):
        self.url=url
        self.raw=raw
        self.status=200
        self.headers={"Content-Type":mime}
    def geturl(self):return self.url
    def read(self,n):return self.raw[:n]
    def __enter__(self):return self
    def __exit__(self,*args):return False

class Tests(unittest.TestCase):
    def test_script_extract_and_duplicates(self):
        p=ScriptParser()
        p.feed('<script src="/node/static/a.js"></script><script src="/node/static/a.js"></script>')
        self.assertEqual(len(p.scripts),2)
    def test_static_allowlist(self):
        self.assertEqual(asset_url("/node/static/app.js"),"https://www.biosino.org/node/static/app.js")
        self.assertIsNone(asset_url("https://attacker.invalid/app.js"))
        self.assertIsNone(asset_url("//attacker.invalid/a.js"))
        self.assertIsNone(asset_url("/node/a.js?token=secret"))
        self.assertIsNone(asset_url("/node/%2e%2e/app.js"))
    def test_js_endpoint_hints_do_not_verify_source(self):
        x=inspect_js(b'const a="/node/api/project/detail";OEP000321;')
        self.assertEqual(x["public_api_path_literal_candidates_n"],1)
        self.assertFalse(x["executable_api_metadata_call_verified"])
    def test_script_only_does_not_claim_rna(self):
        html=b'<html><script src="/node/static/app.js"></script></html>'
        js=b'var path="/node/api/project/detail";var name="OEP000321";'
        calls=[]
        def opener(req,timeout=None):
            calls.append(req.full_url)
            if req.full_url==PROJECT:return Response(PROJECT,html,"text/html")
            return Response(req.full_url,js,"application/javascript")
        x=audit(opener=opener)
        self.assertEqual(x["status"],"PUBLIC_JS_ROUTE_CENSUS_COMPLETED")
        self.assertEqual(len(calls),2)
        self.assertFalse(x["rna_original_expression_file_verified"])
        self.assertEqual(review(x)["G30_D"],"STOP_NO_SAME_BIOPSY_RNA_PROTEOME_JOIN")
    def test_untrusted_redirect_rejected(self):
        def opener(req,timeout=None):
            return Response("https://attacker.invalid/",b"<html></html>","text/html")
        x=audit(opener=opener)
        self.assertEqual(x["status"],"SOURCE_ACCESS_HOLD")
    def test_unapproved_go_rejected(self):
        with self.assertRaises(ValueError):
            review({"status":"PUBLIC_JS_ROUTE_CENSUS_COMPLETED",
                    "rna_original_expression_file_verified":True,
                    "same_patient_PDC_RNA_join_verified":False})
    def test_foreign_script_dropped(self):
        html=b'<html><script src="//evil.invalid/foo.js"></script></html>'
        def opener(req,timeout=None):return Response(PROJECT,html,"text/html")
        a=audit(opener=opener)
        self.assertEqual(a["landing"]["disallowed_script_assets_n"],1)
        self.assertEqual(a["assets"],[])

if __name__=="__main__":unittest.main()
