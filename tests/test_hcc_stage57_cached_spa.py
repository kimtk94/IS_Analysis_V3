#!/usr/bin/env python3
"""Offline Stage57 tests; all sample identifiers below are synthetic."""
import pathlib
import sys
import unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage57_cached_spa import (
    ORIGIN,LANDING,safe_static_path,safe_dynamic_path,declared_script_paths,
    source_cache,get_public,js_evidence,static_audit,review
)

class Fake:
    def __init__(self,url,body,mime="application/javascript",code=206,total=None):
        self.url=url
        self.body=body
        self.status=code
        total=len(body) if total is None else total
        self.headers={"Content-Type":mime,
                      "Content-Range":"bytes 0-"+str(len(body)-1)+"/"+str(total),
                      "Content-Length":str(len(body))}
    def geturl(self):return self.url
    def read(self,n):return self.body[:n]
    def __enter__(self):return self
    def __exit__(self,*args):return False

def mocked(sources):
    def opener(req,timeout=None):
        url=req.full_url
        body,mime=sources[url]
        return Fake(url,body,mime)
    return opener

def sources():
    html=b'<html><script src="/dict.js"></script><script src="/rna.js"></script><script src="/static/js/app.abc.js"></script></html>'
    return {
        LANDING:(html,"text/html"),
        ORIGIN+"/dict.js":(b"var d=1;","application/javascript"),
        ORIGIN+"/rna.js":(b"var x=2;","application/javascript"),
        ORIGIN+"/static/js/app.abc.js":(b'const a="/repository";const b="/download";',"application/javascript"),
    }

class Tests(unittest.TestCase):
    def test_script_source_exact_origin_only(self):
        self.assertEqual(safe_static_path("/dict.js"),"/dict.js")
        self.assertEqual(safe_static_path("/static/js/app.abc.js"),"/static/js/app.abc.js")
        self.assertIsNone(safe_static_path("https://evil.invalid/app.js"))
        self.assertIsNone(safe_static_path("/api/private/a.js"))
        self.assertIsNone(safe_static_path("/static/js/app.js?token=secret"))
        self.assertIsNone(safe_static_path("/static/%2e%2e/private/a.js"))
    def test_dynamic_script_must_be_static_js(self):
        self.assertEqual(safe_dynamic_path("/static/js/123.abcd.js"),"/static/js/123.abcd.js")
        self.assertIsNone(safe_dynamic_path("/dict.js"))
        self.assertIsNone(safe_dynamic_path("/download/HCCDB25.patient.txt"))
        self.assertIsNone(safe_dynamic_path("https://evil.invalid/static/js/app.js"))
    def test_only_html_declarations_no_guessing(self):
        html=(b'<script src="/dict.js"></script><script src="/dict.js"></script>'
              b'<script src="//evil.invalid/attack.js"></script>'
              b'<script src="/static/js/app.abc.js"></script>')
        self.assertEqual(declared_script_paths(html),["/dict.js","/static/js/app.abc.js"])
    def test_bounded_http_content_range_accepts_complete_only(self):
        body,meta=get_public(LANDING,10000,mocked(sources()))
        self.assertEqual(len(body),len(sources()[LANDING][0]))
        self.assertTrue(meta["full_source"])
        def wrong(req,timeout=None):
            b,m=sources()[LANDING]
            return Fake(req.full_url,b,m,total=len(b)+10)
        with self.assertRaises(ValueError):
            get_public(LANDING,10000,wrong)
    def test_source_cached_and_hashed_once(self):
        cache=source_cache(mocked(sources()))
        self.assertEqual(len(cache["paths"]),3)
        self.assertEqual(len(cache["cache"]),3)
        self.assertTrue(all(v["status"]=="SOURCE_COMPLETE" for v in cache["source_checks"]))
        self.assertEqual(len(cache["source_checks"][0]["sha256"]),64)
    def test_missing_public_js_is_not_fabricated(self):
        sample=sources()
        del sample[ORIGIN+"/rna.js"]
        c=source_cache(mocked(sample))
        self.assertEqual(len(c["cache"]),2)
        self.assertEqual(sum(x["status"]=="SOURCE_HOLD" for x in c["source_checks"]),1)
    def test_static_source_inventory_stays_source_hold_for_rna(self):
        s=static_audit(mocked(sources()))
        self.assertTrue(s["all_declared_source_js_complete"])
        self.assertEqual(s["complete_source_js_n"],3)
        self.assertFalse(s["review"]["original_HCCDB25_triplet_verified"])
        self.assertTrue(s["review"]["G30_D"].startswith("STOP"))
    def test_rna_filename_reference_not_contents(self):
        a=js_evidence(b'HCCDB25_mRNA_level3.txt HCCDB25.sample.txt HCCDB25.patient.txt',
                      "/static/js/chunk.example.js")
        self.assertEqual(a["triplet_file_name_mentions"]["HCCDB25.patient.txt"],1)
        self.assertFalse(a["downloaded_RNA_file"])
        self.assertFalse(a["patient_identity_join_verified"])
    def test_private_accession_path_must_be_redacted(self):
        a=js_evidence(b'"/api/download/OED123456" "/api/download/metadata"',
                      "/static/js/normal.js")
        self.assertEqual(a["safe_literal_route_paths"],["/api/download/metadata"])
    def test_false_biology_go_is_never_allowed(self):
        with self.assertRaises(ValueError):
            review({"rna_original_downloaded":True,
                    "pdc_official_patient_join_verified":False})
        with self.assertRaises(ValueError):
            review({"rna_original_downloaded":False,
                    "pdc_official_patient_join_verified":True})

if __name__=="__main__":
    unittest.main()
