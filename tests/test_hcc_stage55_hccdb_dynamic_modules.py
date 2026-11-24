#!/usr/bin/env python3
import pathlib,sys,unittest,json
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage55_hccdb_dynamic_modules import (
    ROOT,HTMLSource,allowed_url,safe_manifest_entries,public_static_routes,
    get_limited,audit,review
)
class Fake:
    def __init__(self,url,content,mime):
        self.url=url;self.content=content;self.status=206
        self.headers={"Content-Type":mime,
                      "Content-Range":f"bytes 0-{len(content)-1}/{len(content)}"}
    def geturl(self):return self.url
    def read(self,n):return self.content[:n]
    def __enter__(self):return self
    def __exit__(self,*args):return False
def fetch(sources):
    def opener(req,timeout=None):
        data,mime=sources[req.full_url]
        return Fake(req.full_url,data,mime)
    return opener
class TestStage55(unittest.TestCase):
    def test_allow_only_literal_public_static_js(self):
        self.assertEqual(allowed_url("/static/js/app.js",".js"),
                         "http://lifeome.net:809/static/js/app.js")
        self.assertIsNone(allowed_url("http://evil.invalid/a.js",".js"))
        self.assertIsNone(allowed_url("http://lifeome.net:809/a.js?token=abc",".js"))
        self.assertIsNone(allowed_url("/static/%2e%2e/a.js",".js"))
    def test_route_context_from_explicit_script_only(self):
        src=b'var route="/repository", f="chunk-repository123", api="/api/download/archive";'
        summary=public_static_routes(src)
        self.assertEqual(summary["route_context"]["/repository"]["literal_occurrences_n"],1)
        self.assertGreater(summary["route_context"]["/repository"]["first_source_locations"][0]["module_hints_near_route"].__len__(),0)
        self.assertFalse(summary["public_download_file_verified"])
    def test_triplet_file_names_are_hints_only(self):
        summary=public_static_routes(b'HCCDB25_mRNA_level3.txt HCCDB25.sample.txt HCCDB25.patient.txt')
        self.assertEqual(summary["expected_triplet_exact_filename_occurrences"]["HCCDB25.sample.txt"],1)
        self.assertFalse(summary["official_RNA_patient_sample_mapping_verified"])
    def test_manifest_source_only(self):
        files=safe_manifest_entries({"files":{"app":"/static/js/a.js"},"unexpected":["/static/js/b.js"]})
        self.assertEqual(files,["/static/js/a.js"])
    def test_mock_linked_manifest_extra_chunk(self):
        manifest="/manifest.json"
        js="/static/js/app.js"
        extra="/static/js/chunk-repository.js"
        html=b'<html><link rel="manifest" href="/manifest.json"><script src="/static/js/app.js"></script></html>'
        sources={ROOT:(html,"text/html"),
                 "http://lifeome.net:809"+manifest:(json.dumps({"files":{"repo":extra}}).encode(),"application/json"),
                 "http://lifeome.net:809"+js:(b'var route="/download"',"application/javascript"),
                 "http://lifeome.net:809"+extra:(b'var p="/repository"',"application/javascript")}
        v=audit(fetch(sources))
        self.assertEqual(v["status"],"PUBLIC_HCCDB_FRONTEND_SOURCE_INVENTORY_DONE")
        self.assertEqual(v["full_js_sources_n"],2)
        self.assertEqual(v["route_occurrences_total"]["/repository"],1)
        self.assertEqual(v["route_occurrences_total"]["/download"],1)
        self.assertFalse(v["RNA_original_matrix_downloaded"])
    def test_missing_manifest_does_not_forge_chunk(self):
        sources={ROOT:(b'<html><script src="/static/js/a.js"></script></html>',"text/html"),
                 "http://lifeome.net:809/static/js/a.js":(b'var a="chunk-repository";',"application/javascript")}
        v=audit(fetch(sources))
        self.assertEqual(len(v["asset_checks"]),1)
        self.assertEqual(v["linked_public_manifest_checks"],[])
    def test_untrusted_redirect_is_rejected(self):
        with self.assertRaises(Exception):
            get_limited("http://evil.invalid/script.js",100,opener=lambda *args:None)
    def test_forced_go_blocked(self):
        with self.assertRaises(ValueError):
            review({"RNA_original_matrix_downloaded":True,
                    "official_NODE_PDC_same_biopsy_crosswalk":False})
    def test_no_template_url_when_only_chunk_name(self):
        item=public_static_routes(b'const m="chunk-downloadabc"')
        self.assertEqual(item["explicit_static_js_asset_paths_n"],0)
    def test_sample_identifier_in_api_rejected(self):
        item=public_static_routes(b'var z="/api/download/OED001111";')
        self.assertEqual(item["source_literal_api_path_samples"],[])
if __name__=="__main__":unittest.main()
