#!/usr/bin/env python3
"""Offline positive / negative HTTP range and conservative RNA evidence tests."""
import pathlib,sys,re,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage50_node_hccdb25_range import (
    SITES,CHUNK,same_origin_asset,read_range,read_asset,summarize_js,inspect,review
)

class Response:
    def __init__(self,url,body,status,headers):
        self.url=url;self.body=body;self.status=status;self.headers=headers
    def geturl(self):return self.url
    def read(self,n):return self.body[:n]
    def __enter__(self):return self
    def __exit__(self,*exc):return False

def public_opener(files,ignore_range=False,bad_range=False,redirect=False):
    def open_request(request,timeout=None):
        url=request.full_url
        source=files[url]
        media="text/html" if b"<html" in source[:100] else "application/javascript"
        if ignore_range:
            return Response(url,source,200,{"Content-Type":media,"Content-Length":str(len(source))})
        m=re.fullmatch(r"bytes=(\d+)-(\d+)",request.get_header("Range"))
        start,end=map(int,m.groups())
        part=source[start:end+1]
        if not part:raise ValueError("synthetic range outside file")
        header=f"bytes {start}-{start+len(part)-1}/{len(source)}"
        if bad_range:header=f"bytes 12-{start+len(part)-1}/{len(source)}"
        return Response("https://evil.invalid/" if redirect else url,part,206,
                        {"Content-Type":media,"Content-Range":header})
    return open_request

class Stage50Tests(unittest.TestCase):
    def test_strict_same_origin(self):
        site=SITES["node"]
        self.assertIsNotNone(same_origin_asset("/node/js/main.js",site))
        self.assertIsNone(same_origin_asset("https://evil.invalid/code.js",site))
        self.assertIsNone(same_origin_asset("//evil.invalid/code.js",site))
        self.assertIsNone(same_origin_asset("/node/%2e%2e/code.js",site))
        self.assertIsNone(same_origin_asset("/node/code.js?api_token=abc",site))
    def test_hccdb_port_809_source(self):
        site=SITES["hccdb25"]
        self.assertEqual(same_origin_asset("/static/js/app.js",site),
                         "http://lifeome.net:809/static/js/app.js")
        self.assertIsNone(same_origin_asset("http://lifeome.net/static/js/app.js",site))
    def test_multi_range_full_file(self):
        url="https://www.biosino.org/node/main.js"
        raw=b"123abc"* (CHUNK//6+20)
        data,meta=read_asset(url,opener=public_opener({url:raw}))
        self.assertEqual(data,raw)
        self.assertTrue(meta["full_source_read_verified"])
        self.assertFalse(meta["bounded_source_prefix_only"])
        self.assertEqual(meta["claimed_total_bytes"],len(raw))
    def test_prefix_cap_never_reports_complete(self):
        url="https://www.biosino.org/node/main.js"
        raw=b"r"* (CHUNK+100)
        data,meta=read_asset(url,limit=CHUNK,opener=public_opener({url:raw}))
        self.assertEqual(len(data),CHUNK)
        self.assertFalse(meta["full_source_read_verified"])
        self.assertTrue(meta["bounded_source_prefix_only"])
    def test_http_200_range_ignored_is_partial(self):
        url="https://www.biosino.org/node/main.js"
        raw=b"a"*(CHUNK+10)
        data,meta=read_asset(url,opener=public_opener({url:raw},ignore_range=True))
        self.assertLessEqual(len(data),CHUNK)
        self.assertFalse(meta["full_source_read_verified"])
    def test_invalid_range_denied(self):
        url="https://www.biosino.org/node/main.js"
        raw=b"ab"*128
        with self.assertRaises(ValueError):
            read_range(url,0,CHUNK+1,opener=public_opener({url:raw}))
        with self.assertRaises(ValueError):
            read_asset(url,opener=public_opener({url:raw},bad_range=True))
    def test_cross_host_redirect_denied(self):
        url="https://www.biosino.org/node/main.js"
        with self.assertRaises(ValueError):
            read_asset(url,opener=public_opener({url:b"abc"},redirect=True))
    def test_html_js_aggregate_stays_hold(self):
        site=SITES["node"]
        root=site["url"];asset="https://www.biosino.org/node/main.js"
        files={root:b'<html><script src="/node/main.js"></script></html>',
               asset:b'var path="/node/api/project/detail"; const n="HCCDB25";'}
        source=inspect("node",opener=public_opener(files))
        self.assertEqual(source["status"],"PUBLIC_STATIC_SOURCE_AUDITED")
        self.assertEqual(source["explicit_literal_api_path_candidates_n"],1)
        self.assertFalse(source["RNA_expression_matrix_retrieved"])
        self.assertTrue(review(source)["G30_C"].startswith("HOLD"))
    def test_source_strings_not_matrix_access(self):
        metrics=summarize_js(b'const z="HCCDB25 RSEM expression";var url="/node/api/sample/metadata"')
        self.assertEqual(metrics["source_HCCDB25_literal_mentions_n"],1)
        self.assertFalse(metrics["expression_matrix_content_verified"])
        self.assertFalse(metrics["public_endpoint_runtime_verified"])
    def test_unverified_join_cannot_be_promoted(self):
        with self.assertRaises(ValueError):
            review({"site":"node","RNA_expression_matrix_retrieved":True,
                    "NODE_to_PDC_case_sample_aliquot_join_verified":False})

if __name__=="__main__":unittest.main()
