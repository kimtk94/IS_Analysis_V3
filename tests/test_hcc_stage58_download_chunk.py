#!/usr/bin/env python3
import hashlib
import pathlib
import sys
import unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
import audit_stage58_download_chunk as m

class Fake:
    def __init__(self,url,data):
        self.url=url
        self.body=data
        self.status=206
        self.headers={"Content-Type":"application/javascript",
                      "Content-Range":"bytes 0-"+str(len(data)-1)+"/"+str(len(data))}
    def geturl(self):return self.url
    def read(self,n):return self.body[:n]
    def __enter__(self):return self
    def __exit__(self,*args):return False

class Tests(unittest.TestCase):
    def test_only_browser_observed_paths_allowed(self):
        self.assertIn("repository.509a9848.js",m.allowed_path(m.CHUNKS["repository"]["path"]))
        with self.assertRaises(ValueError):
            m.allowed_path("/static/js/unverified-chunk.js")
        with self.assertRaises(ValueError):
            m.allowed_path("/api/download/OED012345")
    def test_pinned_observed_source_metadata(self):
        self.assertEqual(m.CHUNKS["repository"]["bytes"],108379)
        self.assertEqual(m.CHUNKS["download"]["bytes"],11335)
        self.assertTrue(all(len(a["sha256"])==64 for a in m.CHUNKS.values()))
    def test_safe_source_path_hints(self):
        self.assertEqual(m.safe_literal("/api/repository/download"),"/api/repository/download")
        self.assertEqual(m.safe_literal("download/bulk"),"/download/bulk")
        self.assertEqual(m.safe_literal("/api/file?id=1"),"/api/file")
        self.assertIsNone(m.safe_literal("/api/download/OED123456"))
        self.assertIsNone(m.safe_literal("https://evil.invalid/download/file"))
        self.assertIsNone(m.safe_literal("/api/private/patient"))
    def test_report_only_aggregate_source_structure(self):
        x=m.inspect_chunk(b'"/api/download/metadata" axios.get("/download/bulk"); HCCDB25.sample.txt',
                          "download")
        self.assertGreaterEqual(x["source_keyword_counts"]["download"],2)
        self.assertFalse(x["verified_RNA_original_matrix_content"])
        self.assertFalse(x["verified_official_NODE_or_PDC_sample_identity"])
    def test_detect_hccdb_triplet_but_never_claim_download(self):
        t=b'HCCDB25_mRNA_level3.txt HCCDB25.sample.txt HCCDB25.patient.txt'
        x=m.inspect_chunk(t,"repository")
        self.assertTrue(all(x["expected_RNA_triplet_file_references"][name]==1
                        for name in m.TARGET_FILES))
        self.assertFalse(x["verified_RNA_original_matrix_content"])
    def test_matching_size_but_wrong_sha_rejected(self):
        data=b"a"*m.CHUNKS["download"]["bytes"]
        def opener(req,timeout=None):return Fake(req.full_url,data)
        with self.assertRaisesRegex(ValueError,"SHA256"):
            m.exact_public_js("download",opener)
    def test_mismatched_size_rejected(self):
        def opener(req,timeout=None):return Fake(req.full_url,b'abc')
        with self.assertRaisesRegex(ValueError,"BYTES_MISMATCH"):
            m.exact_public_js("repository",opener)
    def test_biology_promotion_is_rejected(self):
        with self.assertRaises(ValueError):
            m.review({"rna_matrix_downloaded":True,"official_RNA_PDC_biopsy_identity":False})
        with self.assertRaises(ValueError):
            m.review({"rna_matrix_downloaded":False,"official_RNA_PDC_biopsy_identity":True})
    def test_no_source_file_is_hold(self):
        def opener(req,timeout=None):raise OSError("network")
        x=m.audit("download",opener)
        self.assertEqual(x["status"],"PUBLIC_ROUTE_CHUNK_SOURCE_HOLD")
        self.assertTrue(m.review(x)["G30_C"].startswith("HOLD"))
if __name__=="__main__":unittest.main()
