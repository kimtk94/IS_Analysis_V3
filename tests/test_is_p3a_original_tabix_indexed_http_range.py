"""Synthetic verification of a single original Tabix index and BGZF Range source retrieval.

No live EBI network, no original donor IDs, no GWAS/Coloc assumptions.
"""
import gzip
import hashlib
import importlib
import io
import json
import struct
import unittest
import zlib
from unittest.mock import patch

M=importlib.import_module("scripts.is.qtl_original_indexed_http_range")
URL="https://ftp.ebi.ac.uk"+M.QTL_PATH


def make_bgzf(payload):
    obj=zlib.compressobj(level=6,wbits=-15)
    compressed=obj.compress(payload)+obj.flush()
    h=(b"\x1f\x8b\x08\x04"+b"\x00"*4+b"\x00\xff"+
       struct.pack("<H",6)+b"BC"+struct.pack("<H",2)+b"\x00\x00")
    trailer=struct.pack("<II",zlib.crc32(payload),len(payload))
    length=len(h)+len(compressed)+len(trailer)
    assert 20<=length<=65536
    h=h[:-2]+struct.pack("<H",length-1)
    return h+compressed+trailer


def make_tbi(region=None,wrong_header=False,col_seq=2):
    region=region or {"start":101_854_259,"end":105_353_833}
    name=b"10\x00"
    beg=region["start"]-1
    selected_bin=4681+(beg>>14)
    first=(1_000_000<<16)
    last=(1_000_020<<16)|20
    body=b"BAD!" if wrong_header else b"TBI\x01"
    body+=struct.pack("<8i",1,0,col_seq,3,3,ord("#"),1,len(name))+name
    body+=struct.pack("<i",2)
    # The official ROSMAP index contains bin 37450: it is a pseudo-bin
    # with metadata pair values where y<x; NEVER treat as virtual chunks.
    body+=struct.pack("<IiQQQQ",37450,2,0,1000,14_000_000,0)
    body+=struct.pack("<IiQQ",selected_bin,1,first,last)
    n_intv=7000
    body+=struct.pack("<i",n_intv)
    body+=struct.pack("<"+"Q"*n_intv,*([first]*n_intv))
    return gzip.compress(body)


class FakeHTTPResponse:
    def __init__(self,body,status=200,url=URL,content_range=None):
        self.body=body;self.status=status;self.url=url
        self.headers={}
        if content_range:
            self.headers["Content-Range"]=content_range
    def read(self,n=-1):
        return self.body if n<0 else self.body[:n]
    def geturl(self):
        return self.url
    def __enter__(self):return self
    def __exit__(self,*args):return False


class OriginalTabixIndexAndBGZFRangeTests(unittest.TestCase):
    def test_tbi_v1_source_index_and_10_contig(self):
        i=M.decode_original_tbi(make_tbi())
        self.assertIn("10",i["chromosomes"])
        self.assertEqual(i["col_seq"],2)
        self.assertEqual(i["col_beg"],3)

    def test_real_tabix_pseudo_bin_is_metadata_not_virtual_offset(self):
        i=M.decode_original_tbi(make_tbi())
        self.assertNotIn(37450,i["chromosomes"]["10"]["bins"])
        self.assertIn(4681+((101_854_259-1)>>14),
                      i["chromosomes"]["10"]["bins"])

    def test_partial_unneeded_trailing_bgzf_block_is_not_fatal(self):
        a=make_bgzf(b"gene_id\tan\nENSG00000148842\t1120\n")
        b=make_bgzf(b"irrelevant\tdata\n")
        partial=a+b[:-6]
        decoded=M.parse_bgzf_range(
            partial,0,required_last_block_relative_byte_offset=0)
        self.assertEqual(decoded["source_BGZF_blocks_decompressed"],1)

    def test_partial_last_required_bgzf_block_fails_closed(self):
        a=make_bgzf(b"gene_id\tan\nENSG00000148842\t1120\n")
        b=make_bgzf(b"irrelevant\tdata\n")
        with self.assertRaises(ValueError):
            M.parse_bgzf_range(a+b[:-6],0,
                required_last_block_relative_byte_offset=len(a))

    def test_corrupted_index_fails_closed(self):
        with self.assertRaises(ValueError):
            M.decode_original_tbi(make_tbi(wrong_header=True))

    def test_index_schema_wrong_coordinate_column_rejected(self):
        i=M.decode_original_tbi(make_tbi(col_seq=1))
        with self.assertRaises(ValueError):
            M.merged_one_range_from_index(i,{"start":101_854_259,"end":105_353_833})

    def test_one_contiguous_compressed_region_covers_chunk_bounds(self):
        reg={"start":101_854_259,"end":105_353_833}
        result=M.merged_one_range_from_index(M.decode_original_tbi(make_tbi()),reg)
        self.assertEqual(result["compressed_range_start"],1_000_000)
        self.assertEqual(result["compressed_range_end"],1_000_020+65535)
        self.assertTrue(result["one_contiguous_range_spans_all_candidate_tabix_chunks"])

    def test_tbi_bin_hierarchy_includes_selected_leaf_bin(self):
        reg={"start":101_854_259,"end":105_353_833}
        b=M.tbi_bins_for_region(reg["start"],reg["end"])
        self.assertIn(0,b)
        self.assertIn(4681+((reg["start"]-1)>>14),b)

    def test_tbi_query_invalid_coordinate_rejected(self):
        for x,y in [(0,1),(100,99),(1,2**29)]:
            with self.subTest(x=x,y=y),self.assertRaises(ValueError):
                M.tbi_bins_for_region(x,y)

    def test_bgzf_single_block_parses_exact_bytes(self):
        raw="gene_id\tan\nENSG00000148842\t1120\n".encode()
        blob=make_bgzf(raw)
        res=M.parse_bgzf_range(blob,0)
        self.assertEqual(res["data_lines"].read(),raw.decode())
        self.assertEqual(res["source_BGZF_blocks_decompressed"],1)

    def test_bgzf_skips_first_virtual_inside_block(self):
        raw=b"ignored\tline\nENSG00000148842\t1120\n"
        blob=make_bgzf(raw)
        result=M.parse_bgzf_range(blob,len(b"ignored\tline\n"))
        self.assertEqual(result["data_lines"].read(),"ENSG00000148842\t1120\n")

    def test_bgzf_multi_block_no_raw_snp_row_artifact(self):
        a=make_bgzf(b"ENSG00000065613\t1120\n")
        b=make_bgzf(b"ENSG00000148842\t1120\n")
        result=M.parse_bgzf_range(a+b,0)
        self.assertEqual(result["source_BGZF_blocks_decompressed"],2)
        self.assertFalse(result["source_total_QTL_full_file_downloaded"])

    def test_bgzf_malformed_extra_field_rejected(self):
        bad=bytearray(make_bgzf(b"gene\t1120\n"))
        bad[12:14]=b"ZZ"
        with self.assertRaises(ValueError):
            M.parse_bgzf_range(bytes(bad),0)

    def test_bgzf_truncated_block_rejected(self):
        with self.assertRaises(ValueError):
            M.parse_bgzf_range(make_bgzf(b"gene\t1120\n")[:-2],0)

    def test_publisher_index_source_exact_url_and_HTTP_200(self):
        index=make_tbi()
        def fetch(req,timeout):
            self.assertTrue(req.full_url.endswith(".all.tsv.gz.tbi"))
            return FakeHTTPResponse(index)
        self.assertEqual(M.fetch_original_index(URL,fetch),index)
        with self.assertRaises(ValueError):
            M.fetch_original_index(URL,lambda req,timeout: FakeHTTPResponse(
                index,status=403))

    def test_publisher_source_redirect_rejected(self):
        with self.assertRaises(ValueError):
            M.fetch_original_index(URL,lambda req,timeout:FakeHTTPResponse(
                make_tbi(),url="https://example.org/fake.tbi"))

    def test_exact_range_HTTP_206_and_content_range_required(self):
        buf=b"\x1f\x8b"+b"\x00"*48
        limit={"compressed_range_start":100,"compressed_range_end":149}
        def fetch(req,timeout):
            self.assertEqual(req.headers.get("Range"),"bytes=100-149")
            return FakeHTTPResponse(buf,status=206,
                                    content_range="bytes 100-149/100000")
        output,provenance=M.fetch_contiguous_range(URL,limit,fetch)
        self.assertEqual(output,buf)
        self.assertEqual(provenance["publisher_server_content_range"],
                         "bytes 100-149/100000")
        with self.assertRaises(ValueError):
            M.fetch_contiguous_range(URL,limit,lambda req,timeout:
                FakeHTTPResponse(buf,status=200,content_range="bytes 100-149/100000"))

    def test_range_source_cross_host_rejected(self):
        with self.assertRaises(ValueError):
            M.fetch_contiguous_range("https://example.com"+M.QTL_PATH,
                                     {"compressed_range_start":100,
                                      "compressed_range_end":149})

    def test_indexed_source_two_requests_no_retry_without_network(self):
        blob=make_bgzf(b"ENSG00000148842\t1120\n")
        with patch.object(M,"fetch_original_index",return_value=make_tbi()) as index, \
             patch.object(M,"merged_one_range_from_index",return_value={
                 "compressed_range_start":100,
                 "compressed_range_end":100+len(blob)-1,
                 "compressed_range_bytes":len(blob),
                 "first_virtual_offset_within_bgzf_block":0
             }) as indexed, \
             patch.object(M,"fetch_contiguous_range",return_value=(
                 blob,{"publisher_compressed_range_SHA256":hashlib.sha256(blob).hexdigest()}
             )) as source:
            lines,provenance=M.source_indexed_range(
                URL,{"start":101_854_259,"end":105_353_833},interval_seconds=0)
        self.assertEqual(lines.read(),"ENSG00000148842\t1120\n")
        self.assertEqual(provenance["requests_to_original_EBI_source"],2)
        self.assertEqual(provenance["automatic_retry_count"],0)
        self.assertEqual(index.call_count,1)
        self.assertEqual(indexed.call_count,1)
        self.assertEqual(source.call_count,1)


if __name__=="__main__":
    unittest.main()
