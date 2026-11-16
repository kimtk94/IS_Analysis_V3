import csv
import gzip
import hashlib
import tempfile
import unittest
from pathlib import Path

from scripts.ckm_stroke import glgc_chr_shard as shard

HEADER=("rsID\tCHROM\tPOS_b37\tREF\tALT\tN\tN_studies\tPOOLED_ALT_AF\t"
        "EFFECT_SIZE\tSE\tpvalue_neg_log10\tpvalue\tpvalue_neg_log10_GC\tpvalue_GC")
SNP1="rs1\t22\t20000018\tG\tA\t1000\t2\t0.3\t0.1\t0.1\t0.498\t0.31731050786291415\t0.498\t0.31731050786291415"
SNP2="rs2\t22\t20000030\tA\tT\t1000\t2\t0.5\t-0.2\t0.1\t1.34\t0.04550026389635842\t1.34\t0.04550026389635842"
INDEL="rs3\t22\t20000031\tG\tAC\t1000\t2\t0.2\t0.1\t0.1\t0.498\t0.31731050786291415\t0.498\t0.31731050786291415"
SYMBOLIC="NA\t22\t20000032\tR\tI\t1000\t2\t0.2\t0.1\t0.1\t0.498\t0.31731050786291415\t0.498\t0.31731050786291415"

class FakeTabix:
    def __init__(self,url): self.url=url
    def __enter__(self): return self
    def __exit__(self,*args): return False
    def fetch(self,chrom):
        assert chrom=="22"
        return iter([SNP1,SNP2,INDEL,SYMBOLIC])

class DuplicateTabix(FakeTabix):
    def fetch(self,chrom):
        return iter([SNP1,SNP1,SNP2])

class InvalidTabix(FakeTabix):
    def fetch(self,chrom):
        return iter([SNP1,SNP2.replace("\t0.04550026389635842\t","\tNA\t")])

class MissingNTabix(FakeTabix):
    def fetch(self,chrom):
        return iter([SNP1.replace("\t1000\t2\t","\tNA\t2\t"),SNP2])

class MissingLog10PTabix(FakeTabix):
    def fetch(self,chrom):
        return iter([SNP1.replace("\t0.498\t0.31731050786291415\t","\tNA\t0.31731050786291415\t"),SNP2])

def header(url):
    return {"status":"PASS","header":HEADER,"origin":"official_gzip_prefix"}

class TestGLGCChromShard(unittest.TestCase):
    def test_small_real_format_shard(self):
        with tempfile.TemporaryDirectory() as d:
            r=shard.chromosome_qc("LDL_C","22",factory=FakeTabix,
                                  header_fetcher=header,work_dir=d)
        self.assertEqual(r["status"],"PASS")
        self.assertEqual(r["counts"]["original_rows"],4)
        self.assertEqual(r["counts"]["retained_autosomal_snp"],2)
        self.assertEqual(r["counts"]["non_snp"],1)
        self.assertEqual(r["counts"]["symbolic_or_ambiguous_non_snp"],1)
        self.assertEqual(r["counts"]["palindromic_ambiguous_af_0_42_to_0_58"],1)
        self.assertEqual(r["variant_qc"]["duplicate_excess_rows"],0)
        self.assertFalse(r["full_file_source_sha256_verified"])

    def test_duplicate_locus_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            r=shard.chromosome_qc("TG","22",factory=DuplicateTabix,
                                  header_fetcher=header,work_dir=d)
        self.assertEqual(r["status"],"HOLD")
        self.assertEqual(r["variant_qc"]["duplicate_excess_rows"],1)

    def test_missing_allele_frequency_is_rejected(self):
        bad_header=HEADER.replace("POOLED_ALT_AF\t","")
        r=shard.chromosome_qc("HDL_C","22",factory=FakeTabix,
                              header_fetcher=lambda url: {"status":"PASS","header":bad_header})
        self.assertEqual(r["status"],"FAIL")

    def test_unavailable_original_header_holds(self):
        r=shard.chromosome_qc("HDL_C","22",factory=FakeTabix,
                              header_fetcher=lambda url: {"status":"HOLD","reason":"no header"})
        self.assertEqual(r["status"],"HOLD")

    def test_source_snp_export_retains_sample_n_eaf_rsid_and_logp(self):
        with tempfile.TemporaryDirectory() as d:
            dest=Path(d)/"export"
            r=shard.chromosome_qc("LDL_C","22",factory=FakeTabix,
                                  header_fetcher=header,work_dir=d,export_dir=dest)
            self.assertEqual(r["status"],"PASS")
            self.assertTrue(r["normalized_source_snp_data_persisted"])
            self.assertFalse(r["analysis_ready_data_persisted"])
            self.assertEqual(r["source_snp_export_harmonization"],"NOT_RUN")
            target=Path(r["source_snp_export_path"])
            self.assertTrue(target.is_file())
            with gzip.open(target,"rt",encoding="utf-8") as fh:
                rows=list(csv.DictReader(fh,delimiter="\t"))
            self.assertEqual(len(rows),2)
            self.assertEqual(rows[0]["RSID"],"rs1")
            self.assertEqual(rows[0]["N"],"1000")
            self.assertEqual(float(rows[0]["EAF"]),0.3)
            self.assertEqual(float(rows[0]["NEG_LOG10P"]),0.498)
            self.assertEqual(hashlib.sha256(target.read_bytes()).hexdigest(),
                             r["source_snp_export_sha256"])

    def test_duplicate_snp_can_never_be_exported(self):
        with tempfile.TemporaryDirectory() as d:
            dest=Path(d)/"export"
            r=shard.chromosome_qc("TG","22",factory=DuplicateTabix,
                                  header_fetcher=header,work_dir=d,export_dir=dest)
            self.assertEqual(r["status"],"HOLD")
            self.assertFalse(r["normalized_source_snp_data_persisted"])
            self.assertFalse(dest.exists())

    def test_missing_n_does_not_silently_export_for_downstream(self):
        with tempfile.TemporaryDirectory() as d:
            dest=Path(d)/"export"
            r=shard.chromosome_qc("TG","22",factory=MissingNTabix,
                                  header_fetcher=header,work_dir=d,export_dir=dest)
            self.assertEqual(r["status"],"PASS")
            self.assertEqual(r["source_snp_export_status"],"HOLD")
            self.assertEqual(r["counts"]["auxiliary_missing_or_invalid_n_rows"],1)
            self.assertFalse(dest.exists())

    def test_missing_minus_log10_p_blocks_optional_export(self):
        with tempfile.TemporaryDirectory() as d:
            dest=Path(d)/"export"
            r=shard.chromosome_qc("TG","22",factory=MissingLog10PTabix,
                                  header_fetcher=header,work_dir=d,export_dir=dest)
            self.assertEqual(r["status"],"PASS")
            self.assertEqual(r["counts"]["auxiliary_missing_or_invalid_log10p_rows"],1)
            self.assertEqual(r["source_snp_export_status"],"HOLD")
            self.assertFalse(dest.exists())

    def test_invalid_p_is_recorded(self):
        with tempfile.TemporaryDirectory() as d:
            r=shard.chromosome_qc("HDL_C","22",factory=InvalidTabix,
                                  header_fetcher=header,work_dir=d)
        self.assertEqual(r["status"],"HOLD")
        self.assertEqual(r["counts"]["invalid_effect_or_se"],1)

if __name__=="__main__":
    unittest.main()
