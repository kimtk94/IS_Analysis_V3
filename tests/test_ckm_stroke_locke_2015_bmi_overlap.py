"""Offline original Locke 2015 BMI SNP/N/allele QC, no network."""
import gzip
import io
import hashlib
from unittest.mock import patch
import unittest

from scripts.ckm_stroke.audit_locke_2015_bmi_overlap import (
    HEAD,DigestReader,compare_source,audit, SOURCE_URL,SOURCE_SIZE
)


def raw_rows():
    return "\t".join(HEAD)+"\n"+"\n".join([
      "rs1\tA\tC\t0.6\t0.1\t0.05\t0.00001\t233000",
      "rs2\tC\tA\t0.4\t0.1\t0.05\t0.01\t230000",
      "rs3\tT\tG\t0.5\t-0.1\t0.05\t0.01\t50000",
      "rs4\tG\tT\t0.1\t0.1\t0.05\t0.01\t100000",
      "rs5\tG\tA\t0.2\t0.1\t0.05\t0.01\t120000",
      "rs6\tA\tT\t0.3\t0.1\t0.05\t0.01\t125000",
      "rs7\tA\tC\t0.2\t0.1\t0.05\t0.01\t155000",
    ])+"\n"


class Response(io.BytesIO):
    def __init__(self, raw, length=None, url=SOURCE_URL):
        super().__init__(raw)
        self.status=200
        self.headers={"Content-Length":str(len(raw) if length is None else length)}
        self._url=url
    def geturl(self):return self._url


class Locke2015Tests(unittest.TestCase):
    def test_all_four_direction_classes_and_variable_n(self):
        p={"rs1":("A","C"),"rs2":("A","C"),"rs3":("A","C"),
           "rs4":("A","C"),"rs5":("A","C"),"rs8":("A","C")}
        r=compare_source(iter(raw_rows().splitlines(keepends=True)),p,max_rows=100)
        c=r["counts"]
        self.assertEqual(c["original_rows"],7)
        self.assertEqual(c["candidate_same"],1)
        self.assertEqual(c["candidate_swapped"],1)
        self.assertEqual(c["candidate_complement_same"],1)
        self.assertEqual(c["candidate_complement_swapped"],1)
        self.assertEqual(c["candidate_allele_incompatible"],1)
        self.assertEqual(c["candidate_absent_from_2015_original"],1)
        self.assertEqual(c["candidate_harmonizable"],4)
        self.assertEqual(r["harmonizable_variant_N_summary"]["max"],233000)

    def test_duplicate_rsid_quarantines_all_rows_for_that_id(self):
        original=raw_rows().splitlines(keepends=True)
        result=compare_source(iter(original+[original[1]]),{"rs1":("A","C"),"rs2":("A","C")})
        c=result["counts"]
        self.assertEqual(c["original_duplicate_rsid_rows"],1)
        self.assertEqual(c["original_distinct_ambiguous_duplicate_rsids"],1)
        self.assertEqual(c["candidate_duplicate_rsid_excluded"],1)
        self.assertEqual(c["candidate_harmonizable"],1)
        self.assertEqual(c["candidate_swapped"],1)
        self.assertEqual(c.get("candidate_same",0),0)

    def test_zero_se_rejected(self):
        source="\t".join(HEAD)+"\nrs1\tA\tC\t0.6\t0.1\t0\t0.05\t233000\nrs2\tA\tC\t0.6\t0.1\t0.02\t0.05\t233000\n"
        r=compare_source(iter(source.splitlines(keepends=True)),{"rs1":("A","C"),"rs2":("A","C")})
        self.assertEqual(r["counts"]["candidate_match_with_invalid_original_numeric"],1)
        self.assertEqual(r["counts"]["candidate_harmonizable"],1)
    def test_invalid_numeric_has_no_silent_success(self):
        source="\t".join(HEAD)+"\nrs1\tA\tC\t0.6\t0.1\t0\t0.05\t233000\n"
        with self.assertRaisesRegex(ValueError,"no harmonizable"):
            compare_source(iter(source.splitlines(keepends=True)),{"rs1":("A","C")})

    def test_fractional_sample_N_and_missing_hapmap_frequency_are_preserved(self):
        original="\t".join(HEAD)+"\nrs999\tA\tC\tNA\t0.003\t0.008\t0.9502\t72852.4\n"
        stats=compare_source(iter(original.splitlines(keepends=True)),{"rs999":("A","C")})
        c=stats["counts"]
        self.assertEqual(c["candidate_harmonizable"],1)
        self.assertEqual(c["candidate_fractional_N_rows"],1)
        self.assertEqual(c["candidate_hapmap_CEU_frequency_missing_or_invalid"],1)
        self.assertEqual(stats["harmonizable_variant_N_summary"]["median"],72852.4)

    def test_source_raw_sha_and_integrity(self):
        block=gzip.compress(raw_rows().encode("ascii"))
        rd=DigestReader(Response(block))
        with gzip.GzipFile(fileobj=rd,mode="rb") as f:
            while f.read(512):pass
        while rd.read(1024):pass
        self.assertEqual(rd.n,len(block))
        self.assertEqual(rd.sha.hexdigest(),hashlib.sha256(block).hexdigest())

    def test_invalid_original_header_rejected(self):
        with self.assertRaisesRegex(ValueError,"header"):
            compare_source(iter(["SNP A1 A2\n"]),{"rs1":("A","C")})

    def test_bad_length_is_blocking(self):
        with patch("scripts.ckm_stroke.audit_locke_2015_bmi_overlap.panel_from_candidate",return_value={"rs1":("A","C")}),patch(
            "scripts.ckm_stroke.audit_locke_2015_bmi_overlap.SOURCE_SIZE",35):
            def fetch(req,timeout):
                return Response(gzip.compress(raw_rows().encode()),length=34)
            result=audit(None,None,fetcher=fetch)
        self.assertEqual(result["technical_status"],"HOLD_SOURCE_OR_MATCH_QC")
        self.assertFalse(result["ready_for_ldsc"])
        self.assertIn("compressed size",result["error_detail"])


if __name__=="__main__":
    unittest.main()
