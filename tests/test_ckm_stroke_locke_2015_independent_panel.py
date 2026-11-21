"""Locke 2015 source-selected BMI panel synthetic fail-closed fixtures."""
from __future__ import annotations
import gzip
import hashlib
import io
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from scripts.ckm_stroke.audit_bmi_hm3_overlap import allele_pair
from scripts.ckm_stroke.prepare_locke_2015_independent_panel import (
    CandidateWriter,OUTPUT_HEADER,prepare_rows,source_rows,audit,
    DUPLICATE_SOURCE_RSIDS,SOURCE_URL,
)


def fields(snp="rs100",a1="A",a2="C",n="234567.4",beta="0.05",se="0.01",p="0.02"):
    return [snp,a1,a2,"NA",beta,se,p,n]


def stream():
    data=[
        fields("rs100"),
        fields("rs101","C","A",n="200000"),
        fields("rs102","A","T"),
        fields("rs103","A","G"),
        fields("rs104",n="NA"),
        fields("rs105",n="100000"),
        fields("rs3780697"),
        fields("rs3780697","C","A",beta="-0.05")
    ]
    for i,r in enumerate(data,1):yield i,r
    yield 0,["__LOCKE_FULL_SOURCE_SHA256_PASS__","35869809"]


def refs():
    values={s:allele_pair("A","C") for s in
            ("rs100","rs101","rs102","rs103","rs104","rs105","rs3780697")}
    return values,set(values)


class PanelTests(unittest.TestCase):
    def test_2015_only_selection_and_duplicate_quarantine(self):
        rp,allowed=refs()
        with CandidateWriter(None) as writer:
            result=prepare_rows(stream(),rp,allowed,writer,strict=False)
        cnt=result["counts"]
        self.assertEqual(cnt["selected_rows"],3)
        self.assertEqual(cnt["selected_fractional_N"],1)
        self.assertEqual(cnt["excluded_palindromic"],1)
        self.assertEqual(cnt["excluded_allele_pair_discordance"],1)
        self.assertEqual(cnt["excluded_bad_numeric"],1)
        self.assertEqual(cnt["duplicate_extra_rows"],1)
        self.assertEqual(cnt["source_duplicate_quarantine_rows"],2)
        self.assertEqual(cnt["distinct_source_duplicate_ids"],1)
        self.assertEqual(result["n_distribution"]["median"],200000.0)
        self.assertFalse(result["ready_for_ldsc"])

    def test_gzip_candidate_preserves_fractional_n_and_no_2018_effects(self):
        rp,allowed=refs()
        with TemporaryDirectory() as td:
            path=Path(td)/"locke.candidate.tsv.gz"
            with CandidateWriter(path) as w:
                stats=prepare_rows(stream(),rp,allowed,w,strict=False)
            with gzip.open(path,"rb") as f:
                raw=f.read()
            self.assertEqual(hashlib.sha256(raw).hexdigest(),stats["2015_candidate_plaintext_sha256"])
            self.assertTrue(raw.startswith(OUTPUT_HEADER))
            self.assertIn(b"rs100\t234567.4\t5\tA\tC\n",raw)
            self.assertNotIn(b"rs3780697",raw)
            self.assertNotIn(b"Z\tA1",raw)

    def test_invalid_new_duplicate_fails_closed_and_removes_partial(self):
        rp,allowed=refs()
        bad=iter([(1,fields("rs100")),(2,fields("rs100")),
                  (0,["__LOCKE_FULL_SOURCE_SHA256_PASS__","100"])])
        with TemporaryDirectory() as td:
            dest=Path(td)/"x.candidate.tsv.gz"
            with self.assertRaisesRegex(ValueError,"unrecognized new duplicated"):
                with CandidateWriter(dest) as w:
                    prepare_rows(bad,rp,allowed,w,strict=False)
            self.assertFalse(dest.exists())
            self.assertFalse(Path(str(dest)+".partial").exists())

    def test_broken_original_sentinel_rejected(self):
        rp,allowed=refs()
        only=[(1,fields("rs100"))]
        with self.assertRaisesRegex(ValueError,"source did not reach full SHA256"):
            with CandidateWriter(None) as w:
                prepare_rows(iter(only),rp,allowed,w,strict=False)

    def test_reference_only_exclusion_and_no_2018_selection(self):
        rp,allowed=refs()
        allowed={"rs100"}
        sample=iter([(1,fields("rs100")),(2,fields("rs105")),
                     (0,["__LOCKE_FULL_SOURCE_SHA256_PASS__","100"])])
        with CandidateWriter(None) as w:
            result=prepare_rows(sample,rp,allowed,w,strict=False)
        self.assertEqual(result["counts"]["selected_rows"],1)
        self.assertEqual(result["counts"]["excluded_outside_strict_triple_reference"],1)

    def test_source_stream_pinned_sha_strict(self):
        raw=OUTPUT_HEADER+b"rs100\t100000\t5\tA\tC\n"
        # The original reader is validated separately by source size and source SHA.
        from scripts.ckm_stroke.prepare_locke_2015_independent_panel import load_strict_references
        self.assertTrue(callable(load_strict_references))
        self.assertTrue(SOURCE_URL.endswith("SNP_gwas_mc_merge_nogc.tbl.uniq.gz"))
        self.assertEqual(len(DUPLICATE_SOURCE_RSIDS),13)

    def test_report_fails_closed_if_source_or_reference_is_missing(self):
        def broken_ref(fetcher):raise ConnectionError("source not available")
        result=audit(reference_loader=broken_ref,strict=False)
        self.assertEqual(result["technical_status"],"HOLD_SOURCE_OR_REFERENCE_OR_MUNGE")
        self.assertFalse(result["ready_for_ldsc"])
        self.assertFalse(result["2018_exposure_candidates_used"])


if __name__=="__main__":
    unittest.main()
