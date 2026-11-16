"""Fail-closed deterministic BMI three-reference derived candidate; offline only."""
import gzip
import hashlib
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from scripts.ckm_stroke.audit_bmi_hm3_overlap import allele_pair
from scripts.ckm_stroke.prepare_bmi_ldsc_candidate import (
    OutputWriter, prepare_rows, CANDIDATE_COLS, audit
)


def source_row(snp="rs1",beta="0.1",se="0.05",p="0.01",alleles=("A","C")):
    fields=["1","100",snp,alleles[0],alleles[1],"0.3",beta,se,p,"600000"]
    return fields


class CandidateTests(unittest.TestCase):
    def panels(self):
        mrc={"rs1":allele_pair("A","C"),
             "rs2":allele_pair("A","T"),
             "rs11576885":allele_pair("A","G"),
             "rs3":allele_pair("A","C")}
        return mrc,dict(mrc),set(mrc)

    def test_deterministic_candidate_signed_z_and_zero_p(self):
        m,z,ld=self.panels()
        fixture=[
            (1,source_row("rs1",beta="-0.4",se="0.1",p="0")),
            (2,source_row("rs2",alleles=("A","T"))),
            (3,source_row("rs11576885",alleles=("A","G"))),
            (4,source_row("rs11576885",alleles=("A","G"))),
            (5,source_row("rs3",alleles=("G","T"))),
            (0,["__SOURCE_COMPLETE_SHA256_PASS__","999"])
        ]
        with TemporaryDirectory() as d:
            dest=Path(d)/"bmi.sumstats.gz"
            with OutputWriter(dest) as w:
                result=prepare_rows(iter(fixture),m,z,ld,w,assert_source_profile=False)
            with gzip.open(dest,"rb") as f:
                text=f.read()
        self.assertTrue(text.startswith(CANDIDATE_COLS))
        self.assertIn(b"rs1\t600000\t-4\tA\tC\n",text)
        self.assertNotIn(b"rs2\t",text)
        self.assertNotIn(b"rs3\t",text)
        self.assertEqual(result["counts"]["candidate_rows"],2)
        self.assertEqual(result["counts"]["duplicate_extra_rows"],1)
        self.assertEqual(result["counts"]["candidate_zero_p_signed_z"],1)
        self.assertEqual(result["counts"]["excluded_palindromic"],1)
        self.assertEqual(result["counts"]["excluded_bmi_allele_mismatch"],1)
        self.assertEqual(result["candidate_plaintext_sha256"],hashlib.sha256(text).hexdigest())
        self.assertFalse(result["ready_for_ldsc"])

    def test_discordant_duplicate_aborts(self):
        m,z,ld=self.panels()
        fixture=[
            (1,source_row("rs11576885",beta="0.1",alleles=("A","G"))),
            (2,source_row("rs11576885",beta="0.2",alleles=("A","G"))),
            (0,["__SOURCE_COMPLETE_SHA256_PASS__","999"])
        ]
        with OutputWriter(None) as writer:
            with self.assertRaisesRegex(ValueError,"pinned duplicate differs"):
                prepare_rows(iter(fixture),m,z,ld,writer,assert_source_profile=False)

    def test_partial_file_removed_on_error(self):
        with TemporaryDirectory() as d:
            target=Path(d)/"bad.sumstats.gz"
            with self.assertRaisesRegex(RuntimeError,"intentional"):
                with OutputWriter(target) as writer:
                    writer.write(CANDIDATE_COLS)
                    raise RuntimeError("intentional")
            self.assertFalse(target.exists())
            self.assertFalse(target.with_name(target.name+".partial").exists())

    def test_empty_common_panel_holds(self):
        m,z,ld=self.panels()
        with OutputWriter(None) as writer:
            with self.assertRaisesRegex(ValueError,"strict EUR reference panel"):
                prepare_rows(iter([]),m,z,ld,writer,assert_source_profile=True)

    def test_unverified_output_in_tracked_repo_is_rejected_without_download(self):
        repo=Path(__file__).resolve().parents[1]
        dest=repo/"docs"/"unsafe.sumstats.gz"
        with patch("scripts.ckm_stroke.prepare_bmi_ldsc_candidate.load_references") as downloader:
            report=audit(out_candidate=dest)
        self.assertEqual(report["technical_status"],"HOLD_CANDIDATE_INPUT_OR_TRANSFORM_ERROR")
        self.assertFalse(downloader.called)
        self.assertFalse(report["ready_for_ldsc"])

    def test_alleles_and_sign_not_flipped_for_harmonized_pair(self):
        m,z,ld=self.panels()
        fixture=[(1,source_row("rs1",beta="-1",se="0.25",alleles=("C","A"))),
                 (0,["__SOURCE_COMPLETE_SHA256_PASS__","100"])]
        with TemporaryDirectory() as d:
            dest=Path(d)/"sign.sumstats.gz"
            with OutputWriter(dest) as w:
                report=prepare_rows(iter(fixture),m,z,ld,w,assert_source_profile=False)
            with gzip.open(dest,"rt",encoding="ascii") as f:
                lines=f.read().splitlines()
        self.assertEqual(lines[1],"rs1\t600000\t-4\tC\tA")
        self.assertFalse(report["ready_for_ldsc"])


if __name__=="__main__":
    unittest.main()
