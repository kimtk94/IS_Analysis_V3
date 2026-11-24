#!/usr/bin/env python3
"""Synthetic contract tests: public source artifact results are never fixtures."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from aggregate_gwas_1000g_all22 import (
    EXPECTED_GWAS_ROWS, valid_source, summarize, run,
)


def fixtures(root,ancestry):
    folder=root/f"is-GWAS-1000G-all22-{ancestry}"
    folder.mkdir(parents=True)
    by={}
    for c in range(1,23):
        by[str(c)]={
            "gwas_rows":120,"bim_rows":200,
            "chr_pos_overlap":110,
            "chr_pos_overlap_pct":round(100*110/120,4),
            "UNORDERED_ALLELE_SET_EXACT":90,
            "PALINDROMIC_AMBIGUOUS_NO_EAF_ORIENTATION":20,
            "COORDINATE_NOT_IN_LD_BIM":10,
        }
    d={
        "ancestry":ancestry,
        "status":"ALL_22_AUTOSOME_ORIGINAL_MD5_AND_BIM_OVERLAP_PASS",
        "all_22_complete":True,
        "previous_chr1_2_regression":{"status":"PASS"},
        "GWAS_checksum":{"publisher_MD5_verified":True},
        "LD_checksum":{"publisher_MD5_verified":True},
        "coloc_ready":False,
        "allele_effect_direction_harmonized":False,
        "LD_matrix_computed":False,
        "by_chr":by,
        "overall":{"all_22_autosome_gwas_rows":22*120,
                   "all_22_chr_pos_overlap":22*110}
    }
    p=folder/"all22_audit.json"
    p.write_text(json.dumps(d),encoding="utf-8")
    return p


class All22AncestrySynthesis(unittest.TestCase):
    def test_complete_reconciled_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            fixtures(root,"EUR")
            fixtures(root,"EAS")
            with patch.dict(EXPECTED_GWAS_ROWS,{"EUR":22*120,"EAS":22*120}):
                er,es,_=valid_source(root/"is-GWAS-1000G-all22-EUR"/"all22_audit.json","EUR")
                ar,as_,_=valid_source(root/"is-GWAS-1000G-all22-EAS"/"all22_audit.json","EAS")
                cmp=summarize(er,ar)
                self.assertEqual(len(cmp),22)
                self.assertEqual(es["rows"],22*120)
                self.assertEqual(as_["coordinate_found"],22*110)

    def test_outputs_no_etio_claims(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/"source";out=Path(tmp)/"out"
            fixtures(root,"EUR")
            fixtures(root,"EAS")
            with patch.dict(EXPECTED_GWAS_ROWS,{"EUR":22*120,"EAS":22*120}):
                report=run(root,out)
            self.assertEqual(report["coloc"],"NOT_RUN")
            self.assertEqual(report["MAGMA"],"NOT_RUN")
            self.assertEqual(report["association_effect_orientation"],"NOT_HARMONIZED")
            self.assertTrue((out/"all22_ancestry_by_chromosome.tsv").exists())

    def test_missing_eas_is_a_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fixtures(root,"EUR")
            with self.assertRaisesRegex(ValueError,"Missing true-GWAS"):
                valid_source(root/"is-GWAS-1000G-all22-EAS"/"all22_audit.json","EAS")

    def test_false_coloc_ready_is_a_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=fixtures(root,"EUR")
            d=json.loads(p.read_text(encoding="utf-8"));d["coloc_ready"]=True
            p.write_text(json.dumps(d),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"Unverified"):
                valid_source(p,"EUR")

    def test_missing_chromosome_is_a_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=fixtures(root,"EUR")
            d=json.loads(p.read_text(encoding="utf-8"));d["by_chr"].pop("22")
            p.write_text(json.dumps(d),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"Incomplete chromosome"):
                valid_source(p,"EUR")

    def test_inconsistent_alleles_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=fixtures(root,"EUR")
            d=json.loads(p.read_text(encoding="utf-8"))
            d["by_chr"]["5"]["PALINDROMIC_AMBIGUOUS_NO_EAF_ORIENTATION"]=40
            p.write_text(json.dumps(d),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"Unreconciled"):
                valid_source(p,"EUR")

    def test_reproducing_chr1_2_not_optional(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=fixtures(root,"EUR")
            d=json.loads(p.read_text(encoding="utf-8"))
            d["previous_chr1_2_regression"]["status"]="HOLD"
            p.write_text(json.dumps(d),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"Unverified"):
                valid_source(p,"EUR")

if __name__=="__main__":
    unittest.main()
