import json
import tempfile
import unittest
from pathlib import Path

from scripts.ckm_stroke import summarize_glgc_chr as summary


def make(trait, chromosome="22", status="PASS", duplicate=0, wrong_count=False):
    counts={"original_rows":3,"retained_autosomal_snp":2,
            "non_snp":1,"malformed_columns":0,"wrong_chromosome":0,
            "invalid_coordinate":0,"invalid_alleles":0,
            "invalid_effect_or_se":0,"p_underflow_zero":0,
            "invalid_p":0,"invalid_allele_frequency":0,
            "palindromic_ambiguous_af_0_42_to_0_58":1}
    if wrong_count:
        counts["original_rows"]=100
    return {
        "scope":summary.SCOPE,
        "source":"GLGC_Graham2021","ancestry":"EUR",
        "chromosome":chromosome,"trait":trait,
        "source_commit":"abcdef123",
        "status":status,"row_cap_reached":False,
        "header_evidence":{"status":"PASS","origin":"official_gzip_prefix"},
        "counts":counts,
        "variant_qc":{"status":"HOLD" if duplicate else "PASS","rows":2,
                      "duplicate_excess_rows":duplicate,
                      "p_consistency_sampled":0,
                      "p_consistency_mismatches":0}
    }


def write(root,trait,data):
    (root/(trait+".json")).write_text(json.dumps(data),encoding="utf-8")


class TestGLGCChromAggregate(unittest.TestCase):
    def test_three_traits_pass_only_technical_chromosome(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            for t in summary.TRAITS:
                write(folder,t,make(t))
            result=summary.summarize(folder,"22")
        self.assertEqual(result["status"],"PASS")
        self.assertEqual(result["scientific_status"],"HOLD")
        self.assertEqual(result["snp_rows_total"],6)
        self.assertEqual(result["source_gwas_whole_genome_qc"],"NOT_RUN")

    def test_missing_trait_keeps_hold(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            write(folder,"LDL_C",make("LDL_C"))
            result=summary.summarize(folder,"22")
        self.assertEqual(result["status"],"HOLD")
        self.assertIn("TG",result["missing_traits"])

    def test_duplicate_variant_keys_keep_hold(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            for t in summary.TRAITS:
                write(folder,t,make(t,duplicate=1 if t=="TG" else 0))
            result=summary.summarize(folder,"22")
        self.assertEqual(result["status"],"HOLD")

    def test_wrong_row_accounting_keep_hold(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            for t in summary.TRAITS:
                write(folder,t,make(t,wrong_count=t=="HDL_C"))
            result=summary.summarize(folder,"22")
        self.assertEqual(result["status"],"HOLD")

    def test_duplicated_source_artifact_is_not_silently_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            for t in summary.TRAITS:
                write(folder,t,make(t))
            (folder/"duplicate.json").write_text(json.dumps(make("LDL_C")),encoding="utf-8")
            result=summary.summarize(folder,"22")
        self.assertEqual(result["status"],"HOLD")
        self.assertTrue(any("duplicate artifact" in e for e in result["errors"]))

    def test_different_input_commits_are_not_aggregated(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            for t in summary.TRAITS:
                row=make(t)
                if t=="TG":row["source_commit"]="different"
                write(folder,t,row)
            result=summary.summarize(folder,"22")
        self.assertEqual(result["status"],"HOLD")
        self.assertTrue(any("commits" in e for e in result["errors"]))


if __name__=="__main__":
    unittest.main()
