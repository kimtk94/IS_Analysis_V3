import json
import tempfile
import unittest
from pathlib import Path

from scripts.ckm_stroke import audit_glgc_autosomes as audit

SHA="a"*40
def chrom_row(chrom, commit=SHA):
    traits={}
    for trait in audit.TRAITS:
        traits[trait]={
            "source_commit":commit,
            "technical_status":"PASS",
            "checks":{key:True for key in audit.MUST_CHECK},
            "original_rows":110,
            "valid_snp_rows":100,
            "excluded_non_snp":7,
            "excluded_symbolic_or_ambiguous_non_snp":3,
            "excluded_other":0,
            "palindromic_ambiguous_af":1,
            "duplicate_keys":0,
        }
    return {
        "scope":audit.SCOPE,"chromosome":str(chrom),
        "source_commit":commit,"status":"PASS","scientific_status":"HOLD",
        "missing_traits":[],"errors":[],"traits":traits,
        "snp_rows_total":300,
    }

def write(root, chrom, report=None, filename=None):
    (root/(filename or f"glgc_chr{chrom}-aggregate.json")).write_text(
        json.dumps(report if report is not None else chrom_row(chrom)),
        encoding="utf-8")


class TestGLGCAutosomeCoverage(unittest.TestCase):
    def test_empty_is_truthful_partial(self):
        with tempfile.TemporaryDirectory() as d:
            result=audit.assess(d)
        self.assertEqual(result["status"],"PARTIAL")
        self.assertEqual(result["observed_chromosome_pass"],0)
        self.assertEqual(result["missing_chromosome_count"],22)
        self.assertEqual(result["scientific_status"],"HOLD")

    def test_two_pass_chromosomes_is_not_genome_wide(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            write(folder,22)
            write(folder,21)
            (folder/"batch_plan.json").write_text('{"status":"READY"}',encoding="utf-8")
            result=audit.assess(folder)
        self.assertEqual(result["status"],"PARTIAL")
        self.assertEqual(result["observed_chromosome_pass"],2)
        self.assertEqual(result["trait_chromosome_pass"],6)
        self.assertEqual(result["totals_verified_pass_only"]["retained_snp_rows"],600)
        self.assertEqual(result["missing_chromosome_count"],20)
        self.assertFalse(result["whole_source_checksum_verified"])
        self.assertIn("batch_plan.json",result["ignored_files"])

    def test_22_pass_is_technical_only(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            for chrom in range(1,23):
                write(folder,chrom)
            result=audit.assess(folder)
        self.assertEqual(result["status"],"PASS")
        self.assertEqual(result["observed_chromosome_pass"],22)
        self.assertEqual(result["trait_chromosome_pass"],66)
        self.assertEqual(result["scientific_status"],"HOLD")
        self.assertEqual(result["ldsc_ready"],"NOT_ASSESSED")
        self.assertEqual(result["mr_ready"],"NOT_ASSESSED")

    def test_failed_trait_forces_hold(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            bad=chrom_row(21)
            bad["traits"]["TG"]["checks"]["no_duplicate_variant_keys"]=False
            write(folder,21,bad)
            result=audit.assess(folder)
        self.assertEqual(result["status"],"HOLD")
        self.assertEqual(result["coverage"][20]["traits"]["TG"],"HOLD")

    def test_missing_trait_forces_hold(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            bad=chrom_row(22)
            del bad["traits"]["HDL_C"]
            write(folder,22,bad)
            result=audit.assess(folder)
        self.assertEqual(result["status"],"HOLD")

    def test_bad_total_forces_hold(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            bad=chrom_row(20)
            bad["snp_rows_total"]=100
            write(folder,20,bad)
            result=audit.assess(folder)
        self.assertEqual(result["status"],"HOLD")
        self.assertTrue(any("SNP row total" in err for err in result["errors"]))

    def test_duplicate_chromosome_evidence_not_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            write(folder,20)
            write(folder,20,filename="another_chr20.json")
            result=audit.assess(folder)
        self.assertEqual(result["status"],"HOLD")
        self.assertTrue(any("duplicate chromosome" in e for e in result["errors"]))

    def test_cross_chromosome_code_revision_is_reviewed(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            write(folder,22)
            write(folder,21,chrom_row(21,commit="b"*40))
            result=audit.assess(folder)
        self.assertEqual(result["status"],"PARTIAL")
        self.assertTrue(result["mixed_code_commits_review_required"])
        self.assertEqual(len(result["source_commits"]),2)

    def test_corrupt_json_is_hold(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            (folder/"bad.json").write_text("{nope}",encoding="utf-8")
            result=audit.assess(folder)
        self.assertEqual(result["status"],"HOLD")

if __name__=="__main__":
    unittest.main()
