import json
import tempfile
import unittest
from pathlib import Path

from scripts.ckm_stroke import verify_glgc_previous_batch as gate


COMMIT="a"*40
def aggregate(chrom):
    traits={}
    for name in gate.TRAITS:
        traits[name]={
            "source_commit":COMMIT,
            "technical_status":"PASS",
            "checks":{k:True for k in gate.REQUIRED_CHECKS},
            "original_rows":110,
            "valid_snp_rows":100,
            "duplicate_keys":0,
        }
    return {"scope":gate.SCOPE,"chromosome":str(chrom),
            "status":"PASS","scientific_status":"HOLD",
            "source_commit":COMMIT,
            "missing_traits":[],"errors":[],
            "traits":traits,"snp_rows_total":300}


def populate(root):
    for chrom in (21,20,19):
        (root/f"glgc_chr{chrom}-aggregate.json").write_text(
            json.dumps(aggregate(chrom)),encoding="utf-8")


class TestPreviousBatchGate(unittest.TestCase):
    def test_three_prior_chromosomes_unlock_next_only(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            populate(folder)
            report=gate.verify_previous(folder,"21,20,19","18,17,16")
        self.assertEqual(report["status"],"PASS")
        self.assertEqual(report["scientific_status"],"HOLD")
        self.assertEqual(len(report["matrix"]["include"]),9)
        self.assertEqual(report["prior_chr_snp_totals"]["21"],300)

    def test_missing_prior_chromosome_blocks(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);populate(folder)
            (folder/"glgc_chr20-aggregate.json").unlink()
            with self.assertRaisesRegex(ValueError,"missing aggregate"):
                gate.verify_previous(folder,"21,20,19","18,17,16")

    def test_prior_hold_blocks(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);populate(folder)
            bad=aggregate(20);bad["status"]="HOLD"
            (folder/"glgc_chr20-aggregate.json").write_text(json.dumps(bad))
            with self.assertRaisesRegex(ValueError,"not technically verified"):
                gate.verify_previous(folder,"21,20,19","18,17,16")

    def test_mixed_source_commits_block(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);populate(folder)
            bad=aggregate(19)
            bad["source_commit"]="b"*40
            for record in bad["traits"].values():
                record["source_commit"]="b"*40
            (folder/"glgc_chr19-aggregate.json").write_text(json.dumps(bad))
            with self.assertRaisesRegex(ValueError,"mixed code revisions"):
                gate.verify_previous(folder,"21,20,19","18,17,16")

    def test_empty_technical_check_does_not_pass(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);populate(folder)
            bad=aggregate(21)
            bad["traits"]["TG"]["checks"]={}
            (folder/"glgc_chr21-aggregate.json").write_text(json.dumps(bad))
            with self.assertRaisesRegex(ValueError,"mandatory technical check"):
                gate.verify_previous(folder,"21,20,19","18,17,16")

    def test_skip_or_oversized_batch_is_blocked(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);populate(folder)
            for next_batch in ("17,16,15","18,18,17","18,17,16,15","22"):
                with self.assertRaises(ValueError):
                    gate.verify_previous(folder,"21,20,19",next_batch)

    def test_unexpected_artifact_blocks(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);populate(folder)
            (folder/"glgc_chr22-aggregate.json").write_text(json.dumps(aggregate(22)))
            with self.assertRaisesRegex(ValueError,"unexpected chromosome"):
                gate.verify_previous(folder,"21,20,19","18,17,16")


if __name__=="__main__":
    unittest.main()
