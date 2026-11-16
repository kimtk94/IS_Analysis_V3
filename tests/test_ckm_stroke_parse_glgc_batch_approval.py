import copy
import unittest
from scripts.ckm_stroke import parse_glgc_batch_approval as p


def good():
    return {
        "schema":p.SCHEMA,"scientific_status":"HOLD",
        "prior_run_id":43865163287,
        "prior_chromosomes":["18","17","16"],
        "next_chromosomes":["15","14","13"],
        "expected_prior_source_commit":"a"*40,
        "prior_aggregate_statuses":{"18":"PASS","17":"PASS","16":"PASS"},
        "expected_prior_snp_rows":{"18":3649239,"17":3659008,"16":4237341},
        "max_parallel_jobs":3,"expected_next_matrix_jobs":9
    }


class TestGLGCApproval(unittest.TestCase):
    def test_valid_contiguous_source_provenance_plan(self):
        d=p.parse_approval(good())
        self.assertEqual(d["status"],"READY_FOR_REAL_EVIDENCE_GATE_ONLY")
        self.assertEqual(d["next_chromosomes"],"15,14,13")
        self.assertEqual(d["prior_chromosomes"],"18,17,16")
        self.assertEqual(d["scientific_status"],"HOLD")
        self.assertEqual(d["prior_run_id"],"43865163287")

    def test_skip_an_autosome_blocks(self):
        d=good();d["next_chromosomes"]=["14","13","12"]
        with self.assertRaises(ValueError):p.parse_approval(d)

    def test_duplicated_or_noncontiguous_input_blocks(self):
        for nums in (["15","15"],["15","13"],["15","14","13","12"],["22"]):
            d=good();d["next_chromosomes"]=nums
            d["expected_next_matrix_jobs"]=len(nums)*3
            with self.assertRaises(ValueError):p.parse_approval(d)

    def test_no_unverified_source_sha(self):
        d=good();d["expected_prior_source_commit"]="arbitrary"
        with self.assertRaisesRegex(ValueError,"40 hex"):p.parse_approval(d)

    def test_no_holds_or_partial_predecessor(self):
        d=good();d["prior_aggregate_statuses"]["16"]="HOLD"
        with self.assertRaisesRegex(ValueError,"HOLD"):p.parse_approval(d)

    def test_source_counts_exact(self):
        for counts in ({"18":3,"17":3},{"18":3,"17":0,"16":3},{"18":3,"17":True,"16":3}):
            d=good();d["expected_prior_snp_rows"]=counts
            with self.assertRaises(ValueError):p.parse_approval(d)

    def test_no_misleading_scientific_ready(self):
        d=good();d["scientific_status"]="PASS"
        with self.assertRaisesRegex(ValueError,"HOLD"):p.parse_approval(d)

    def test_no_extra_parallel_jobs(self):
        d=good();d["max_parallel_jobs"]=9
        with self.assertRaisesRegex(ValueError,"parallelism"):p.parse_approval(d)

    def test_must_match_matrix_count(self):
        d=good();d["expected_next_matrix_jobs"]=63
        with self.assertRaisesRegex(ValueError,"matrix"):p.parse_approval(d)

    def test_positive_numeric_run_id_not_bool(self):
        for run_id in (False,0,-1,"43865163287"):
            d=good();d["prior_run_id"]=run_id
            with self.assertRaises(ValueError):p.parse_approval(d)


if __name__=="__main__":
    unittest.main()
