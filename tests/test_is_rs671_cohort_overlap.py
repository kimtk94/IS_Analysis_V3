"""Verify cohort-overlap audit never promotes overlap as independent MR cohorts."""
import csv,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_rs671_cohort_overlap import build
class Cohorts(unittest.TestCase):
    def test_japanese_exposure_BBJ_and_heldout_PGS_divergence(self):
        with tempfile.TemporaryDirectory() as td:
            r=build(Path(td))
            self.assertEqual(r["koyanagi_BBJ_n_as_published"],134993)
            self.assertEqual(r["koyanagi_total_drinking_n"],175672)
            self.assertAlmostEqual(r["BBJ_share_of_drinking_sample_count_approx"],.768438,places=5)
            self.assertFalse(r["two_sample_MR_independence_attested"])
            self.assertTrue(r["gigastroke_BBJ_PGS_evaluation_held_out"])
            self.assertFalse(r["held_out_PGS_implies_GIGASTROKE_EAS_summary_excludes_all_BBJ"])
            with (Path(td)/"IS_RS671_PAIRWISE_COHORT_OVERLAP_GATE.tsv").open() as f:
                p=list(csv.DictReader(f,delimiter="\t"))
            self.assertEqual(len(p),7)
            shared=[x for x in p if x["published_membership_risk"]=="SHARED_BBJ_COHORT"]
            self.assertGreaterEqual(len(shared),2)
            self.assertTrue(all(x["source_cohorts_completely_independent_confirmed"]=="NO" for x in shared))
if __name__=="__main__":unittest.main()
