import tempfile
import unittest
from pathlib import Path
from scripts.ckm_stroke import pipeline

ROOT=Path(__file__).resolve().parents[1]
FIX=ROOT/"tests/fixtures/ckm_stroke"

class StagedTests(unittest.TestCase):
    def test_manifest_schema_but_science_hold(self):
        result=pipeline.p0(ROOT)
        self.assertEqual(result["status"],"PASS")
        self.assertEqual(result["scientific_status"],"HOLD")
        self.assertTrue(result["unresolved"])
    def test_normalized_gwas_fixture(self):
        self.assertEqual(pipeline.gwas_qc(FIX/"gwas.tsv")["status"],"PASS")
    def test_partial_is_never_complete(self):
        self.assertEqual(pipeline.gwas_qc(FIX/"gwas.tsv",1)["status"],"HOLD")
    def test_unapproved_ldsc_is_hold(self):
        self.assertEqual(pipeline.p2_to_p5("p2",ROOT,False,10)["status"],"HOLD")
    def test_no_real_replication_is_hold(self):
        self.assertEqual(pipeline.p6(ROOT)["status"],"HOLD")
    def test_math_fixture_not_science(self):
        with tempfile.TemporaryDirectory() as d:
            dest=Path(d)/"data/ckm_stroke"
            dest.mkdir(parents=True)
            (dest/"replication.csv").write_bytes((FIX/"replication.csv").read_bytes())
            result=pipeline.p6(Path(d))
            self.assertEqual(result["status"],"PASS")
            self.assertEqual(result["scientific_status"],"REVIEW_REQUIRED")
    def test_circular_exposure_rejected(self):
        # Review the condition used by the manifest gate.
        rows=pipeline.manifest(ROOT)
        rows[0]["phenotype"]="STROKE_CKM_STAGE"
        # The real gate always reads the current manifest and therefore
        # does not mutate persistent records in this test.
        self.assertEqual(pipeline.p0(ROOT,rows)["status"],"FAIL")

if __name__=="__main__":
    unittest.main()
