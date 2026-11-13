import json
import tempfile
import unittest
from pathlib import Path

from scripts.ckm_stroke import plan_glgc_batch as p


def gate():
    return {
        "status":"PASS","chromosome":"22","scope":p.SCOPE,
        "scientific_status":"HOLD","missing_traits":[],"errors":[],
        "traits":{t:{"technical_status":"PASS","original_rows":100,
                     "valid_snp_rows":90,"source_commit":"a"*40,
                     "checks":{"observed_header":True,"row_accounting":True}}
                  for t in p.TRAITS},
    }


class TestGLGCBatchPlan(unittest.TestCase):
    def test_chromosome_input_is_small_and_unique(self):
        self.assertEqual(p.chromosomes("21, 20, 19"),["21","20","19"])
        with self.assertRaises(ValueError): p.chromosomes("21,21")
        with self.assertRaises(ValueError): p.chromosomes("1,2,3,4")
        with self.assertRaises(ValueError): p.chromosomes("22")
        with self.assertRaises(ValueError): p.chromosomes("X")
        with self.assertRaises(ValueError): p.chromosomes("0")

    def test_real_chr22_gate_unlocks_only_three_by_three(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"gate.json"
            path.write_text(json.dumps(gate()),encoding="utf-8")
            out=p.plan("21,20,19",path)
            self.assertEqual(len(out["matrix"]["include"]),9)
            self.assertEqual(out["status"],"READY_FOR_TECHNICAL_QC_ONLY")
            self.assertEqual(out["scientific_status"],"NOT_RUN")
            self.assertEqual(out["required_chr22_gate_commit"],"a"*40)

    def test_missing_chromosome_evidence_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"gate.json"
            data=gate()
            del data["traits"]["HDL_C"]
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError,"3 GLGC traits"):
                p.plan("21,20,19",path)

    def test_failed_technical_gate_cannot_start_batch(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"gate.json"
            data=gate()
            data["status"]="HOLD"
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError,"PASS not established"):
                p.plan("21",path)

    def test_missing_source_commit_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"gate.json"
            data=gate()
            data["traits"]["TG"]["source_commit"]=""
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError,"source commit missing"):
                p.plan("21",path)


if __name__=="__main__":
    unittest.main()
