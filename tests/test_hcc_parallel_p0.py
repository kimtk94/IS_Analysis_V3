"""Offline regression: HCC parallel lanes, task isolation, fail-closed rollup."""
import json
import pathlib
import sys
import tempfile
import unittest
from types import SimpleNamespace
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from run_parallel_p0_track import ALL_TRACKS, TRACKS, run_track, tasks_for
from summarize_parallel_p0 import summarize, as_markdown

class Tests(unittest.TestCase):
    def test_six_parallel_lanes_have_sixteen_tasks(self):
        self.assertEqual(len(ALL_TRACKS),6)
        for name in ("ng-public","pdc-mapping","rna-access",
                     "geo-clinical","xing-clinical","regression"):
            self.assertIn(name,TRACKS)
        self.assertEqual(sum(len(tasks_for(t,pathlib.Path("/tmp"))) for t in ALL_TRACKS),16)

    def test_failure_isolation_inside_lane(self):
        with tempfile.TemporaryDirectory() as d:
            out=pathlib.Path(d)
            calls=[]
            def fake(command,**kwargs):
                i=len(calls)
                calls.append(command)
                (out/f"{TRACKS['rna-access'][i][0]}.json").write_text("{}")
                return SimpleNamespace(returncode=1 if i==1 else 0)
            rep=run_track("rna-access",out,runner=fake)
            self.assertEqual(len(calls),3)
            self.assertEqual(rep["successful_task_count"],2)
            self.assertEqual(rep["status"],"FAIL_EXECUTION")
            self.assertEqual(rep["tasks"][2]["result"],"PASS_EXECUTION")
            self.assertFalse(rep["biological_claim_supported"])
            self.assertTrue((out/"track_status.json").exists())

    def test_missing_json_fail_closed(self):
        with tempfile.TemporaryDirectory() as d:
            rep=run_track("geo-clinical",d,
                runner=lambda *a,**k:SimpleNamespace(returncode=0))
            self.assertEqual(rep["status"],"FAIL_EXECUTION")
            self.assertEqual(rep["tasks"][0]["error_class"],"missing_audit_output")

    def test_aggregator_guards_missing_or_unreviewed(self):
        with tempfile.TemporaryDirectory() as d:
            root=pathlib.Path(d)
            dest=root/"hcc-p0-xing-clinical"
            dest.mkdir()
            doc={"schema_version":"hcc-p0-parallel-v1",
                 "track":"xing-clinical","status":"PASS_EXECUTION",
                 "task_count":1,"tasks":[{"task":"a","result":"PASS_EXECUTION"}],
                 "biological_claim_supported":False}
            (dest/"track_status.json").write_text(json.dumps(doc))
            partial=summarize(root,expected=("xing-clinical","ng-public"))
            self.assertEqual(partial["overall_ci_status"],"FAIL_OR_MISSING")
            self.assertEqual(partial["lane_success_count"],1)
            self.assertIn("V02",as_markdown(partial))
            self.assertEqual(summarize(root,expected=("xing-clinical",))
                             ["overall_ci_status"],"PASS_EXECUTION")
            doc["biological_claim_supported"]=True
            (dest/"track_status.json").write_text(json.dumps(doc))
            self.assertEqual(summarize(root,expected=("xing-clinical",))
                             ["tracks"][0]["status"],"MALFORMED_ARTIFACT")

if __name__=="__main__":
    unittest.main()
