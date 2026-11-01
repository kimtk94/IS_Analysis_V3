"""Regression tests of the HCC multi-job source-evidence batch summary."""
import json
import pathlib
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "hcc"))
from aggregate_parallel_p0 import aggregate, to_markdown
from run_parallel_p0_track import TRACKS, tasks_for


def create_track(root, name, failed=(), damage=None):
    tasks = tasks_for(name, pathlib.Path("unused"))
    status = []
    for task, _, expected in tasks:
        bad = task in failed
        status.append({
            "task": task,
            "result": "FAIL_EXECUTION" if bad else "PASS_EXECUTION",
            "exit_code": 1 if bad else 0,
            "artifact_present": expected is not None and not bad,
            "seconds_elapsed": 0.03,
            "error_class": "nonzero_exit" if bad else None
        })
    doc = {
        "schema_version": "hcc-p0-parallel-v1",
        "track": name,
        "task_count": len(tasks),
        "executed_task_count": len(tasks),
        "successful_task_count": len(tasks)-len(failed),
        "status": "FAIL_EXECUTION" if failed else "PASS_EXECUTION",
        "tasks": status,
    }
    if damage == "bad_identity":
        doc["tasks"][0]["task"] = "unauthorized_task"
    if damage == "bad_count":
        doc["successful_task_count"] = -99
    d = root / f"hcc-p0-{name}"
    d.mkdir(parents=True)
    (d / "track_status.json").write_text(json.dumps(doc), encoding="utf-8")


class ParallelAggregationTests(unittest.TestCase):
    def test_all_tracks_need_full_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            for name in TRACKS:
                create_track(root, name)
            result = aggregate(root)
            self.assertTrue(result["all_execution_tracks_passed"])
            self.assertEqual(result["tracks_passed_execution"], len(TRACKS))
            self.assertEqual(result["tasks_expected"], sum(len(tasks_for(n, root)) for n in TRACKS))
            self.assertEqual(result["tasks_passed_execution"], result["tasks_expected"])
            self.assertFalse(result["biological_gate_pass_claimed"])
            self.assertEqual(result["V02_matched_rna_protein_quantitative_join"], "BLOCKED")
            self.assertIn("PASS_EXECUTION", to_markdown(result))

    def test_missing_lane_fails_even_if_all_returned_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            for name in list(TRACKS)[:-1]:
                create_track(root, name)
            result = aggregate(root)
            self.assertFalse(result["all_execution_tracks_passed"])
            self.assertEqual(result["tracks_failed_or_missing"], 1)
            self.assertEqual(result["tracks"][-1]["status"], "MISSING")

    def test_one_failed_task_keeps_other_tasks(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            for name in TRACKS:
                tasks=tasks_for(name,root)
                failed={tasks[0][0]} if name=="rna-access" else set()
                create_track(root, name, failed=failed)
            result=aggregate(root)
            self.assertFalse(result["all_execution_tracks_passed"])
            lane=next(x for x in result["tracks"] if x["track"]=="rna-access")
            self.assertEqual(lane["status"],"PARTIAL_EXECUTION_FAILURE")
            self.assertEqual(lane["successful_task_count"],len(TRACKS["rna-access"])-1)

    def test_tampered_task_name_or_count_fails_closed(self):
        for damage in ("bad_identity","bad_count"):
            with self.subTest(damage=damage),tempfile.TemporaryDirectory() as tmp:
                root=pathlib.Path(tmp)
                for name in TRACKS:
                    create_track(root,name,damage=damage if name=="ng-public" else None)
                result=aggregate(root)
                self.assertFalse(result["all_execution_tracks_passed"])
                self.assertEqual(result["tracks"][0]["status"],"INVALID_STATUS_FILE")

    def test_green_status_missing_claimed_output_is_invalid(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=pathlib.Path(tmp)
            for name in TRACKS:
                create_track(root,name)
            p=root/"hcc-p0-ng-public"/"track_status.json"
            data=json.loads(p.read_text())
            data["tasks"][0]["artifact_present"]=False
            p.write_text(json.dumps(data))
            result=aggregate(root)
            self.assertEqual(result["tracks"][0]["status"],"INVALID_STATUS_FILE")

if __name__=="__main__":
    unittest.main()
