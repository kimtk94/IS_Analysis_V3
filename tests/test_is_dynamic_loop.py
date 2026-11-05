"""Offline safety/contract checks for the standalone IS Harness dynamic loop."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/harness/is_dynamic_loop.py"
spec = importlib.util.spec_from_file_location("is_dynamic_loop", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

class ISDynamicLoopTests(unittest.TestCase):
    def test_cycle_mapping(self):
        self.assertEqual(module.review_path(4).name, "stage_024_autonomous_review.json")
        self.assertEqual(module.review_path(5).name, "stage_026_autonomous_review.json")
        self.assertEqual(module.work_record_path(5).name, "stage_025_autonomous.json")

    def test_result_requires_novel_metrics(self):
        prior = {"metric_keys": ["same"]}
        valid = {"status": "ANALYSIS_COMPLETE", "metrics": {"new": 1.2},
                 "limitations": ["Exploratory only"], "source_files": ["data/handoff.md"]}
        self.assertIsNone(module.validate_result(valid, prior))
        valid["metrics"] = {"same": 2.0}
        self.assertEqual(module.validate_result(valid, prior), "NO_NEW_METRIC_KEYS")

    def test_nonfinite_values_blocked(self):
        prior = {"metric_keys": []}
        result = {"status": "ANALYSIS_COMPLETE", "metrics": {"ratio": float("nan")},
                  "limitations": ["Data uncertain"], "source_files": ["data/handoff.md"]}
        self.assertEqual(module.validate_result(result, prior), "NONFINITE_OR_NONJSON_RESULT")

    def test_generated_source_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "experiment.py"
            path.write_text('import json\nfrom pathlib import Path\nprint("ok")\n')
            self.assertIsNone(module.source_gate(path))
            path.write_text('import subprocess\n')
            self.assertEqual(module.source_gate(path), "DISALLOWED_IMPORT")
            path.write_text('eval("2+2")\n')
            self.assertEqual(module.source_gate(path), "DISALLOWED_CALL")

if __name__ == "__main__":
    unittest.main()
