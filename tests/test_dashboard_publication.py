"""Dashboard publication reads only explicit synthetic result files."""
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests/fixtures/operations"
SCRIPT = ROOT / "server/publish_research_dashboard.py"
SPEC = importlib.util.spec_from_file_location("dashboard_publication_test", SCRIPT)
dashboard = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(dashboard)


class DashboardPublicationTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.work = Path(self.directory.name)
        self.results = self.work / "results"
        self.repo = self.work / "dashboard"
        self.public = self.repo / "research-dashboard/public"
        self.public.mkdir(parents=True)

    def source(self, fixture, relative):
        target = self.results / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(FIXTURES / fixture, target)
        return target

    def load(self):
        return json.loads((self.public / "data/research-manifest.json").read_text())

    def test_missing_results_never_publish_historical_metrics_as_live(self):
        dashboard.run(self.results, self.repo, False)
        manifest = self.load()
        for study in manifest["studies"].values():
            self.assertEqual(study["status"], "pending")
            self.assertEqual(study["metrics"], [])
        self.assertEqual(manifest["results_root"], str(self.results))

    def test_separate_results_root_loads_real_fixture_metrics(self):
        self.source("ckd_susie.tsv", "ckd/stage2c/susie_results/STAGE2C_SUSIE_DEFAULT.tsv")
        checkout = self.work / "code/IS_Analysis_V3"
        checkout.mkdir(parents=True)
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--results-root", str(self.results),
             "--dashboard-repo", str(self.repo)], cwd=checkout, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        ckd = self.load()["studies"]["ckd"]
        self.assertEqual(ckd["status"], "live")
        self.assertEqual([metric["value"] for metric in ckd["metrics"]], ["1", "2"])

    def test_work_root_environment_and_legacy_source_resolve_to_results(self):
        self.source("ckd_susie.tsv", "ckd/STAGE2C_SUSIE_DEFAULT.tsv")
        with patch.dict(os.environ, {"IS_ANALYSIS_ROOT": str(self.work)}):
            subprocess.run([sys.executable, str(SCRIPT), "--dashboard-repo", str(self.repo)],
                           check=True, capture_output=True, text=True)
        self.assertEqual(self.load()["results_root"], str(self.results))
        subprocess.run([sys.executable, str(SCRIPT), "--source", str(self.work),
                        "--dashboard-repo", str(self.repo)], check=True, capture_output=True, text=True)
        self.assertEqual(self.load()["studies"]["ckd"]["status"], "live")

    def test_muscle_stage_reflects_loaded_input_instead_of_a_fixed_stage(self):
        self.source("muscle_stage1.json", "muscle/stage1/MUSCLE_STAGE1_SUMMARY.json")
        first = dashboard.muscle_manifest(self.results, self.public)
        self.assertEqual(first["status"], "stage1")
        self.source("muscle_stage2a.json", "muscle/stage2a/MUSCLE_STAGE2A_SUMMARY.json")
        second = dashboard.muscle_manifest(self.results, self.public)
        self.assertEqual(second["status"], "stage2a")
        self.assertIn({"label": "Variants with muscle sQTL", "value": "0"}, second["metrics"])

    def test_unreadable_or_wrong_shape_summary_is_pending(self):
        self.source("invalid_summary.json", "muscle/MUSCLE_STAGE2A_SUMMARY.json")
        self.source("invalid_summary.json", "ckd/STAGE4_PUBLIC_PROTOTYPE.json")
        self.assertEqual(dashboard.muscle_manifest(self.results, self.public)["status"], "pending")
        self.assertEqual(dashboard.ckd_manifest(self.results, self.public)["status"], "pending")

    def test_empty_existing_study_directory_does_not_imply_live_data(self):
        (self.results / "is").mkdir(parents=True)
        self.assertEqual(dashboard.generic_manifest(self.results, self.public, "is", "is")["status"], "pending")

    def test_missing_figures_do_not_make_git_add_reference_nonexistent_directory(self):
        completed = subprocess.CompletedProcess([], 0, "", "")
        with patch.object(dashboard.subprocess, "run", return_value=completed) as git, \
             patch("sys.stdout", new_callable=io.StringIO):
            dashboard.run(self.results, self.repo, True)
        add = git.call_args_list[0].args[0]
        self.assertEqual(add, ["git", "add", "research-dashboard/public/data/research-manifest.json"])


if __name__ == "__main__":
    unittest.main()
