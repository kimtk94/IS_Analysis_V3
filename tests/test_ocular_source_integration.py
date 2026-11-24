"""Offline checks for the extracted ocular source and explicit raw-data opt-in."""

import gzip
import importlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import yaml

try:
    import numpy as np
except ModuleNotFoundError as error:
    if error.name != "numpy":
        raise
    np = None


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts" / "ocular"
FIXTURES = ROOT / "tests" / "fixtures" / "ocular"


@unittest.skipIf(np is None, "NumPy is required for ocular source integration fixtures")
class OcularSourceIntegrationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Import the full dependency chain with all source network access forbidden.
        sys.path.insert(0, str(SCRIPTS))
        try:
            with mock.patch("urllib.request.urlopen", side_effect=AssertionError("offline test")):
                cls.full = importlib.import_module("ocular_p4_fullcohort")
                cls.pilot = importlib.import_module("ocular_p4_pseudobulk_pilot")
                cls.integrity = importlib.import_module("ocular_p1_integrity")
        finally:
            sys.path.remove(str(SCRIPTS))

    def setUp(self):
        self.network = mock.patch("urllib.request.urlopen", side_effect=AssertionError("offline test"))
        self.network.start()
        self.addCleanup(self.network.stop)
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.expected = json.loads((FIXTURES / "tiny_counts_expected.json").read_text())

    def compressed_fixture(self, filename, replacement=None):
        data = (FIXTURES / filename).read_bytes()
        if replacement:
            old, new = replacement
            self.assertIn(old, data)
            data = data.replace(old, new)
        path = self.directory / (filename + ".gz")
        path.write_bytes(gzip.compress(data, mtime=0))
        return path

    def stream(self, path):
        return self.full.matrix_stream(
            path,
            self.expected["sample_columns"],
            np.array(self.expected["group_for_column"], dtype=np.int32),
            self.expected["gene_count"],
            self.expected["group_count"],
            bytes_per_chunk=16,
        )

    def test_full_cohort_dependency_imports_are_complete(self):
        self.assertIs(self.full.matrix_stream, self.pilot.matrix_stream)
        self.assertIs(self.full.download, self.integrity.download)
        for name in ("ocular_p1_integrity.py", "ocular_p4_pseudobulk_pilot.py",
                     "ocular_p4_fullcohort.py", "ocular_gateb_edger_cutoff_sensitivity.py"):
            with self.subTest(script=name):
                result = subprocess.run(
                    [sys.executable, str(SCRIPTS / name), "--help"],
                    cwd=ROOT, capture_output=True, text=True, timeout=15,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("usage:", result.stdout)

    def test_fixture_counts_are_conserved_with_excluded_cells_and_chunk_boundaries(self):
        counts, metadata = self.stream(self.compressed_fixture("tiny_counts.mtx"))
        np.testing.assert_array_equal(counts, self.expected["pooled_counts"])
        for key in ("nnz", "selected_nnz", "selected_umis"):
            self.assertEqual(metadata[key], self.expected[key])
        self.assertEqual(metadata["gzip_crc"], "PASS")
        self.assertEqual(int(counts.sum()), self.expected["selected_umis"])

    def test_declared_nnz_mismatch_is_rejected(self):
        path = self.compressed_fixture("tiny_counts.mtx", (b"3 4 5\n", b"3 4 6\n"))
        with self.assertRaisesRegex(ValueError, "NNZ_MISMATCH"):
            self.stream(path)

    def test_out_of_range_cell_is_rejected(self):
        path = self.compressed_fixture("tiny_counts.mtx", (b"2 4 5\n", b"2 5 5\n"))
        with self.assertRaisesRegex(ValueError, "SPARSE_ENTRY_OUT_OF_RANGE"):
            self.stream(path)

    def test_small_metadata_fixture_preserves_cutoff_and_paired_donor_selection(self):
        path = self.compressed_fixture("tiny_metadata.tsv")
        for cutoff, expected in ((0.50, ["fixture_A", "fixture_B"]), (0.75, ["fixture_A"])):
            with self.subTest(cutoff=cutoff):
                cells, _ = self.pilot.read_metadata(path, cutoff)
                donors, _ = self.pilot.eligible_donors(cells, minimum_cells=1)
                self.assertEqual(donors, expected)
                self.assertEqual(set(cells["vehicle"].values()), set(expected))
                self.assertEqual(set(cells["dex"].values()), set(expected))

    def test_full_cohort_does_not_accept_a_tiny_fixture_as_actual_source(self):
        path = self.compressed_fixture("tiny_metadata.tsv")
        with self.assertRaisesRegex(ValueError, "UNEXPECTED_METADATA_COUNT_OR_LIBRARY_COUNT"):
            self.full.collect_metadata(path, 0.75)

class OcularMetadataAndWorkflowTest(unittest.TestCase):
    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is required for actual R parser/metadata evaluation")
    def test_p5_success_metadata_uses_the_selected_cutoff(self):
        result = subprocess.run(
            ["Rscript", str(ROOT / "tests" / "test_ocular_p5_summary.R"),
             str(SCRIPTS / "ocular_p5_edger.R"), str(FIXTURES / "cutoff_summary.tsv")],
            cwd=ROOT, capture_output=True, text=True, timeout=30,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("P5_SOURCE_PARSE_AND_CUTOFF_METADATA_PASS", result.stdout)

    def test_raw_data_job_requires_explicit_dispatch_and_offline_validation(self):
        path = ROOT / ".github" / "workflows" / "ocular-gateb-raw-score-sensitivity.yml"
        workflow = yaml.load(path.read_text(), Loader=yaml.BaseLoader)
        self.assertEqual(set(workflow["on"]), {"pull_request", "workflow_dispatch"})
        jobs = workflow["jobs"]
        raw = jobs["score-050-independent-refit"]
        self.assertEqual(raw["if"], "github.event_name == 'workflow_dispatch'")
        self.assertEqual(raw["needs"], "offline-source-validation")
        offline = jobs["offline-source-validation"]
        offline_commands = "\n".join(step.get("run", "") for step in offline["steps"])
        self.assertIn("test_ocular_source_integration.py", offline_commands)
        self.assertIn("numpy==1.26.4", offline_commands)
        self.assertIn("command -v Rscript", offline_commands)
        for job_name, job in jobs.items():
            commands = "\n".join(step.get("run", "") for step in job["steps"])
            if any(marker in commands for marker in (
                "python -u scripts/ocular/ocular_p4_fullcohort.py",
                "Rscript scripts/ocular/ocular_p5_edger.R",
                "gh run download",
            )):
                self.assertEqual(job_name, "score-050-independent-refit")
        self.assertNotIn("GH_TOKEN", json.dumps(offline))
        paths = workflow["on"]["pull_request"]["paths"]
        for name in ("ocular_p1_integrity.py", "ocular_p4_pseudobulk_pilot.py",
                     "ocular_p4_fullcohort.py", "ocular_p5_edger.R"):
            self.assertIn("scripts/ocular/" + name, paths)


if __name__ == "__main__":
    unittest.main()
