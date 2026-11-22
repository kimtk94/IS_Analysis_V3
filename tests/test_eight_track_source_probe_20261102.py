"""CI-only unit tests for bounded source probes and scientific HOLD boundaries."""
import csv
import importlib.util
import io
import json
import pathlib
import tempfile
import unittest
from unittest.mock import patch

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "live_source", ROOT / "scripts/research/probe_eight_tracks_20261102.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class SourceGateTests(unittest.TestCase):
    def test_exactly_eight_named_tracks(self):
        self.assertEqual(set(mod.TRACKS),
                         {"is", "glaucoma", "hcc", "radiomics",
                          "ocular", "muscle", "oral", "skin"})

    def test_pdc_parent_mapping_uniqueness(self):
        rows = [
            {"aliquot_id": "a1", "sample_id": "s1", "case_id": "c1"},
            {"aliquot_id": "a2", "sample_id": "s1", "case_id": "c1"},
        ]
        mapping, stats = mod.validate_pdc(rows)
        self.assertEqual((len(mapping), stats["samples"], stats["cases"]), (2, 1, 1))

    def test_pdc_aliquot_conflict_fails(self):
        rows = [
            {"aliquot_id": "a1", "sample_id": "s1", "case_id": "c1"},
            {"aliquot_id": "a1", "sample_id": "s2", "case_id": "c2"},
        ]
        with self.assertRaises(mod.IntegrityFail):
            mod.validate_pdc(rows)

    def test_pdc_sample_conflict_fails(self):
        rows = [
            {"aliquot_id": "a1", "sample_id": "s1", "case_id": "c1"},
            {"aliquot_id": "a2", "sample_id": "s1", "case_id": "c2"},
        ]
        with self.assertRaises(mod.IntegrityFail):
            mod.validate_pdc(rows)

    def test_pdc_duplicate_aliquot_rows_fail(self):
        rows = [
            {"aliquot_id": "a1", "sample_id": "s1", "case_id": "c1"},
            {"aliquot_id": "a1", "sample_id": "s1", "case_id": "c1"},
        ]
        with self.assertRaises(mod.IntegrityFail):
            mod.validate_pdc(rows)

    def ena_text(self, duplicate=False):
        buf = io.StringIO()
        w = csv.DictWriter(buf, fieldnames=[
            "run_accession", "sample_accession", "sample_alias", "fastq_md5"],
            delimiter="\t")
        w.writeheader()
        for i in range(46):
            w.writerow({"run_accession": "SRR0" if duplicate else "SRR" + str(i),
                        "sample_accession": "SAMN" + str(i),
                        "sample_alias": "TM" + str(i) + " 48h",
                        "fastq_md5": "a" * 32 + ";" + "b" * 32})
        return buf.getvalue()

    def test_ena_bounded_run_identity(self):
        out = mod.ena_summary(self.ena_text())
        self.assertEqual(out["unique_runs"], 46)
        self.assertEqual(out["fastq_md5_entries"], 92)
        self.assertEqual(out["timepoint_unique_text_matches"], 46)

    def test_ena_duplicate_runs_fail(self):
        with self.assertRaises(mod.IntegrityFail):
            mod.ena_summary(self.ena_text(duplicate=True))

    def test_metadata_pass_never_promotes_scientific_gate(self):
        with tempfile.TemporaryDirectory() as d:
            with patch.dict(mod.PROBES, {"is": lambda: ("PASS_STUDY_ID_ONLY", {}, "LD missing")}):
                self.assertEqual(mod.run_one("is", pathlib.Path(d)), 0)
            result = json.loads((pathlib.Path(d) / "is_20261102.json").read_text())
            self.assertEqual(result["scientific_gate"], "HOLD")
            self.assertEqual(result["probe_status"], "PASS_STUDY_ID_ONLY")
            self.assertTrue(result["no_biological_inference"])

    def test_source_hold_is_recorded_without_fake_pass(self):
        def blocked():
            raise mod.Hold("MISSING_SOURCE")
        with tempfile.TemporaryDirectory() as d:
            with patch.dict(mod.PROBES, {"hcc": blocked}):
                self.assertEqual(mod.run_one("hcc", pathlib.Path(d)), 0)
            result = json.loads((pathlib.Path(d) / "hcc_20261102.json").read_text())
            self.assertEqual(result["probe_status"], "HOLD_SOURCE_OR_SCHEMA")
            self.assertEqual(result["scientific_gate"], "HOLD")

    def test_integrity_conflict_fails_ci(self):
        def conflict():
            raise mod.IntegrityFail("BAD_MAPPING")
        with tempfile.TemporaryDirectory() as d:
            with patch.dict(mod.PROBES, {"hcc": conflict}):
                self.assertEqual(mod.run_one("hcc", pathlib.Path(d)), 1)
            result = json.loads((pathlib.Path(d) / "hcc_20261102.json").read_text())
            self.assertEqual(result["probe_status"], "FAIL_SOURCE_INTEGRITY")


if __name__ == "__main__":
    unittest.main()
