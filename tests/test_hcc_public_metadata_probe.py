"""Offline tests for bounded, public GEO supplement schema inspection."""
import gzip
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "hcc"))
from probe_public_metadata import inspect, SOURCES  # noqa: E402


class PublicMetadataProbeTests(unittest.TestCase):
    def test_schema_detection_not_patient_claim(self):
        data = b"sample\tpatient_id\trecurrence_event\trecurrence_time\nS1\tP1\t1\t50\nS2\tP1\t0\t70\n"
        report = inspect(data)
        self.assertEqual(report["row_count"], 2)
        self.assertEqual(report["column_count"], 4)
        self.assertTrue(any(x["name"] == "recurrence_event" for x in report["candidate_columns"]))
        self.assertNotIn("verified_unique_patients", report)

    def test_empty_metadata(self):
        self.assertEqual(inspect(b"")["row_count"], 0)

    def test_sources_small_public_supplements_only(self):
        self.assertEqual(len(SOURCES), 2)
        self.assertTrue(all(url.startswith("https://ftp.ncbi.nlm.nih.gov/geo/") for _, _, url in SOURCES))
        self.assertTrue(all(".txt.gz" in filename for _, filename, _ in SOURCES))


if __name__ == "__main__":
    unittest.main()
