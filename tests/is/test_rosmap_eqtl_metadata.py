#!/usr/bin/env python3
"""Metadata contract tests; no restricted ROSMAP files are read."""
import sys
import unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "is"))
from audit_rosmap_eqtl_metadata import SOURCES, parse_entity  # noqa: E402


class RosmapMetadataTests(unittest.TestCase):
    def test_four_specific_public_entity_ids(self):
        self.assertEqual(set(SOURCES.values()), {
            "syn52363813", "syn52363851", "syn52363818", "syn52363834"
        })

    def test_file_entity_metadata(self):
        result = parse_entity(
            {"id": "syn52363851", "name": "subtype-eqtl-sumstats.Mic.2.tsv.gz",
             "dataFileHandleId": "1234", "contentSize": 123},
            "syn52363851"
        )
        self.assertEqual(result["synapse_id"], "syn52363851")
        self.assertTrue(result["file_handle_present"])
        self.assertEqual(result["content_size_bytes"], 123)

    def test_mismatched_entity_id_rejected(self):
        with self.assertRaises(ValueError):
            parse_entity({"id": "syn99999"}, "syn52363813")


if __name__ == "__main__":
    unittest.main()
