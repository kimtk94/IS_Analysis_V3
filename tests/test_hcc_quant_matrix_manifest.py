"""Offline tests of safe HCC public quantitative-file manifest inventory."""
import io
import json
import pathlib
import sys
import unittest
import zipfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "hcc"))
from probe_quant_matrix_manifest import workbook_shape, read_bounded, NG_FILES, PDC_IDS, PRIDE_IDS


class HCCMatrixManifestTests(unittest.TestCase):
    def test_xlsx_worksheet_dimensions(self):
        data = ('<?xml version="1.0"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
                '<dimension ref="A1:AZ100"/><sheetData><row r="1"><c r="A1"><v>1</v></c></row>'
                '<row r="2"><c r="B2"><v>2</v></c></row></sheetData></worksheet>').encode()
        f = io.BytesIO()
        with zipfile.ZipFile(f, "w") as z:
            z.writestr("xl/worksheets/sheet1.xml", data)
        summary = workbook_shape(f.getvalue())
        self.assertEqual(summary[0]["xml_row_count"], 2)
        self.assertEqual(summary[0]["shape_ref"], "A1:AZ100")
        self.assertNotIn("sample_ids", summary[0])

    def test_invalid_file_fails_closed(self):
        with self.assertRaises(ValueError):
            workbook_shape(b"not-zip")

    def test_public_sources_only(self):
        self.assertEqual(set(PDC_IDS), {"PDC000198", "PDC000199"})
        self.assertEqual(set(PRIDE_IDS), {"PXD025705", "PXD025836"})
        self.assertEqual(len(NG_FILES), 4)


if __name__ == "__main__":
    unittest.main()
