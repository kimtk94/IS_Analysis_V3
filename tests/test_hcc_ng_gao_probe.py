"""Offline bounded xlsx parser tests: do not need openpyxl or sample IDs."""
import io
import pathlib
import sys
import unittest
import zipfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "hcc"))
from probe_ng_gao_public import inspect_xlsx, _metadata_pair_summary, NG_URLS, PDC_QUERY


def make_xlsx(rows):
    ss = '<?xml version="1.0" encoding="UTF-8"?><sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">' + ''.join(
        '<si><t>' + str(s) + '</t></si>' for s in ["patient ID", "RNA", "Proteome"]) + '</sst>'
    xml = '<?xml version="1.0" encoding="UTF-8"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>' + ''.join(
        '<row r="' + str(i+1) + '">' + ''.join('<c t="s"><v>' + str(j) + '</v></c>' for j in row) + '</row>'
        for i, row in enumerate(rows)) + '</sheetData></worksheet>'
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w") as z:
        z.writestr("xl/sharedStrings.xml", ss)
        z.writestr("xl/worksheets/sheet1.xml", xml)
    return output.getvalue()


class HCCPublicPairedProbeTests(unittest.TestCase):
    def test_headers_can_be_inspected_without_patient_leaks(self):
        data = make_xlsx([[0, 1, 2], [0, 1, 2]])
        result = inspect_xlsx(data)
        self.assertEqual(result["worksheet_count"], 1)
        self.assertEqual(result["sheets"][0]["rows_in_xml"], 2)
        self.assertEqual(result["sheets"][0]["header_labels"], ["patient ID", "RNA", "Proteome"])
        self.assertFalse(result["sheets"][0]["expression_matrix_verified"])

    def test_metadata_pairing_counts_biopsies_and_unique_patients(self):
        rows = [
            {0:"tumor_biopsyID",1:"Patient_ID",2:"RNAseq",3:"Proteomics",4:"Phosphoproteomics"},
            {0:"T01",1:"P01",2:"Yes",3:"Yes",4:"Yes"},
            {0:"T02",1:"P01",2:"Yes",3:"Yes",4:"Yes"},
            {0:"T03",1:"P02",2:"Yes",3:"No",4:"No"},
            {0:"T04",1:"P03",2:"Yes",3:"Yes",4:"No"},
        ]
        result = _metadata_pair_summary(rows)
        self.assertEqual(result["rna_and_proteome_yes_biopsy_rows"], 3)
        self.assertEqual(result["rna_and_proteome_yes_unique_patient_ids"], 2)
        self.assertEqual(result["rna_and_proteome_and_phosphoproteome_yes_biopsy_rows"], 2)
        self.assertFalse(result["matrix_join_verified"])
        self.assertNotIn("P01", str(result))

    def test_non_xlsx_fails_closed(self):
        with self.assertRaises(ValueError):
            inspect_xlsx(b"not an xlsx")

    def test_only_public_https_sources(self):
        self.assertTrue(all(x.startswith("https://") for x in NG_URLS))
        self.assertIn("PDC000198", PDC_QUERY)
        self.assertIn("PDC000199", PDC_QUERY)


if __name__ == "__main__":
    unittest.main()
