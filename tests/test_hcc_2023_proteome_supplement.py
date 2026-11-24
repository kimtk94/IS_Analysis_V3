"""P0 external public 2023 HCC protein supplement ZIP manifest tests."""
import io
import pathlib
import sys
import unittest
import zipfile
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
import probe_2023_hcc_proteome_supplement as m
class Tests(unittest.TestCase):
    def test_ooxml_workbook(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,"w") as z:
            z.writestr("xl/workbook.xml",b'<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheets><sheet name="Proteomics" sheetId="1"/></sheets></workbook>')
            z.writestr("xl/worksheets/sheet1.xml",b"<worksheet/>")
        a=m.workbook_metadata(b.getvalue())
        self.assertEqual(a["sheet_count"],1)
        self.assertEqual(a["sheet_names"],["Proteomics"])
        self.assertFalse(a["values_or_patient_ids_logged"])
    def test_html_response_is_not_downloaded_spreadsheet(self):
        with self.assertRaises(ValueError):
            m.workbook_metadata(b"<html>not a real Excel workbook</html>")
    def test_repository_only_candidates(self):
        self.assertEqual(len(m.CANDIDATES),2)
        self.assertTrue(all("pmc.ncbi.nlm.nih.gov" in x for x in m.CANDIDATES))
if __name__=="__main__":
    unittest.main()
