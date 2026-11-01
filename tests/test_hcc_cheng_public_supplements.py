"""No data-leak and XLSX schema tests for Cheng 2023 public supplements."""
import io
import pathlib
import sys
import unittest
import zipfile
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
import probe_cheng_public_supplements as m

class Tests(unittest.TestCase):
    def test_publisher_allowlist(self):
        self.assertTrue(m.accept_url("https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bib/24/2/10.1093_bib_bbac601/1/"+m.FILE))
        self.assertFalse(m.accept_url("http://random.org/"+m.FILE))
        self.assertFalse(m.accept_url("https://evil.org/"+m.FILE))

    def test_link_discovery_does_not_persist_signature(self):
        html=b'<a href="https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bib/24/2/10.1093_bib_bbac601/1/table_s1-s11_bbac601.xlsx?Expires=1&amp;Signature=something">'
        c=m.discover_links(html)
        self.assertTrue(len(c)>=1)
        self.assertTrue(c[0].startswith("https://"))
        self.assertEqual(m.FILE, c[0].split("?")[0].split("/")[-1])

    def test_fake_xlsx_manifest_sheets(self):
        buf=io.BytesIO()
        w=b'<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheets><sheet name="Table S1" sheetId="1"/></sheets></workbook>'
        s=b'<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><dimension ref="A1:C3"/><sheetData><row r="1"/><row r="2"/><row r="3"/></sheetData></worksheet>'
        with zipfile.ZipFile(buf,"w") as z:
            z.writestr("xl/workbook.xml",w)
            z.writestr("xl/worksheets/sheet1.xml",s)
        v=m.xlsx_manifest(buf.getvalue())
        self.assertEqual(v["sheet_count"],1)
        self.assertEqual(v["sheets"][0]["name"],"Table S1")
        self.assertEqual(v["xml_sheet_schemas"][0]["dimension"],"A1:C3")
        self.assertEqual(v["xml_sheet_schemas"][0]["row_elements"],3)
        self.assertFalse(v["patient_ids_exposed"])

    def test_html_is_not_workbook(self):
        with self.assertRaises(ValueError):
            m.xlsx_manifest(b"<html>fake</html>")

if __name__=="__main__":
    unittest.main()
