"""PMC Cloud 2026 metadata-first file source audit; no guessed URLs."""
import io
import json
import pathlib
import sys
import unittest
from unittest.mock import patch
import zipfile
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
import probe_pmc2026_cloud_proteomics as m

class CloudTests(unittest.TestCase):
    def test_accept_official_metadata_and_media_only(self):
        self.assertTrue(m.official("https://pmc-oa-opendata.s3.amazonaws.com/metadata/PMC10783603.1.json"))
        self.assertTrue(m.official("https://pmc-oa-opendata.s3.amazonaws.com/PMC10783603.1/mmc3.xlsx?md5=aabb"))
        self.assertFalse(m.official("https://evil.example.org/PMC10783603.1/mmc3.xlsx"))
        self.assertFalse(m.official("http://pmc-oa-opendata.s3.amazonaws.com/PMC10783603.1/mmc3.xlsx"))

    def test_metadata_true_file_and_metadata_no_file(self):
        data=json.dumps({"pmcid":m.PMCID,"doi":m.DOI,"version":1,"is_pmc_openaccess":True,
                         "license_code":"CC-BY","media_urls":[
                         "https://pmc-oa-opendata.s3.amazonaws.com/PMC10783603.1/mmc3.xlsx?md5=abc",
                         "https://pmc-oa-opendata.s3.amazonaws.com/PMC10783603.1/mmc2.xlsx",
                         "https://evil.com/mmc3.xlsx"]}).encode()
        v=m.parse_article_metadata(data)
        self.assertEqual(v["matching_file_count"],1)
        self.assertEqual(v["media_url_count"],3)
        self.assertNotIn("abc",str({k:x for k,x in v.items() if k!="candidate_urls_in_memory"}))
        w=m.parse_article_metadata(json.dumps({"pmcid":m.PMCID,"media_urls":[]}).encode())
        self.assertEqual(w["matching_file_count"],0)

    def test_official_pubmed_filename_with_mmc3_suffix(self):
        meta={"pmcid":m.PMCID,"doi":m.DOI,"version":1,
              "media_urls":["https://pmc-oa-opendata.s3.amazonaws.com/PMC10783603.1/1-s2.0-S2666379123013154-mmc3.xlsx?md5=abc",
                            "https://pmc-oa-opendata.s3.amazonaws.com/PMC10783603.1/1-s2.0-S2666379123013154-mmc2.xlsx?md5=def"]}
        result=m.parse_article_metadata(json.dumps(meta).encode())
        self.assertEqual(result["matching_file_count"],1)
        self.assertIn("1-s2.0-S2666379123013154-mmc3.xlsx",result["official_xlsx_media_filenames"])
        self.assertNotIn("md5=",str({k:v for k,v in result.items() if k!="candidate_urls_in_memory"}))

    def test_official_s3_uri_converts_to_approved_https(self):
        s3="s3://pmc-oa-opendata/PMC10783603.1/1-s2.0-S2666379123013154-mmc3.xlsx?md5=beef"
        expected="https://pmc-oa-opendata.s3.amazonaws.com/PMC10783603.1/1-s2.0-S2666379123013154-mmc3.xlsx"
        self.assertEqual(m.resolve_cloud_uri(s3),expected)
        meta={"pmcid":m.PMCID,"doi":m.DOI,"version":1,"media_urls":[s3]}
        x=m.parse_article_metadata(json.dumps(meta).encode())
        self.assertEqual(x["matching_file_count"],1)
        self.assertEqual(x["candidate_urls_in_memory"],[expected])
        self.assertIn("xlsx",x["media_file_extension_counts"]) if "media_file_extension_counts" in x else None
        with self.assertRaises(ValueError):
            m.resolve_cloud_uri("s3://untrusted-bucket/PMC10783603.1/mmc3.xlsx")

    def test_wrong_pmc_is_rejected(self):
        with self.assertRaises(ValueError):
            m.parse_article_metadata(b'{"pmcid":"PMC12345","media_urls":[]}')

    def test_workbook_is_structurally_checked_without_values(self):
        buf=io.BytesIO()
        with zipfile.ZipFile(buf,"w") as z:
            z.writestr("xl/workbook.xml",b'<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheets><sheet name="Protein profiles" sheetId="1"/></sheets></workbook>')
        r=m.xlsx_only_manifest(buf.getvalue())
        self.assertEqual(r["sheet_count"],1)
        self.assertFalse(r["patient_expression_values_or_raw_ids_logged"])

    def test_anonymous_sheet_numeric_structure_only(self):
        buf=io.BytesIO()
        wb=b'<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheets><sheet name="1. 6512 proteins" sheetId="1"/></sheets></workbook>'
        ws=(b'<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            b'<dimension ref="A1:E3"/><sheetData>'
            b'<row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c></row>'
            b'<row r="2"><c r="A2" t="s"><v>2</v></c><c r="B2"><v>1.5</v></c><c r="C2"><v>2.5</v></c></row>'
            b'</sheetData></worksheet>')
        with zipfile.ZipFile(buf,"w") as z:
            z.writestr("xl/workbook.xml",wb)
            z.writestr("xl/worksheets/sheet1.xml",ws)
        profile=m.xlsx_only_manifest(buf.getvalue())
        self.assertEqual(profile["sheet_count"],1)
        d=profile["worksheet_structures"][0]
        self.assertEqual(d["dimension_ref"],"A1:E3")
        self.assertEqual(d["row_summaries"][1]["numeric_cells"],2)
        self.assertTrue(profile["patient_expression_values_or_raw_ids_logged"] is False)

    def test_tn_patient_labels_and_unclassified_are_aggregate_only(self):
        self.assertEqual(m.classify_patient_label("P001_T"),("p001","T",1))
        self.assertEqual(m.classify_patient_label("P001_N"),("p001","N",1))
        self.assertEqual(m.classify_patient_label("T_002"),("002","T",0))
        self.assertEqual(m.classify_patient_label("unknown"),None)

    def test_realistic_sharedstring_header_hides_identifiers(self):
        buf=io.BytesIO()
        wb=b'<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheets><sheet name="1. 6512 proteins" sheetId="1"/></sheets></workbook>'
        ws=(b'<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            b'<dimension ref="A1:E3"/><sheetData>'
            b'<row r="2"><c r="A2" t="s"><v>0</v></c><c r="B2" t="s"><v>1</v></c><c r="C2" t="s"><v>2</v></c>'
            b'<c r="D2" t="s"><v>3</v></c><c r="E2" t="s"><v>4</v></c></row>'
            b'</sheetData></worksheet>')
        strings=b'<sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'+b''.join(
            b'<si><t>'+x+b'</t></si>' for x in (b"gene",b"protein",b"annotation",b"Patient007_T",b"Patient007_N"))+b'</sst>'
        with zipfile.ZipFile(buf,"w") as z:
            z.writestr("xl/workbook.xml",wb)
            z.writestr("xl/worksheets/sheet1.xml",ws)
            z.writestr("xl/sharedStrings.xml",strings)
        with zipfile.ZipFile(io.BytesIO(buf.getvalue())) as z:
            result=m.sample_header_audit(z,"xl/worksheets/sheet1.xml",2)
        self.assertEqual(result["candidate_sample_columns"],2)
        self.assertEqual(result["candidate_keys_with_T_and_N"],1)
        self.assertTrue(result["matches_expected_304_columns"])
        self.assertFalse(result["patient_ids_or_labels_logged"])
        self.assertFalse(result["tumor_normal_pairing_verified"])
        self.assertNotIn("Patient007",str(result))

    def test_metadata_http404_is_hold_not_missing_science(self):
        import urllib.error
        def missing(*args,**kwargs):raise urllib.error.HTTPError("https://pmc-oa-opendata.s3.amazonaws.com/",404,"not found",None,None)
        with patch.object(m,"bounded",side_effect=missing):
            report=m.audit()
        self.assertEqual(report["status"],"PMC_2026_CLOUD_NOT_VERIFIED")
        self.assertFalse(report["download_verified"])
        self.assertFalse(report["same_patient_rna_protein_join"])

if __name__=="__main__":unittest.main()
