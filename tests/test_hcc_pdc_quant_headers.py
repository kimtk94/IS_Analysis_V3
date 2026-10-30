"""Offline tests of bounded public PDC matrix header probes."""
import pathlib
import sys
import unittest
from unittest.mock import patch

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "hcc"))
from probe_pdc_quant_headers import (
    check_download_host, sniff_tsv, get_file_link, PDC_FILES, READ_LIMIT
)

class TestPDCQuantHeaders(unittest.TestCase):
    def test_reject_untrusted_or_http_url(self):
        for url in ("http://pdcdatastore.s3.amazonaws.com/a", "https://evil.example/f",
                    "https://pdcdatastore.s3.amazonaws.com.evil.example/x"):
            with self.assertRaises(ValueError):
                check_download_host(url)

    def test_valid_public_s3_url(self):
        self.assertTrue(check_download_host("https://pdcdatastore.s3.amazonaws.com/path"))

    def test_schema_no_patient_id_leak(self):
        row=b"Protein ID\tGene Symbol\tPatient123\tPatient456\nP123\tAPOE\t1\t2\n"
        o=sniff_tsv(row)
        self.assertTrue(o["header_complete_within_range"])
        self.assertEqual(o["header_delimiter_tabs"],3)
        self.assertNotIn("Patient123",str(o))
        self.assertEqual(o["first_25_columns_matching_biological_field_names"],["Protein ID","Gene Symbol"])

    def test_signed_url_kept_in_memory_only(self):
        name=PDC_FILES["PDC000198"][0]
        with patch("probe_pdc_quant_headers.request_graphql",return_value={
            "filesPerStudy":[{"file_name":name,"file_size":"134","file_format":"TSV",
                              "data_category":"Protein Assembly",
                              "signedUrl":[{"url":"https://pdcdatastore.s3.amazonaws.com/secret?Signature=abc"}]}]
        }):
            file,link=get_file_link("PDC000198",name)
        self.assertEqual(file["file_name"],name)
        self.assertTrue(link.startswith("https://"))
        self.assertEqual(READ_LIMIT,32768)

if __name__=="__main__":
    unittest.main()
