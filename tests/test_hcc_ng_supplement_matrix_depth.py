"""Synthetic Ng XLSX layout tests; raw values and identifiers never released."""
import io,json,pathlib,sys,unittest,zipfile
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
from probe_ng_supplement_matrix_depth import (
    header_categories,parse_worksheets,inspect_file,FILES
)
XMLNS="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
def make_xlsx(bounds="A1:Z1000",first="Gene",second="RSEM TPM"):
    xml=(f'<worksheet xmlns="{XMLNS}"><dimension ref="{bounds}"/>'
         f'<sheetData><row r="1"><c r="A1" t="inlineStr"><is><t>{first}</t></is></c>'
         f'<c r="B1" t="inlineStr"><is><t>{second}</t></is></c></row>'
         f'<row r="2"><c r="A2"><v>2</v></c></row></sheetData></worksheet>')
    bio=io.BytesIO()
    with zipfile.ZipFile(bio,"w") as z:z.writestr("xl/worksheets/sheet1.xml",xml)
    return bio.getvalue()

class Tests(unittest.TestCase):
    def test_matrix_geometry_not_patient_join(self):
        x=inspect_file(make_xlsx())
        self.assertEqual(x["sheets_total"],1)
        a=x["sheet_summaries"][0]
        self.assertTrue(a["candidate_gene_by_sample_matrix_geometry"])
        self.assertIn("gene_feature",a["header_family_tags_allowlist"])
        self.assertIn("expression_measure",a["header_family_tags_allowlist"])
        self.assertFalse(a["sample_identifier_columns_verified"])
        self.assertFalse(x["same_patient_rna_protein_quantitative_join_verified"])
        self.assertNotIn("RSEM TPM",json.dumps(x))
    def test_narrow_summary_is_not_expression_matrix(self):
        r=inspect_file(make_xlsx("A1:I650"))["sheet_summaries"][0]
        self.assertEqual(r["excel_column_count"],9)
        self.assertFalse(r["candidate_gene_by_sample_matrix_geometry"])
    def test_invalid_archive_rejected(self):
        with self.assertRaises(ValueError):inspect_file(b"not an excel file")
    def test_source_registry_covers_supplements_2_to_14(self):
        self.assertEqual(len(FILES),14)
        self.assertIn("Supplementary_Data_13",FILES)
        self.assertNotIn("Supplementary_Data_1",FILES)
        self.assertEqual(FILES["Supplementary_Data_2"],"41467_2022_29960_MOESM5_ESM.xlsx")

if __name__=="__main__":unittest.main()
