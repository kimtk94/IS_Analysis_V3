"""Stage26 synthetic error contracts: no exact paper 5481 claim from aggregate sheet counts."""
import pathlib,sys,unittest,io,zipfile
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from probe_ng_stage26_paper_genes import ENSG,FILES,workbook_profile
from compare_ng_stage26_gene_membership import compute
from verify_ng_stage26_source_and_lineage import validate_lineage,validate_paper,SHA,COMMIT

class Stage26Tests(unittest.TestCase):
    def test_sheet_names_follow_relationships_not_sheet_id(self):
        wb=('<?xml version="1.0"?><workbook xmlns="http://schemas.openxmlformats.org/'
            'spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/'
            'officeDocument/2006/relationships"><sheets>'
            '<sheet name="Fig2b" sheetId="1" r:id="rId1"/>'
            '<sheet name="Fig2a" sheetId="2" r:id="rId2"/>'
            '</sheets></workbook>')
        rel=('<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/'
             'package/2006/relationships">'
             '<Relationship Id="rId1" Target="worksheets/sheet2.xml" Type="worksheet"/>'
             '<Relationship Id="rId2" Target="worksheets/sheet1.xml" Type="worksheet"/>'
             '</Relationships>')
        sheet=('<worksheet xmlns="http://schemas.openxmlformats.org/'
               'spreadsheetml/2006/main"><sheetData><row r="1">'
               '<c r="A1" t="inlineStr"><is><t>header</t></is></c>'
               '</row></sheetData></worksheet>')
        buf=io.BytesIO()
        with zipfile.ZipFile(buf,"w") as z:
            z.writestr("xl/workbook.xml",wb)
            z.writestr("xl/_rels/workbook.xml.rels",rel)
            z.writestr("xl/worksheets/sheet1.xml",sheet)
            z.writestr("xl/worksheets/sheet2.xml",sheet)
        sheets=workbook_profile(buf.getvalue())
        self.assertEqual([v["workbook_sheet_label"] for v in sheets],["Fig2a","Fig2b"])
    def test_exact_overlap_count_not_identical_paper_reproduction(self):
        orig={"ENSG00000141510","ENSG00000111111"}
        cur={"ENSG00000141510","ENSG00000222222"}
        q=compute(orig,cur)
        self.assertEqual(q["exact_Ensembl_gene_ID_overlap_n"],1)
        self.assertEqual(q["original_2022_genes_not_in_current_strict_n"],1)
        self.assertFalse(q["original_gene_set_semantics_reviewer_verified"])
        self.assertFalse(q["2022_Ensembl_vs_current_UniProt_mapping_versions_identical"])
    def test_original_paper_filename_provenance(self):
        self.assertTrue(FILES["publisher_Source_Data"].endswith("MOESM19_ESM.xlsx"))
        self.assertEqual(ENSG.findall("ENSG00000141510.12 | ENST01"),["ENSG00000141510.12"])
    def test_cannot_claim_exact_original_5481_by_count(self):
        sources={k:{"status":"UNVERIFIED"} for k in FILES}
        r={"schema":"hcc-stage26-ng-paper-gene-source-probe-v1",
           "paper_doi":"10.1038/s41467-022-29960-8","sources":sources,
           "exact_published_5481_gene_membership_verified":False,
           "patient_gene_and_protein_identifiers_exported":False,
           "raw_original_XLSX_in_artifact":False}
        self.assertEqual(validate_paper(r)["original_published_5481_genes_exact_membership"],"HOLD")
        r["exact_published_5481_gene_membership_verified"]=True
        with self.assertRaises(ValueError):validate_paper(r)
    def test_matrix_source_provenance_fails_on_wrong_sha(self):
        layer=lambda k:{"modality":k,"samples":50,"features":100,"both_finite_n":5000,
           "source_full_features_n":110,"source_fullz_features_n":100,
           "source_full_samples_n":50,"source_fullz_samples_n":50,
           "source_common_feature_n":100,"sample_set_identical":True,
           "feature_set_identical":False,
           "row_and_column_names_identical":False,
           "row_z_frac_abs_lt_1e_5":1.0,"col_z_frac_abs_lt_1e_5":0.01,
           "row_z_median_abs_residual":0.0,"col_z_median_abs_residual":1.0}
        r={"schema":"hcc-stage26-full-vs-fullz-lineage-v1",
           "source_commit":COMMIT,"full_sha256":SHA,"fullz_sha256":"f"*64,
           "layer_statistics":[layer(k) for k in ("M","CN","T","P","PP")],
           "original_ng_RSEM_TMM_edgeR_reproduced":False,
           "original_ng_5481_gene_universe_reproduced":False,
           "source_patient_and_gene_identifiers_exported":False}
        self.assertEqual(validate_lineage(r)["V03"],"PARTIAL")
        r["full_sha256"]="0"*64
        with self.assertRaises(ValueError):validate_lineage(r)
    def test_gene_row_candidate_does_not_imply_exact_identity(self):
        sources={k:{"status":"UNVERIFIED"} for k in FILES}
        sources["publisher_Source_Data"]={"status":"PUBLIC_XLSX_CONTENT_AUDITED",
          "source_url":"https://media.springernature.com/original/x",
          "sha256":"a"*64,"sheet_count_inspected":1,"sheets":[{
           "status":"SCAN_COMPLETED","unique_ENSG_like_identifiers_n":5481,
           "candidate_exact_5481_count_only":True}]}
        r={"schema":"hcc-stage26-ng-paper-gene-source-probe-v1",
          "paper_doi":"10.1038/s41467-022-29960-8","sources":sources,
          "exact_published_5481_gene_membership_verified":False,
          "patient_gene_and_protein_identifiers_exported":False,
          "raw_original_XLSX_in_artifact":False}
        self.assertEqual(validate_paper(r)["original_published_5481_genes_exact_membership"],"HOLD")
if __name__=="__main__":unittest.main()
