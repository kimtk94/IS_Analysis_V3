"""Stage26 synthetic error contracts: no exact paper 5481 claim from aggregate sheet counts."""
import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from probe_ng_stage26_paper_genes import ENSG,FILES
from verify_ng_stage26_source_and_lineage import validate_lineage,validate_paper,SHA,COMMIT

class Stage26Tests(unittest.TestCase):
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
           "row_and_column_names_identical":True,
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
