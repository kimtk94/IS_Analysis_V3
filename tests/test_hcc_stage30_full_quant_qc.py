"""Source Gao PDC full-matrix aggregate numerical QC and biological HOLD guards."""
import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage30_gao_direct_quant_full_matrix import audit_matrix,verify

def fixture(nrows=11):
    header=["Gene/Aliquot"]+["PRIVATE_SAMPLE_"+str(i) for i in range(330)]
    return [header]+[["PRIVATE_GENE_"+str(j)]+
       [str(round((j-i)/12,5)) if (i%11) else "NA" for i in range(330)]
       for j in range(nrows)]
class MatrixTests(unittest.TestCase):
    def test_numeric_fraction_and_anonymized_labels(self):
        q=audit_matrix(fixture())
        self.assertEqual(q["quant_source_anonymous_sample_columns_n"],330)
        self.assertEqual(q["quant_source_matrix_row_n"],11)
        self.assertEqual(q["finite_numeric_cell_n"],11*300)
        self.assertEqual(q["nonfinite_missing_or_nonnumeric_cell_n"],11*30)
        self.assertFalse(q["gene_or_patient_names_or_individual_values_exported"])
        self.assertNotIn("PRIVATE_",str(q))
    def test_duplicate_header_fails(self):
        x=fixture();x[0][1]=x[0][2]
        with self.assertRaises(ValueError):audit_matrix(x)
    def test_missing_row_width_fails(self):
        x=fixture();x[1]=x[1][:-1]
        with self.assertRaises(ValueError):audit_matrix(x)
    def test_source_quant_ok_not_biological_claim(self):
        n=audit_matrix(fixture())
        studies=[{"pdc_study_id":p,"status":"FULL_PDC_DIRECT_QUANT_TABLE_NUMERIC_QC_SUCCEEDED",
          "unredacted_source_saved_to_public_artifacts":False,
          "same_patient_RNA_protein_join_verified":False,
          "quant_source_anonymous_sample_columns_n":330,
          "quant_source_matrix_row_n":n["quant_source_matrix_row_n"],
          "finite_numeric_cell_n":n["finite_numeric_cell_n"],
          "nonfinite_missing_or_nonnumeric_cell_n":n["nonfinite_missing_or_nonnumeric_cell_n"],
          "total_numeric_candidate_cells_n":n["total_numeric_candidate_cells_n"],
          "gene_or_patient_names_or_individual_values_exported":False,
          "complete_numeric_biological_sample_ID_join_verified":False}
          for p in ("PDC000198","PDC000199")]
        x={"schema":"hcc-stage30-PDC-full-quant-table-numeric-QC-v1","studies":studies}
        self.assertTrue(verify(x)["G30_D"].startswith("STOP"))
        x["studies"][0]["same_patient_RNA_protein_join_verified"]=True
        with self.assertRaises(ValueError):verify(x)
if __name__=="__main__":unittest.main()
