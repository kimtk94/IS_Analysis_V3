"""Synthetic contract tests for Ng/Suter raw-RDS aggregate auditing."""
import pathlib,sys,unittest,copy
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
from verify_ng_public_rds_evidence import verify,EXPECTED,MARKERS,SOURCE_COMMIT

def fixture():
    files={}
    objs={}
    for filename in EXPECTED:
        files[filename]={
            "status":"PUBLIC_PINNED_RDS_READ_IN_MEMORY",
            "sha256":"a"*64}
        objs[filename]={
            "status":"PUBLIC_MULTIMODAL_RDS_SHAPES_AUDITED",
            "transcript_protein_identical_sample_column_labels":True,
            "layers":{k:{"feature_count":25,"sample_columns":49,
                         "numeric_matrix":True,
                         "raw_labels_or_values_exported":False} for k in MARKERS}}
    return {
        "schema_version":"hcc-ng-public-rds-audit-v1",
        "source_commit":SOURCE_COMMIT,
        "raw_identifiers_or_numeric_matrices_in_artifact":False,
        "true_patient_by_gene_and_protein_quantitative_join_scientifically_validated":False,
        "required_three_multimodal_objects_read":True,
        "files":files,"objects":objs}

class Tests(unittest.TestCase):
    def test_readable_quantitative_layers_not_verified_patient_join(self):
        result=verify(fixture())
        self.assertEqual(result["execution"],"VERIFIED_PUBLIC_QUANTITATIVE_RDS_STRUCTURE")
        self.assertEqual(result["V02"],"BLOCKED_CLINICAL_SAMPLE_PROVENANCE_AND_QUANT_UNITS")
        self.assertFalse(result["biological_claim_ready"])
        self.assertEqual(result["files"]["HCCfull.rds"]["layers"]["T"]["sample_columns"],49)
        self.assertNotIn("raw_ids",str(result))

    def test_missing_rna_or_protein_layer_fails(self):
        v=fixture()
        v["objects"]["HCCfull.rds"]["layers"].pop("T")
        with self.assertRaises(ValueError):verify(v)

    def test_mismatch_of_sample_count_fails(self):
        v=fixture()
        v["objects"]["HCCfull.rds"]["layers"]["T"]["sample_columns"]=42
        result=verify(v)
        self.assertFalse(result["files"]["HCCfull.rds"]["rna_protein_equal_column_counts"])
        self.assertEqual(result["files"]["HCCfull.rds"]["status"],"SOURCE_GEOMETRY_OR_LABEL_MISMATCH_REQUIRES_JOIN")
        self.assertFalse(result["biological_claim_ready"])

    def test_unapproved_clinical_claim_fails(self):
        v=fixture()
        v["true_patient_by_gene_and_protein_quantitative_join_scientifically_validated"]=True
        with self.assertRaises(ValueError):verify(v)

    def test_source_commit_and_hash_are_required(self):
        v=fixture()
        v["source_commit"]="main"
        with self.assertRaises(ValueError):verify(v)
        v=fixture()
        v["files"]["HCCdata.rds"]["sha256"]="?"
        with self.assertRaises(ValueError):verify(v)
if __name__=="__main__":unittest.main()
