"""Contract tests for RNA/protein aggregate numeric provenance."""
import copy,pathlib,sys,unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
from assess_ng_rds_quantitative_qc import check

def good():
    files={}
    for name in ("HCCfull.rds","HCCfullz.rds","HCCdata.rds"):
        layers={}
        for k in ("T","P","PP"):
            n=150 if name=="HCCdata.rds" else 15000
            layers[k]={"feature_count":n,"sample_count":50,
                       "finite_numeric_count":n*50-2,"missing_or_nonfinite_count":2,
                       "negative_numeric_count":5,"zero_numeric_count":5,
                       "value_quantiles_global":{"p05":-3,"p25":0,"p50":1,"p75":2,"p95":3},
                       "feature_label_style_counts":{"ensembl_gene":n-1,"uniprot_swissprot":0,"nonempty":n}}
        files[name]={"status":"PUBLIC_MULTIMODAL_RDS_SHAPES_AUDITED",
                     "layers":layers,
                     "transcript_protein_shared_feature_label_count":140}
    return {"source_commit":"438f52271d0a86a6701baa7a96e80c516ba575c2",
            "objects":files}

class Tests(unittest.TestCase):
    def test_valid_summary_does_not_allow_inference(self):
        r=check(good())
        self.assertTrue(r["sources_from_real_rds"])
        self.assertFalse(r["quantitative_scale_reviewed_against_S5_appendix"])
        self.assertEqual(r["V02"],"PARTIAL")
        self.assertTrue(r["RNA_PROTEIN_EFFECTS_NOT_EVALUATED"])

    def test_mismatched_missingness_is_fail_closed(self):
        x=good();x["objects"]["HCCfull.rds"]["layers"]["T"]["missing_or_nonfinite_count"]=3
        with self.assertRaises(ValueError):check(x)

    def test_cannot_promote_unverifiable_feature_mapping(self):
        x=good();x["objects"]["HCCfull.rds"]["transcript_protein_shared_feature_label_count"]=None
        with self.assertRaises(ValueError):check(x)

    def test_wrong_provenance_is_rejected(self):
        x=good();x["source_commit"]="untrusted"
        with self.assertRaises(ValueError):check(x)

if __name__=="__main__":unittest.main()
