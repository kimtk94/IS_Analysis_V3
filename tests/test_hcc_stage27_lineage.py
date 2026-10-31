"""Stage27 source code and source-matrix rank lineages cannot imply replication."""
import copy,pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage27_public_code_lineage import summarize
from audit_ng_stage27_original_figure2_source import header_type,SHEETS
from verify_stage27_lineage import code_review,ranks_review,COMMIT,FULL,FULLZ

def synthetic_rank():
    layers=[]
    for k in ("M","CN","T","P","PP"):
        layers.append({
            "modality":k,"source_common_features_n":100,
            "biopsies_n":50,"rho_finite_features_n":95,
            "undefined_rho_feature_n":5,
            "spearman_rank_correlation_median":0.7,
            "spearman_rank_correlation_p10":0.1,
            "spearman_rank_correlation_p90":0.99,
            "spearman_fraction_ge_0_99":0.05,
            "spearman_fraction_lt_0_9":0.8,
            "spearman_fraction_negative":0.05,
            "pearson_median":0.6,
            "pearson_fraction_ge_0_99":0.02,
            "feature_or_sample_identifiers_exported":False,
        })
    return {
       "schema":"hcc-stage27-full-fullz-paired-rank-v1",
       "source_commit":COMMIT,
       "source_full_sha256":FULL,"source_fullz_sha256":FULLZ,
       "layer_metrics":layers,
       "exact_source_matrix_generation_reproduced":False,
       "original_Figure2b_gene_universe_reproduced":False,
       "original_Ng_2022_tumor_normal_relative_quantitation_reproduced":False,
       "gene_level_FDR_or_biological_significance_tested":False,
       "per_feature_rhos_or_patient_identifiers_exported":False}
class TestStage27(unittest.TestCase):
    def test_published_source_header_redaction(self):
        self.assertEqual(header_type("ENSG00000141510"),"ENS_IDENTIFIER_REDACTED")
        self.assertEqual(header_type("Spearman rho"),"CORRELATION_HEADER")
        self.assertEqual(header_type("P-value"),"SIGNIFICANCE_HEADER")
        self.assertEqual(header_type("TP53"),"UNCLASSIFIED_REDACTED")
        self.assertEqual(SHEETS,{"Fig2a","Fig2b","SuppFig13f","SuppFig13g"})
    def test_raw_code_loader_is_not_generator(self):
        s=summarize('HCCfullz<-readRDS("HCCinputs/HCCfullz.rds")\n')
        self.assertTrue(s["HCCfullz_consumed_in_clustering_script"])
        self.assertEqual(s["HCCfullz_file_writer_mentions"],0)
        self.assertEqual(s["saveRDS_call_mentions"],0)
    def test_raw_code_writer_detected(self):
        s=summarize('saveRDS(HCCfullz,"HCCinputs/HCCfullz.rds")\n')
        self.assertGreater(s["HCCfullz_file_writer_mentions"],0)
    def test_rank_contract_remains_partial(self):
        self.assertEqual(ranks_review(synthetic_rank())["V03"],"PARTIAL")
    def test_rank_missing_source_guard(self):
        d=synthetic_rank();d["source_fullz_sha256"]="0"*64
        with self.assertRaises(ValueError):ranks_review(d)
    def test_disallowed_biological_claim(self):
        d=synthetic_rank();d["original_Figure2b_gene_universe_reproduced"]=True
        with self.assertRaises(ValueError):ranks_review(d)
    def test_disallowed_feature_release(self):
        d=synthetic_rank();d["layer_metrics"][2]["feature_or_sample_identifiers_exported"]=True
        with self.assertRaises(ValueError):ranks_review(d)
    def test_bad_numerical_statistics(self):
        d=synthetic_rank();d["layer_metrics"][3]["spearman_rank_correlation_median"]=1.5
        with self.assertRaises(ValueError):ranks_review(d)
    def test_undefined_feature_count(self):
        d=synthetic_rank();d["layer_metrics"][2]["undefined_rho_feature_n"]=1
        with self.assertRaises(ValueError):ranks_review(d)
if __name__=="__main__":unittest.main()
