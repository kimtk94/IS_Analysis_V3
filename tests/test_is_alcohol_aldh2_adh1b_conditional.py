"""Fail-closed synthetic unit tests for Japanese alcohol SuSiE source alignment."""
import csv,gzip,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from prepare_alcohol_aldh2_adh1b_conditional_susie import source_stats
from audit_alcohol_aldh2_adh1b_conditional_scientific_gates import audit
def tsv(path,rs):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter="\t")
        w.writeheader();w.writerows(rs)
class ConditionalAlcohol(unittest.TestCase):
    def test_exact_original_beta_se_alt_flip(self):
        with tempfile.TemporaryDirectory() as td:
            src=Path(td)/"study.gz"
            with gzip.open(src,"wt") as f:
                w=csv.DictWriter(f,fieldnames=["SNP","CHR","POS","EA","NEA",
                   "EAF","BETA","SE","P","HetP","N"],delimiter="\t")
                w.writeheader()
                w.writerow(dict(SNP="chr4_100239319_T_C",CHR=4,POS=100239319,
                   EA="C",NEA="T",EAF=.25,BETA=.4,SE=.1,P=1e-4,HetP=.1,N=140000))
                w.writerow(dict(SNP="chr12_112241766_G_A",CHR=12,POS=112241766,
                   EA="G",NEA="A",EAF=.75,BETA=.3,SE=.1,P=1e-4,HetP=.1,N=140000))
            observed=source_stats(src,{
               "4:100239319:T:C":{"reference_ALT":"C"},
               "12:112241766:G:A":{"reference_ALT":"A"}})
            self.assertAlmostEqual(observed["4:100239319:T:C"][0],4)
            self.assertAlmostEqual(observed["12:112241766:G:A"][0],-3)
    def fixture(self,root):
        for locus,s,blocked in (("ADH1B",.0133,False),("ALDH2",.0908,True)):
            p=root/locus;p.mkdir()
            (p/"input_qc.json").write_text(json.dumps({
               "status":"ALCOHOL_REGIONAL_REAL_GWAS_EAS_LD_INPUT_PASS",
               "source_beta_SE_per_SNP_verified":True,
               "LD_input_order_exact_variant_table":True,
               "reference_EAS_n":504,"variants":637 if locus=="ADH1B" else 574,
               "source_p_zero_count":442 if locus=="ALDH2" else 0}))
            (p/"ALCOHOL_SUSIE_DIAGNOSTIC_AND_GATE.json").write_text(json.dumps({
               "snps":637 if locus=="ADH1B" else 574,"reference_n":504,
               "LD_mismatch_s":s,"LD_gate":"BLOCKED_MODERATE_GWAS_LD_MISMATCH"
                 if blocked else "EXPLORATORY_MODEL_RUN_ALLOWED_NOT_CAUSAL",
               "converged":not blocked,"n_credible_sets":0 if blocked else 6,
               "final_finemapping_validated":False,
               "causal_alcohol_to_BP_AIS_mediation_estimated":False}))
        filterrows=[]
        for locus in ("ADH1B","ALDH2"):
            for scenario in ("SOURCE_QC_DEFAULT","REF_MAF_10_EAFDELTA_05",
              "REF_MAF_10_EAFDELTA_03","REF_MAF_05_EAFDELTA_03"):
                filterrows.append({"locus":locus,"scenario":scenario,
                  "LD_mismatch_s":.0133 if locus=="ADH1B" else .0908,
                  "causal_SNP_validated":"FALSE","causal_gene_validated":"FALSE"})
        tsv(root/"ALCOHOL_ADH1B_ALDH2_GWAS_LD_MISMATCH_FILTER_SENSITIVITY.tsv",filterrows)
        models=[]
        for name,lead,cs in (("BASELINE","4:100293600:C:T",6),
                             ("RESTRICTED","4:100243310:G:A",4)):
            for L in (3,6,10):
                models.append({"scenario":name,"L_maxcausal_assumed":L,
                  "SNP_count":637 if name=="BASELINE" else 261,
                  "LD_mismatch_s":.0133,"converged":"TRUE",
                  "top_pip_snp":lead,"exploratory_95pct_credible_sets":cs if L!=3 else 3,
                  "biological_target_or_mediation_validated":"FALSE",
                  "cross_model_independent_causal_signals_proven":"FALSE"})
        tsv(root/"ALCOHOL_ADH1B_SUSIE_MODEL_STABILITY.tsv",models)
    def test_invalid_causal_promotion_is_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);self.fixture(root)
            r=audit(root)
            self.assertEqual(r["status"],"ALDH2_BLOCKED_ADH1B_EXPLORATORY_UNSTABLE")
            self.assertFalse(r["ALDH2"]["credible_sets_run"])
            self.assertEqual(r["ADH1B"]["causal_credible_sets_validated"],0)
            p=root/"ALDH2"/"ALCOHOL_SUSIE_DIAGNOSTIC_AND_GATE.json"
            j=json.loads(p.read_text());j["final_finemapping_validated"]=True
            p.write_text(json.dumps(j))
            with self.assertRaisesRegex(ValueError,"Causal status"):
                audit(root)
    def test_missing_stability_scenario_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);self.fixture(root)
            path=root/"ALCOHOL_ADH1B_SUSIE_MODEL_STABILITY.tsv"
            with path.open() as f:listrows=list(csv.DictReader(f,delimiter="\t"))
            tsv(path,listrows[:-1])
            with self.assertRaisesRegex(ValueError,"Incomplete"):
                audit(root)
if __name__=="__main__":unittest.main()
