"""Fail-closed EAS/JPT LD evidence ledger; no automatic causal-IV promotion."""
import csv,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_nonaldh2_real_eas_ld_readiness import build,RS671
from extract_nonaldh2_alcohol_1kg_eas_ld import REGIONS

IDS=[f"{c}:{p}:{ref}:{alt}" for l in REGIONS.values() for c,p,ref,alt in l]
def tsv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader();w.writerows(rows)
class EvidenceGate(unittest.TestCase):
    def fixture(self,root):
        base=root/"result";ld=root/"ld";base.mkdir();ld.mkdir()
        pairs=[{"SNP_a":a,"SNP_b":b,"ALT_dosage_r2":.01}
           for i,a in enumerate(IDS) for b in IDS[i+1:]]
        tsv(ld/"IS_RS671_ALCOHOL_7SNP_REAL_EAS_LD_PAIRS.tsv",pairs)
        state={"reference_samples":504,"exact_verified_variants":7,
               "pairwise_r2_results":21,
               "ALT_EAS_reference_af":{v:{"n_nonmissing":504,
                    "n_missing":0,"1KG_EAS_ALT_allele_freq":.25} for v in IDS}}
        (ld/"IS_RS671_ALCOHOL_7SNP_REAL_EAS_LD_SUMMARY.json").write_text(json.dumps(state))
        for c in REGIONS:
            (ld/f"chr{c}.source_qc.json").write_text(json.dumps({
                "swapped_REF_ALT_source_count":1}))
        tsv(base/"IS_RS671_NON_ALDH2_ALCOHOL_POSITIONAL_GWS_AIS_QC.tsv",
            [{"variant_grch37":v,"alcohol_beta_ALT":-.2,
              "alcohol_se":.01,"alcohol_ALT_eaf":.24} for v in IDS if v!=RS671])
        tsv(base/"G0022_4CS_KOYANAGI2024_ALCOHOL_INTAKE_EXPOSURE.tsv",
            [{"variant_grch37":RS671,"beta_ALT":-1.3,"se":.01,"ALT_EAF":.25}])
        tsv(ld/"IS_RS671_7SNP_EAS_JPT_ALLELE_FREQUENCIES.tsv",
            [{"SNP":v,"JPT_ALT_EAF":.23} for v in IDS])
        (ld/"IS_RS671_7SNP_JPT_SENSITIVITY_SUMMARY.json").write_text(json.dumps({
              "sample_count_JPT":104,"pair_count":21,"max_measurable_JPT_r2":.04}))
        return base,ld,state
    def test_ld_qc_and_jpt_sources_are_not_causal_IVs(self):
        with tempfile.TemporaryDirectory() as td:
            base,ld,st=self.fixture(Path(td))
            r=build(base,ld)
            self.assertEqual(r["pairs_verified"],21)
            self.assertEqual(r["swapped_REF_ALT"],4)
            self.assertTrue(r["all_seven_pairwise_r2_below_0_1"])
            self.assertTrue(r["JPT_104_reference_validation_complete"])
            self.assertEqual(r["independent_genetic_instruments_causally_validated"],0)
    def test_missing_ld_pair_or_frequency_inference_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            base,ld,st=self.fixture(Path(td))
            st["pairwise_r2_results"]=20
            (ld/"IS_RS671_ALCOHOL_7SNP_REAL_EAS_LD_SUMMARY.json").write_text(json.dumps(st))
            with self.assertRaisesRegex(ValueError,"not complete"):
                build(base,ld)
            st["pairwise_r2_results"]=21
            st["ALT_EAS_reference_af"][IDS[0]]["1KG_EAS_ALT_allele_freq"]=.60
            (ld/"IS_RS671_ALCOHOL_7SNP_REAL_EAS_LD_SUMMARY.json").write_text(json.dumps(st))
            with self.assertRaisesRegex(ValueError,"10pp"):
                build(base,ld)
if __name__=="__main__":unittest.main()
