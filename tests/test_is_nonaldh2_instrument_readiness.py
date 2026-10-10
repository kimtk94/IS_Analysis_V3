"""No automatic IV promotion from strong marginal alcohol GWAS F statistics."""
import csv,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_nonaldh2_instrument_readiness import audit,GENE
class InstrumentReadiness(unittest.TestCase):
    def test_strong_exposure_is_not_valid_mr(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            lines=[]
            for i,var in enumerate(GENE):
                available=(i!=3)
                lines.append({
                    "variant_grch37":var,"alcohol_beta_over_se_abs":10,
                    "ais_qc":"ALLELE_HARMONIZED" if available else "NOT_IN_EAS_AIS_SOURCE",
                    "alcohol_ALT_eaf":.25,"ais_ALT_eaf":.255 if available else "",
                    "ais_p":.9 if available else "",
                    "alcohol_beta_ALT":.3,"ais_beta_ALT":.01 if available else ""})
            with (root/"IS_RS671_NON_ALDH2_ALCOHOL_POSITIONAL_GWS_AIS_QC.tsv").open("w") as f:
                w=csv.DictWriter(f,fieldnames=list(lines[0]),delimiter="\t")
                w.writeheader();w.writerows(lines)
            r=audit(root)
            self.assertEqual(r["minimum_marginal_F"],100)
            self.assertEqual(r["valid_independent_instruments_verified"],0)
            self.assertFalse(r["causal_MR_outcome_calculated"])
            self.assertEqual(r["AIS_nominal_p_less_0_05"],0)
            lines[0]["ais_ALT_eaf"]=".75"
            with (root/"IS_RS671_NON_ALDH2_ALCOHOL_POSITIONAL_GWS_AIS_QC.tsv").open("w") as f:
                w=csv.DictWriter(f,fieldnames=list(lines[0]),delimiter="\t")
                w.writeheader();w.writerows(lines)
            with self.assertRaisesRegex(ValueError,"freq"):
                audit(root)
if __name__=="__main__":unittest.main()
