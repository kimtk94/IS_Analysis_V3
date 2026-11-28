import csv,json,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts/is"))
import external_validation_qc as q
PACKET=ROOT/"docs/publication/external_validation_20261125"
class TestExternalAggregateQC(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.path=Path(self.tmp.name)/"synthetic.tsv"
        self.rows=q.read_tsv(PACKET/"aggregate_response_template.tsv",q.FIELDS,96)
        self.policy=Path(self.tmp.name)/"synthetic-policy.json"
        self.policy.write_text(json.dumps({"models":{"SYNTHETIC":{"phenotype":"LAS","imputation_metric":"INFO","min_quality":.8,"min_cases":100,"min_controls":100}}}))
    def save(self):
        with self.path.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=q.FIELDS,delimiter="\t");w.writeheader();w.writerows(self.rows)
    def ok(self,rs="rs7203740"):
        a=next(r for r in self.rows if r["rsid"]==rs and r["phenotype"]=="LAS")
        a.update(returned_build=a["requested_build"],returned_chr=a["requested_chr"],returned_pos=a["requested_pos"],returned_effect_allele=a["requested_effect_allele"],returned_other_allele=a["requested_other_allele"],beta="0.12",se="0.05",p_two_sided="0.0164",effect_allele_frequency="0.2",frequency_population="SYNTHETIC",imputation_metric="INFO",imputation_value="0.95",n_cases="400",n_controls="1600",per_marker_n="1999",model_id="SYNTHETIC",status="OK")
        return a
    def run_qc(self,policy=False):
        self.save();return q.run(PACKET,self.path,self.policy if policy else None)
    def issue(self,reason,policy=False):
        r,n=self.run_qc(policy);self.assertEqual(r["counts"]["ROW_QC_HOLD"],1);self.assertEqual(r["issues"][0]["reason_code"],reason);self.assertEqual(n,[])
    def test_blank_template(self):
        r,n=self.run_qc();self.assertEqual(r["counts"],{"NOT_RETURNED":96});self.assertEqual(n,[]);self.assertEqual(r["independent_validation_status"],"HOLD")
    def test_direct_without_policy(self):
        self.ok();r,n=self.run_qc();self.assertEqual(r["counts"]["STRUCTURAL_OK_POLICY_PENDING"],1);self.assertEqual(n[0]["beta_requested_allele"],.12)
    def test_reversed_exact_pair(self):
        a=self.ok();a["returned_effect_allele"],a["returned_other_allele"]=a["returned_other_allele"],a["returned_effect_allele"]
        r,n=self.run_qc(True);self.assertEqual(r["counts"]["STRUCTURAL_AND_CONFIGURED_QC_OK"],1);self.assertEqual(r["reversed_alleles"],1);self.assertAlmostEqual(n[0]["eaf_requested_allele"],.8);self.assertEqual(n[0]["beta_requested_allele"],-.12);self.assertEqual(n[0]["se"],.05);self.assertEqual(n[0]["p_two_sided"],.0164)
    def test_multiallelic_wrong(self):
        a=self.ok();a["returned_effect_allele"]="A";self.issue("RETURNED_ALLELE_MISMATCH")
    def test_coordinate(self):
        self.ok()["returned_pos"]="81775212";self.issue("RETURNED_COORDINATE_MISMATCH")
    def test_nonfinite(self):
        self.ok()["beta"]="NaN";self.issue("BETA_INVALID")
    def test_zero_p(self):
        self.ok()["p_two_sided"]="0";self.issue("P_INVALID_OR_UNDERFLOW")
    def test_invalid_se(self):
        self.ok()["se"]="0";self.issue("SE_INVALID")
    def test_invalid_n(self):
        self.ok()["per_marker_n"]="2001";self.issue("MARKER_N_EXCEEDS_TOTAL_OR_TOO_SMALL")
    def test_low_quality(self):
        self.ok()["imputation_value"]="0.75";self.issue("QUALITY_BELOW_PRELOCKED_THRESHOLD",True)
    def test_no_implicit_quality_threshold(self):
        self.ok()["imputation_value"]="0.75";r,_=self.run_qc();self.assertEqual(r["counts"]["STRUCTURAL_OK_POLICY_PENDING"],1)
    def test_unavailable(self):
        self.rows[0]["status"]="NOT_AVAILABLE";self.rows[0]["reason_if_unavailable"]="synthetic missing";r,n=self.run_qc();self.assertEqual(r["counts"]["NOT_AVAILABLE"],1);self.assertEqual(n,[])
    def test_unavailable_no_reason(self):
        self.rows[0]["status"]="QC_FAIL";self.issue("UNAVAILABLE_REASON_MISSING")
    def test_unreturned_data(self):
        self.rows[0]["returned_build"]="GRCh38";self.issue("UNRETURNED_WITH_DATA")
    def test_duplicate(self):
        self.rows[1]=dict(self.rows[0]);self.save()
        with self.assertRaisesRegex(ValueError,"UNEXPECTED_OR_DUPLICATE"):q.run(PACKET,self.path)
    def test_missing(self):
        self.rows.pop();self.save()
        with self.assertRaisesRegex(ValueError,"INCOMPLETE_FIXED"):q.run(PACKET,self.path)
    def test_request_changed(self):
        self.rows[0]["requested_pos"]="1";self.save()
        with self.assertRaisesRegex(ValueError,"REQUEST_IDENTITY_DRIFT"):q.run(PACKET,self.path)
    def test_missing_policy_model(self):
        self.ok()["model_id"]="UNKNOWN";r,_=self.run_qc(True);self.assertEqual(r["counts"]["STRUCTURAL_OK_POLICY_PENDING"],1)
if __name__=="__main__":unittest.main()
