"""G0022 Bayesian pilot gate does not promote high LD mismatch or unconverged fits."""
import csv,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_g0022_susie_pilot_gates import build

def save(path,rs):
    with path.open("w") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=list(rs[0]))
        w.writeheader();w.writerows(rs)

class TestPilotGate(unittest.TestCase):
 def test_block_mismatch_and_convergence_preserve_finemap_unresolved(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)
   inputs=[];fits=[];diag=[]
   for n in range(5):
    signal=f"S{n}"
    inputs.append(dict(signal=signal,status="PILOT_INPUT_QC_PASS",
         study="GCST90104544" if n==0 else "GCST90104545",
         trait="AS" if n==0 else "AIS",variant_count=300,reference_samples=504))
    for eff in (20000,256274):
     fits.append(dict(signal=signal,hypothetical_n=eff,converged="FALSE" if n==1 else "TRUE",
         cs_count=2,top_pip_variant="V1"))
     diag.append(dict(signal=signal,hypothetical_n=eff,
         ld_mismatch_s=.28 if n==0 else .02))
   (p/"G0022_SUSIE_PILOT_INPUT_MANIFEST.json").write_text(json.dumps(inputs))
   save(p/"G0022_SUSIE_RSS_PILOT_SUMMARY.tsv",fits)
   save(p/"G0022_GWAS_LD_CONSISTENCY.tsv",diag)
   summary=build(p)
   self.assertEqual(summary["final_fine_mapping_validated"],0)
   self.assertEqual(summary["high_ld_mismatch"],1)
   self.assertEqual(summary["both_n_converged"],4)
   self.assertEqual(summary["exploratory_low_mismatch_converged"],3)
   with (p/"G0022_SUSIE_PILOT_EVIDENCE_GATES.tsv").open() as f:
    rows={r["signal"]:r for r in csv.DictReader(f,delimiter="\t")}
   self.assertEqual(rows["S0"]["scientific_gate"],"BLOCKED_VERY_HIGH_LD_MISMATCH")
   self.assertEqual(rows["S1"]["scientific_gate"],"BLOCKED_NONCONVERGENCE")
   self.assertEqual(rows["S2"]["claim"],"NOT_FINE_MAPPING_VALIDATED")
if __name__=="__main__":unittest.main()
