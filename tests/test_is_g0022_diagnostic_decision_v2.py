"""Final IS diagnostic ledger rejects unsupported causal/fine-mapping promotion."""
import csv,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_g0022_diagnostic_decision_v2 import build

def write(path,rows):
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open("w") as f:
  w=csv.DictWriter(f,delimiter="\t",fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

class DecisionTest(unittest.TestCase):
 def test_five_windows_and_reject_unresolved_as_subset(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td)
   names=["AS1","AS2","AIS1","AIS2","AIS3"]
   neff=[dict(trait="AS",working_n_rounded="98294"),dict(trait="AIS",working_n_rounded="70474")]
   runs=[];gates=[];krig=[];filt=[];comparison=[];random=[]
   for name in names:
    trait="AS" if name.startswith("AS") else "AIS"
    runs.append(dict(signal=name,trait=trait,
      approximate_n_eff="98294" if trait=="AS" else "70474",
      converged="FALSE" if name=="AS1" else "TRUE",
      purity_filtered_cs_count=1,top_pip=0.5))
    gates.append(dict(signal=name,claim="NOT_FINE_MAPPING_VALIDATED",
      scientific_gate="BLOCKED_LD_MISMATCH" if trait=="AS" else "EXPLORATORY"))
    krig.append(dict(signal=name,flagged_loglr2_absz2=0))
    for i in range(3):
     filt.append(dict(signal=name,filter=str(i)))
    if trait=="AS":
     random.append(dict(window=name,full_s="0.2",
       matched_study_shared_s="0.008",shared_n="400",
       random_dropout_s_min="0.17",random_dropout_s_median="0.2",
       random_dropout_replicates=20,random_subset_s_le_matched=0))
     for t in ("AS","AIS"):
      comparison.append(dict(window=name,trait=t,estimated_s="0.008"))
   write(root/"G0022_GIGASTROKE_EAS_EFFECTIVE_N_APPROX.tsv",neff)
   write(root/"effective_n_sensitivity/G0022_SUSIE_EFFECTIVE_N_SENSITIVITY.tsv",runs)
   write(root/"G0022_SUSIE_PILOT_EVIDENCE_GATES.tsv",gates)
   write(root/"G0022_KRIGING_OUTLIER_SUMMARY.tsv",krig)
   write(root/"G0022_LD_MISMATCH_FILTER_SENSITIVITY.tsv",filt)
   write(root/"matched_as_ais_ld_sensitivity/G0022_SHARED_SNP_AS_AIS_LD_DIAGNOSTIC.tsv",comparison)
   pr=root/"matched_as_ais_ld_sensitivity/G0022_AS_MATCHED_SNP_VS_RANDOM_DROPOUT.tsv"
   write(pr,random)
   (root/"G0022_AS_ONLY_SNP_CANONICAL_SOURCE_SUMMARY.json").write_text(json.dumps({
      "as_only_marker_rows":78,
      "source_presence_status":{"NO_AIS_CANONICAL_POSITION":78},
      "as_only_gws_p_le_5e8":0}))
   report=build(root)
   self.assertEqual(report["AS_windows_blocked"],2)
   self.assertEqual(report["AIS_windows_exploratory"],3)
   self.assertEqual(report["fine_mapping_validated"],0)
   random[0]["random_subset_s_le_matched"]=1
   write(pr,random)
   with self.assertRaisesRegex(RuntimeError,"Random SNP deletion"):
    build(root)
if __name__=="__main__":unittest.main()
