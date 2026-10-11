import csv,json,sys,unittest
from pathlib import Path

ROOT=Path("/srv/is-analysis/results/is/audits")
PRE=ROOT/"is_priority4_bbj_eas_cojo_readiness_20261011_v1"
REF=ROOT/"is_priority4_bbj_eas_cojo_reference_20261011_v1"
FINAL=ROOT/"is_priority4_bbj_eas_cojo_final_20261011_v1"
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from prepare_is_priority4_bbj_eas_cojo_readiness import LOCI,MIN_INFO,MIN_FREQ
from finalize_is_priority4_bbj_eas_cojo_reference_qc import TOL

class TestBBJEASCOJOReadiness(unittest.TestCase):
 def check_sources(self):
  if not FINAL.exists():self.skipTest("Frozen server COJO format-preflight not mounted")
 def test_authoritative_reference_gate(self):
  self.check_sources()
  x=json.loads((PRE/"IS_PRIORITY4_BBJ_EAS_COJO_PRELIMINARY_INPUT_MANIFEST.json").read_text())
  self.assertEqual(x["scanned_BBJ_rows"],13435531)
  self.assertEqual(len(x["regions"]),2)
  self.assertTrue(x["no_formal_COJO_executed"])
  self.assertTrue(x["BBJ_IS_N_not_guessed"])
  self.assertIn("BLOCKED",x["full_COJO_scientific_gate"])
  self.assertEqual(x["summary_INFO_threshold"],MIN_INFO)
  self.assertEqual(x["BBJ_MAF_threshold"],MIN_FREQ)
  self.assertEqual([r["BBJ_exact_nonpal_INFO0p8_MAF0p01_N_valid"] for r in x["regions"]],[3513,3759])
  self.assertTrue(all(r["reference_1000G_EAS_samples"]==504 for r in x["regions"]))
 def test_alleles_frequency_and_finalized_ma(self):
  self.check_sources()
  gate=json.loads((FINAL/"IS_PRIORITY4_COJO_EAS504_REFERENCE_MATCH_QC.json").read_text())
  self.assertTrue(gate["no_formal_COJO_run"])
  self.assertTrue(gate["GCTA_official_recommended_min_reference_gt_4000"])
  self.assertTrue(gate["IS_2225_positional_gene_universe_preserved"])
  self.assertEqual({x["locus"]:x["GCTA_format_QC_ready"] for x in gate["results"]},
                   {"BBJ_IS_L001":3374,"BBJ_IS_L002":3447})
  for locus,(chrom,lo,hi,lead) in LOCI.items():
   fp=FINAL/(locus+".EAS504.INFO0p8.MAF0p01.diffFreq0p2.ma")
   with fp.open() as f:
    rr=csv.DictReader(f,delimiter="\t")
    self.assertEqual(rr.fieldnames,["SNP","A1","A2","freq","b","se","p","N"])
    data=list(rr)
   self.assertEqual(len(data),3374 if chrom=="4" else 3447)
   self.assertEqual(sum(r["SNP"]==lead for r in data),1)
   for r in data:
    part=r["SNP"].split(":")
    self.assertEqual(part[0],chrom)
    self.assertEqual(r["A1"],part[3])
    self.assertEqual(r["A2"],part[2])
    self.assertGreaterEqual(float(r["freq"]),0.01)
    self.assertLessEqual(float(r["freq"]),0.99)
    self.assertGreater(float(r["se"]),0)
    self.assertEqual(float(r["N"]),174686)
    self.assertLessEqual(lo,int(part[1]))
    self.assertLessEqual(int(part[1]),hi)
 def test_reference_maf_difference_kept_below_max(self):
  self.check_sources()
  a=json.loads((FINAL/"IS_PRIORITY4_COJO_EAS504_REFERENCE_MATCH_QC.json").read_text())
  self.assertTrue(all(r["max_abs_EAS_BBJ_AF_delta"]<=TOL for r in a["results"]))
  self.assertEqual(sum(r["excluded_ancestry_AF_delta_gt_0p2"] for r in a["results"]),1)
if __name__=="__main__":unittest.main()
