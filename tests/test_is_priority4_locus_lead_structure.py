import csv,unittest
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/audits/is_priority4_regional_leads_20261011_v1")
class TestLocusLeadStructure(unittest.TestCase):
 def test_two_locus_leads_not_four_independent_replications(self):
  p=ROOT/"IS_PRIORITY4_BBJ_GIGA_REGIONAL_LEAD_STRUCTURE.tsv"
  if not p.is_file():self.skipTest("Research server artifact unavailable")
  with p.open() as f: rr={x["locus"]:x for x in csv.DictReader(f,delimiter="\t")}
  self.assertEqual(set(rr),{"BBJ_IS_L001","BBJ_IS_L002"})
  a=rr["BBJ_IS_L001"];b=rr["BBJ_IS_L002"]
  self.assertEqual(a["same_lead_variant"],"True")
  self.assertEqual(b["same_lead_variant"],"False")
  self.assertEqual(a["common_nonpal_allele_exact_snps"],"2697")
  self.assertEqual(b["common_nonpal_allele_exact_snps"],"2626")
  self.assertLess(float(b["EAS504_r2_between_leads"]),.01)
  self.assertTrue(all(v["LD_proxy_not_in_study"]=="True" for v in rr.values()))
  self.assertTrue(all(v["independent_replication_not_proven"]=="True" for v in rr.values()))
if __name__=="__main__":unittest.main()
