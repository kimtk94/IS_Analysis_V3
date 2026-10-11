import csv,json,unittest
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/audits/is_priority4_unrestricted_gwas_leads_20261011_v1")
class TestRegionalFullGWAS(unittest.TestCase):
 def test_original_full_gwas_scan(self):
  p=ROOT/"IS_UNRESTRICTED_GWAS_LEADS_MANIFEST.json"
  if not p.is_file():self.skipTest("Research server output unavailable")
  d=json.loads(p.read_text())
  self.assertEqual(d["BBJ_full_source_rows"],13435531)
  self.assertEqual(d["GIGA_AIS_EAS_full_source_rows"],6725216)
  self.assertTrue(d["QTL_not_used_for_lead_selection"])
  self.assertTrue(d["original_2225_IS_gene_universe_preserved"])
  s={x["chr"]:x for x in d["results"]}
  self.assertEqual(s["4"]["GIGA_minp_variant"],"4:81181072:T:A")
  self.assertEqual(s["10"]["GIGA_minp_variant"],"10:105599770:G:T")
  self.assertEqual(s["4"]["BBJ_minp_variant"],"4:81182554:T:C")
  self.assertEqual(s["10"]["BBJ_minp_variant"],"10:105459834:T:C")
  self.assertTrue(all(x["LD_and_variant_alleles_recheck_required"] for x in s.values()))
