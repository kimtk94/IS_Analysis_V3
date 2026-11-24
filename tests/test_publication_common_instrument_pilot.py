import math
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import publication_common_instrument_pilot as m

class PilotTest(unittest.TestCase):
    def test_orientation(self):
        self.assertEqual(m.orientation("A","C","C","A"),(-1,"swap"))
        self.assertEqual(m.orientation("A","C","T","G"),(1,"complement"))
        self.assertIsNone(m.orientation("A","T","A","T"))
    def test_exact_ivw(self):
        x=[.1,.2,.3]
        r=m.fit(x,[.05,.1,.15],[.01]*3)
        self.assertAlmostEqual(r["beta"],.5)
        self.assertAlmostEqual(r["se"],.01/math.sqrt(.14))
        self.assertAlmostEqual(r["q"],0)
    def test_common_only(self):
        panel=[{"rsid":"rs"+str(i),"bx":.1,"f":40} for i in range(12)]
        out={k:{"rs"+str(i):{"by":.1*(j+1),"sy":.01} for i in range(j,12)} for j,k in enumerate(m.OUTCOMES)}
        r=m.analyze(panel,out)
        self.assertEqual(r["common"],10)
        self.assertEqual(r["family_n"],6)
        self.assertAlmostEqual(r["paired_contrasts"][0]["beta"],-1)
    def test_build_mapping_consistent_and_not_outcome_selected(self):
        panel=[{"rsid":"rs"+str(i),"chrom":i+1,"pos":100+i,"ea":"A","oa":"C"} for i in range(12)]
        def source(build,offset):
            return {r["rsid"]:{"mappings":[{"assembly_name":build,"seq_region_name":str(r["chrom"]),"start":r["pos"]+offset,"end":r["pos"]+offset,"allele_string":"A/C"}]} for r in panel}
        refs={"GRCh37":source("GRCh37",20),"GRCh38":source("GRCh38",0)}
        mapped,diag=m.map_sbp_panel(panel,refs)
        self.assertEqual(diag["source_build"],"GRCh38")
        self.assertEqual(len(mapped),12)
        self.assertEqual(mapped[0]["pos"],120)
        refs["GRCh37"]["rs0"]["mappings"][0].update(start=100,end=100)
        # Now the first source position is identical in both builds, not conflicting.
        mapped,diag=m.map_sbp_panel(panel,refs)
        self.assertEqual(diag["source_build"],"GRCh38")
    def test_missing_common(self):
        self.assertEqual(m.analyze([],{"LAA":{},"CE":{},"SVD":{}})["status"],"HOLD_TOO_FEW_COMMON_SNP")

if __name__=="__main__":unittest.main()
