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
    def test_missing_common(self):
        self.assertEqual(m.analyze([],{"LAA":{},"CE":{},"SVD":{}})["status"],"HOLD_TOO_FEW_COMMON_SNP")

if __name__=="__main__":unittest.main()
