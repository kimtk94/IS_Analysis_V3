"""Synthetic genome-wide non-ALDH2 alcohol GWS + AIS harmonization gate."""
import csv,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_nonaldh2_alcohol_ais_sentinels import candidate,select_bins,match_ais

def alc(ch,pos,ref,alt,b=.2,se=.02,p=1e-12,n=145000,ea=None,nea=None):
    return {"SNP":f"chr{ch}_{pos}_{ref}_{alt}",
      "CHR":str(ch),"POS":str(pos),"EA":ea or alt,"NEA":nea or ref,
      "EAF":".3","BETA":str(b),"SE":str(se),"P":str(p),"N":str(n),
      "HetP":".50"}
def stroke(ch,pos,ea,oa,b=.04,se=.02,p=.06,eaf=.3):
    return {"chr":str(ch),"pos":str(pos),"effect_allele":ea,
      "other_allele":oa,"beta":str(b),"se":str(se),"p":str(p),"eaf":str(eaf)}
class NonALDH2Tests(unittest.TestCase):
    def test_only_gws_nonpalindromic_distant_snps(self):
        x=candidate(alc(4,100239319,"T","C"))
        self.assertEqual(x["variant_grch37"],"4:100239319:T:C")
        self.assertFalse(x["independent_genetic_IV_certified"] if "independent_genetic_IV_certified" in x else False)
        self.assertIsNone(candidate(alc(12,112241766,"G","A")))
        self.assertIsNone(candidate(alc(4,100239319,"A","T")))
        self.assertIsNone(candidate(alc(4,100239319,"T","C",p=.4)))
        self.assertIsNone(candidate(alc(4,100239319,"T","C",n=5000)))
    def test_positional_separation_is_not_ld(self):
        rows=[alc(4,100239319,"T","C",.4,.01),
           alc(4,100300000,"T","G",.2,.01),
           alc(4,39413780,"A","G",-.3,.01),
           alc(2,27730940,"T","C",.2,.01)]
        selected,qc=select_bins(rows)
        self.assertEqual(len(selected),3)
        self.assertTrue(all(x["independent_genetic_IV_certified"] is False for x in selected))
        self.assertTrue(all(x["pleiotropy_exclusion_restriction_certified"] is False for x in selected))
    def test_stroke_allele_orientation_and_exclusion(self):
        selected,_=select_bins([alc(4,100239319,"T","C",.4,.01)])
        ais=[stroke(4,100239319,"T","C",.1,.02,.02,.7)]
        rows,qc=match_ais(selected,ais)
        self.assertEqual(qc["ais_exact_allele_match"],1)
        self.assertAlmostEqual(rows[0]["ais_beta_ALT"],-.1)
        self.assertAlmostEqual(rows[0]["ais_ALT_eaf"],.3)
        self.assertFalse(rows[0]["between_study_cohort_independence_verified"])
        bad,_=match_ais(selected,[stroke(4,100239319,"A","G")])
        self.assertEqual(bad[0]["ais_qc"],"ALLELE_PAIR_DISCORDANT")
        missing,_=match_ais(selected,[])
        self.assertEqual(missing[0]["ais_qc"],"NOT_IN_EAS_AIS_SOURCE")
if __name__=="__main__":unittest.main()
