"""Offline unit tests for GWAS build/strand/allele harmonization safety."""
import unittest
from scripts.ckm_stroke.tpmi_gigastroke_eas_genomewide_overlap import (
    allele_align,translate_snv,is_snv,unique_lift
)
class FakeLifter:
    def __init__(self,results):self.results=results
    def convert_coordinate(self,c,p):return self.results
class TestGenomewideOverlap(unittest.TestCase):
    def test_same_orientation(self):
        r=allele_align("T","C","C",0.21,"C","T",0.20)
        self.assertEqual(r[0],"MATCHED")
        self.assertEqual(r[1],1)
        self.assertAlmostEqual(r[2],0.01)
    def test_swapped_allele(self):
        r=allele_align("T","C","C",0.21,"T","C",0.8)
        self.assertEqual(r[0],"MATCHED")
        self.assertEqual(r[1],-1)
        self.assertAlmostEqual(r[2],0.01)
    def test_disallow_palindromic(self):
        self.assertEqual(allele_align("A","T","T",0.45,"T","A",0.45)[0],
                         "PALINDROMIC_UNRESOLVED")
        self.assertEqual(allele_align("C","G","G",0.45,"G","C",0.45)[0],
                         "PALINDROMIC_UNRESOLVED")
    def test_disallow_mismatch(self):
        self.assertEqual(allele_align("A","G","G",0.4,"A","C",0.4)[0],
                         "ALLELES_MISMATCH")
    def test_negative_chain_complement(self):
        self.assertEqual(translate_snv("T","C","C","-"),("A","G","G"))
        self.assertEqual(allele_align("A","G","G",0.21,"G","A",0.20)[1],1)
    def test_non_snv_ambiguous(self):
        self.assertFalse(is_snv("AT","A","A"))
        self.assertFalse(is_snv("A","G","T"))
        self.assertFalse(is_snv("A","A","A"))
    def test_coordinate_base_conversion(self):
        x=FakeLifter([("chr19",15573378,"+",450)])
        self.assertEqual(unique_lift(x,"19",15462568),
                         ("UNIQUE",("19",15573379,"+")))
    def test_reject_multimapping(self):
        x=FakeLifter([("chr19",15573378,"+",450),("chr19",100,"-",99)])
        self.assertEqual(unique_lift(x,"19",15462568)[0],"MULTIMAPPED")
    def test_reject_transchrom_mapping(self):
        x=FakeLifter([("chr18",15573378,"+",450)])
        self.assertEqual(unique_lift(x,"19",15462568)[0],"TRANSCHROM_OR_INVALID")
if __name__=="__main__":unittest.main()
