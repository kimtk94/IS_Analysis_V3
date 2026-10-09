"""Synthetic allele orientations, ambiguous SNPs, reference conflicts, stable IDs."""
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_gws_gwas_reference_overlap import allele_class,ref_index,pvar_path

class GWASReferenceAuditTest(unittest.TestCase):
    def test_reference_alt_orientation(self):
        ref=("A","G","2:100:A:G")
        self.assertEqual(allele_class(ref,None,"G","A","","",.3),
                         ("MATCH_ALT_EFFECT",1,.3))
        self.assertEqual(allele_class(ref,None,"A","G","","",.3),
                         ("MATCH_REF_EFFECT",-1,.7))
    def test_conflict_palindrome_and_strand_mismatch(self):
        r=("A","T","2:100:A:T")
        self.assertEqual(allele_class(r,None,"A","T","","",.48)[0],
                         "PALINDROMIC_REVIEW")
        self.assertEqual(allele_class(r,None,"A","T","T","A",.2)[0],
                         "SOURCE_REF_ALT_CONFLICT")
        self.assertEqual(allele_class(("A","G","x"),None,"T","C","","",.2)[0],
                         "REFERENCE_ALLELES_DIFFER")
    def test_missing_snp_ids_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp)
            old=base/"old"
            new=base/"new"
            new.mkdir()
            p=new/"G1.stable.pvar"
            p.write_text("#CHROM\tPOS\tID\tREF\tALT\n2\t100\t.\tA\tG\n")
            with self.assertRaisesRegex(ValueError,"not stable"):
                ref_index(["G1"],new,old)
            p.write_text("#CHROM\tPOS\tID\tREF\tALT\n2\t100\t2:100:A:G\tA\tG\n")
            idx=ref_index(["G1"],new,old)
            self.assertEqual(idx["G1"][("2",100)][0][2],"2:100:A:G")
if __name__=="__main__":
    unittest.main()
