"""GTEx v8 DAP-G -> GRCh37 allele-safe liftover; no false coloc claims."""
import sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_g0022_gtex_v8_dapg import variant_liftover,complement,stable,fetch_gene

class MockLiftOver:
    def __init__(self,outputs):self.outputs=outputs
    def convert_coordinate(self,chrom,pos):return self.outputs.get((chrom,pos),[])
class GTExConversionTest(unittest.TestCase):
    def test_liftover_coordinate_zero_based_and_alt_direction(self):
        lookup=MockLiftOver({("chr12",111675703):[("chr12",112133109,"+",900)]})
        val,status=variant_liftover("chr12_111675704_G_A_b38",lookup)
        self.assertEqual(status,"LIFTOVER_ONE_TO_ONE")
        self.assertEqual(val,("12:112133110:G:A","+"))
    def test_negative_strand_complement_and_bad_mapping(self):
        l=MockLiftOver({("chr12",110):[("chr12",200,"-",999)]})
        value,status=variant_liftover("chr12_111_A_G_b38",l)
        self.assertEqual(value,("12:201:T:C","-"))
        self.assertEqual(complement("ACGT"),"TGCA")
        self.assertEqual(variant_liftover("chr12_111_A_T_b37",l)[1],
                         "BAD_GTEX_VARIANT_ID")
        self.assertEqual(variant_liftover("chr12_111_AA_G_b38",l)[1],
                         "CHROM_OR_NON_SNP")
    def test_ambiguous_liftover_fails_closed(self):
        l=MockLiftOver({("chr12",100):[("chr12",200,"+",1),
                                      ("chr12",300,"+",1)]})
        self.assertEqual(variant_liftover("chr12_101_C_T_b38",l)[1],
                         "LIFTOVER_AMBIGUOUS_OR_UNMAPPED")
    def test_stable_ensembl_id_does_not_copy_version(self):
        self.assertEqual(stable("ENSG00000111275.8"),stable("ENSG00000111275.12"))
if __name__=="__main__":unittest.main()
