"""Local-only GTEx ALT-major source and BBJ inverse-chain availability guards."""
import gzip,importlib.util,sys,tempfile,unittest
from pathlib import Path

P=Path(__file__).resolve().parents[1]/"scripts/is/audit_gtex_alt_major_10_snp_native_bbj.py"
spec=importlib.util.spec_from_file_location("eqtl_bbj_inverse",P)
m=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=m
spec.loader.exec_module(m)

class TestQTLBBJInverseChain(unittest.TestCase):
    def test_inverse_liftover_plus_one_based(self):
        blocks=[(111803000,111804000,112240804,1500,"test")]
        x=m.invert_one(blocks,111803962-1)
        self.assertEqual(x["pos37"],112241766)
        self.assertEqual(x["status"],"TOP_CHAIN_SELECTED_UNIQUE_POS")
    def test_wrong_position_unmapped(self):
        self.assertEqual(m.invert_one([],1234)["status"],"CHAIN_NO_MAPPING")
    def test_ambiguous_top_scoring_chain_fail_closed(self):
        blocks=[(100,200,1000,600,"a"),(100,200,2000,600,"b")]
        x=m.invert_one(blocks,150)
        self.assertEqual(x["status"],"CHAIN_AMBIGUOUS_TOP_SCORE")
        self.assertIsNone(x["pos37"])
    def test_parse_local_chain_block_gaps(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"test.chain.gz"
            # target=hg19, query=hg38. Query starts at 50 and source at 100.
            chain=("chain 500 chr12 1000 + 100 200 chr12 900 + 50 150 9\n"+
                   "30 10 20\n"+"40\n\n")
            with gzip.open(p,"wt") as f:f.write(chain)
            blocks=m.inverse_blocks(p,"chr12")
            self.assertEqual(len(blocks),2)
            self.assertEqual(m.invert_one(blocks,50)["pos37"],101)
            self.assertEqual(m.invert_one(blocks,100)["pos37"],141)
            self.assertEqual(m.invert_one(blocks,92)["status"],"CHAIN_NO_MAPPING")
if __name__=="__main__":
    unittest.main()
