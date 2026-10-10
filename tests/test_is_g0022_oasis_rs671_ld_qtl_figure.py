"""Guard descriptive LD vs QTL figure against accidental coloc promotion."""
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from render_g0022_oasis_rs671_ld_qtl_evidence_svg import figure

class TestOasisTagLdFigure(unittest.TestCase):
    def rows(self):
        out=[]
        for gene,cell in [
            ("ALDH2","Mono-L1"),("ALDH2","B_Activated-L2"),
            ("BRAP","Mono-L1"),("BRAP","B_Activated-L2"),
            ("RPH3A","Mono-L1")]:
            for i in range(20):
                out.append(dict(gene=gene,cell=cell,
                    GWAS_SNP_grch37=f"12:{112000000+i}:A:G",
                    EAS504_r2_to_rs671=(i+1)/21,
                    qtl_p=10**(-1-i/3),
                    coloc_claim="NOT_VALID"))
        return out

    def test_svg_contains_explicit_no_coloc(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"study.svg"
            info=figure(self.rows(),p)
            ET.parse(p)
            self.assertEqual(info["groups"],5)
            self.assertIn("NOT colocalization",p.read_text())

    def test_block_causal_promotion(self):
        rows=self.rows()
        rows[0]["coloc_claim"]="VALID"
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError,"Causal claim"):
                figure(rows,Path(tmp)/"bad.svg")
    def test_wrong_source_group_count_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError,"five"):
                figure(self.rows()[:20],Path(tmp)/"bad.svg")

if __name__=="__main__":
    unittest.main()
