import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from plot_ckd_manuscript_fig3 import GENES,load,parse_h4,create

class FigureThreeTests(unittest.TestCase):
    def write_fixture(self,p):
        columns=("gene_symbol","eur_abf_h4","eur_susie_max_h4","eur_stage2_class","eur_comparison_status")
        with p.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=columns,delimiter="\t");w.writeheader()
            for g in GENES:
                w.writerow(dict(gene_symbol=g,eur_abf_h4="0.9",
                    eur_susie_max_h4="NA" if g=="UMOD" else "0.82",
                    eur_stage2_class="susie_unresolved" if g=="UMOD" else "core_shared_signal",
                    eur_comparison_status=""))
    def test_missing_is_not_zero(self):
        self.assertIsNone(parse_h4("NA","UMOD"))
        self.assertEqual(parse_h4("0","negative signal"),0)
    def test_invalid_pph4_rejected(self):
        with self.assertRaises(ValueError):parse_h4("1.2","bad")
    def test_nine_gene_contract(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"matrix.tsv";self.write_fixture(p)
            self.assertEqual(len(load(p)),9)
            p.write_text(p.read_text().replace("UMOD\t","ACP1\t"))
            with self.assertRaisesRegex(ValueError,"exactly 9"):load(p)
    @unittest.skipUnless(importlib.util.find_spec("matplotlib"),"matplotlib optional")
    def test_svg_and_tsv_and_no_overwrite(self):
        import hashlib
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=root/"matrix.tsv";self.write_fixture(p)
            sha=hashlib.sha256(p.read_bytes()).hexdigest()
            out=root/"output";report=create(p,out,sha)
            self.assertEqual(report["abf_nonmissing"],9)
            self.assertEqual(report["susie_nonmissing"],8)
            self.assertIn("<svg",(out/"CKD_FIG3_EUR_ABF_VS_SUSIE.svg").read_text())
            with (out/"CKD_FIG3_EVIDENCE.tsv").open() as f:
                genes=list(csv.DictReader(f,delimiter="\t"))
            self.assertEqual(genes[-1]["eur_susie_pph4"],"NA")
            with self.assertRaisesRegex(ValueError,"refusing overwrite"):create(p,out,sha)
            with self.assertRaisesRegex(ValueError,"SHA256 mismatch"):create(p,root/"another","0"*64)

if __name__=="__main__": unittest.main()
