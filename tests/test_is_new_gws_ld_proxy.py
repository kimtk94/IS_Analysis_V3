"""Regression: missing SNP IDs may never silently produce LD proxies."""
import csv
import tempfile
import unittest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from build_new_gws_ld_proxies import run

class NewGWSLDTest(unittest.TestCase):
    def test_dot_id_blocked_and_stable_id_computed(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)/"result";root.mkdir()
            ref=Path(td)/"ref";ref.mkdir()
            with (root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv").open("w") as f:
                w=csv.DictWriter(f,fieldnames=["group_id","lead_variant"],delimiter="\t")
                w.writeheader();w.writerow({"group_id":"IS_XDATA_G0004","lead_variant":"2:100:G:T"})
            prefix=ref/"IS_XDATA_G0004.stable"
            Path(str(prefix)+".pgen").write_bytes(b"stub")
            Path(str(prefix)+".psam").write_text("#IID\nHG1\n")
            pvar=Path(str(prefix)+".pvar")
            pvar.write_text("#CHROM\tPOS\tID\tREF\tALT\n2\t100\t.\tT\tG\n")
            with self.assertRaisesRegex(ValueError,"stable"):
                run(root,ref,"IS_XDATA_G0004")
            pvar.write_text("#CHROM\tPOS\tID\tREF\tALT\n2\t100\t2:100:T:G\tT\tG\n")
            out=root/"new_gws_ld_proxy";out.mkdir()
            samples=[f"S{i}" for i in range(60)]
            traw=out/"IS_XDATA_G0004_STABLE_DOSAGE.traw"
            with traw.open("w") as f:
                f.write("\t".join(["CHR","SNP","(C)M","POS","COUNTED","ALT"]+samples)+"\n")
                vals=[str(i%3) for i in range(60)]
                f.write("\t".join(["2","2:100:T:G","0","100","T","G"]+vals)+"\n")
                f.write("\t".join(["2","2:110:A:C","0","110","A","C"]+vals)+"\n")
            result=run(root,ref,"IS_XDATA_G0004")
            self.assertEqual(result["lead"],"2:100:T:G")
            self.assertEqual(result["n_proxy_r2_ge_0.8"],1)
            self.assertEqual(result["lead_allele_order"],"REVERSED_ALLELE_ORDER")
if __name__=="__main__":unittest.main()
