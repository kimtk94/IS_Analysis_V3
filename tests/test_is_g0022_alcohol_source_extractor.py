"""Official Koyanagi2024 Japanese alcohol source allele extraction; synthetic fixture."""
import csv,gzip,hashlib,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from extract_g0022_koyanagi_alcohol_rs671 import EXPECTED,read_source,extract,FILENAMES
def row(pos,ea,nea,beta,eaf):
    r={x:'.' for x in EXPECTED}
    r.update(SNP="rs671" if pos==112241766 else "rsX",
       CHR="12",POS=str(pos),EA=ea,NEA=nea,EAF=str(eaf),
       BETA=str(beta),SE=".02",P="1e-8",HetP=".21",N="145500")
    return r
def prep(folder,entries,source="alcohol_intake"):
    folder.mkdir(exist_ok=True,parents=True)
    filename=FILENAMES[source];file=folder/filename
    with gzip.open(file,"wt") as f:
        w=csv.DictWriter(f,fieldnames=list(EXPECTED),delimiter="\t")
        w.writeheader();w.writerows(entries)
    x={"key":filename,"size":file.stat().st_size,
       "checksum":"md5:"+hashlib.md5(file.read_bytes()).hexdigest()}
    (folder/"ZENODO_10038152_METADATA.json").write_text(json.dumps({"files":[x]}))
    return x
class AlcoholExtractionTests(unittest.TestCase):
    def test_ALT_effect_sign_and_frequency(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);file=prep(root,[row(112241766,"A","G",-.50,.25),
                       row(112168009,"G","A",.10,.70)])
            matching,n=read_source(root/FILENAMES["alcohol_intake"],file,minimum_rows=2)
            self.assertEqual(n,2)
            self.assertAlmostEqual(matching[("12","112241766")]["beta_ALT"],-.50)
            self.assertAlmostEqual(matching[("12","112168009")]["beta_ALT"],-.10)
            self.assertAlmostEqual(matching[("12","112168009")]["ALT_EAF"],.30)
            summary=extract(root,root/"out","alcohol_intake",minimum_rows=2)
            self.assertTrue(summary["rs671_present"])
            self.assertEqual(summary["4CS_GWAS_variants_harmonized"],2)
            self.assertFalse(summary["causal_stroke_alcohol_mediation_computed"])
    def test_partial_and_wrong_MD5_never_used(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);m=prep(root,[row(112241766,"A","G",-.50,.25)])
            f=root/FILENAMES["alcohol_intake"]
            with self.assertRaisesRegex(ValueError,"SIZE"):
                read_source(f,dict(m,size=f.stat().st_size+1),minimum_rows=1)
            with self.assertRaisesRegex(ValueError,"MD5"):
                read_source(f,dict(m,checksum="md5:"+"0"*32),minimum_rows=1)
    def test_allele_mismatch_and_duplicate_fail_closed(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);m=prep(root,[row(112241766,"T","C",-.1,.25)])
            with self.assertRaisesRegex(ValueError,"allele discordance"):
                read_source(root/FILENAMES["alcohol_intake"],m,1)
            m=prep(root,[row(112241766,"A","G",-.5,.25),
                row(112241766,"A","G",-.3,.25)])
            with self.assertRaisesRegex(ValueError,"Duplicate"):
                read_source(root/FILENAMES["alcohol_intake"],m,1)
if __name__=="__main__":unittest.main()
