"""Synthetic 504-EAS/104-JPT genotype source orientation + subsampling tests.

No real genotype files or external downloads in GitHub CI.
"""
import csv,json,sys,tempfile,unittest,hashlib
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from prepare_alcohol_jpt104_ld_sensitivity import build as jpt
from prepare_alcohol_eas_n104_reference_controls import run as downsample

def fixture(root):
    base=root/"result"
    base.mkdir()
    samples=[f"HG{i:05d}" for i in range(504)]
    sf=root/"EAS.samples.txt"
    sf.write_text("\n".join(samples)+"\n")
    panel=root/"panel.txt"
    with panel.open("w",newline="") as f:
        w=csv.writer(f,delimiter="\t")
        w.writerow(["sample","pop","super_pop","gender"])
        for i,k in enumerate(samples):
            w.writerow([k,"JPT" if i<104 else "CHB","EAS","female"])
    for locus,ch in (("ADH1B","4"),("ALDH2","12")):
        d=base/locus
        d.mkdir()
        rows=[]
        arr=[]
        for j in range(50):
            pos=10000000+j*1000
            snp=f"{ch}:{pos}:G:A"
            x=np.array([(i+j)%3 for i in range(504)],dtype=float)
            arr.append(x)
            rows.append({"ID":snp,"reference_ALT":"A",
                "reference_EAS_ALT_EAF":x.mean()/2,
                "Japanese_original_ALT_EAF":x.mean()/2,
                "source_z":5,"source_n":154559,
                "source_beta_ALT":.5,"source_se":.1,
                "pos":pos})
        with (d/"variants.tsv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
            w.writeheader();w.writerows(rows)
        A=np.stack(arr)
        XX=(A-A.mean(axis=1,keepdims=True))/A.std(axis=1,ddof=1,keepdims=True)
        R=XX@XX.T/(504-1)
        (d/"ld.f64.rowmajor").write_bytes(R.astype("<f8").tobytes())
        (d/"input_qc.json").write_text(json.dumps({
            "status":"ALCOHOL_REGIONAL_REAL_GWAS_EAS_LD_INPUT_PASS",
            "reference_EAS_n":504,"variants":4}))
        with (d/"reference_eas504.traw").open("w",newline="") as f:
            w=csv.writer(f,delimiter="\t")
            w.writerow(["CHR","SNP","(C)M","POS","COUNTED","ALT"]+
                 [f"{a}_{a}" for a in samples])
            for j,row in enumerate(rows):
                # First 2 variants use ref as counted, others use ALT.
                counted,other=("G","A") if j<2 else ("A","G")
                dosages=2-A[j] if counted=="G" else A[j]
                w.writerow([ch,row["ID"],0,row["pos"],counted,other]+
                           [int(x) for x in dosages])
    return base,sf,panel

class JptEasSubsample(unittest.TestCase):
    def test_exact_allele_signed_ld_rebuild_and_fixed_random_controls(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);base,samples,panel=fixture(root)
            res=jpt(base,panel,samples)
            self.assertEqual(set(res),{"ADH1B","ALDH2"})
            for locus in res:
                self.assertLessEqual(
                    res[locus]["max_reconstructed_existing_EAS_signed_LD_absolute_difference"],1e-12)
                self.assertEqual(res[locus]["JPT_polymorphic_snps"],50)
                self.assertEqual(res[locus]["JPT_n"],104)
            one=downsample(base,samples,reps=3,seed=20261010)
            for locus in one:
                self.assertEqual(len(one[locus]["draws"]),3)
                self.assertTrue(all(x["status"]=="COMPLETE" for x in one[locus]["draws"]))
            p=base/"ADH1B"/"EAS104_random_reference_controls"/"reference_draw_01.ld.f64.rowmajor"
            check=hashlib.sha256(p.read_bytes()).hexdigest()
            downsample(base,samples,reps=3,seed=20261010)
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),check)
    def test_unrecognized_panel_blocks(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);base,samples,panel=fixture(root)
            txt=panel.read_text().replace("\tJPT\tEAS\t", "\tCHB\tEAS\t",1)
            panel.write_text(txt)
            with self.assertRaisesRegex(ValueError,"Official"):
                jpt(base,panel,samples)
    def test_bad_existing_eas_ld_blocks(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);base,samples,panel=fixture(root)
            p=base/"ALDH2"/"ld.f64.rowmajor"
            values=np.fromfile(p,dtype="<f8")
            values[2]+=.1
            p.write_bytes(values.tobytes())
            with self.assertRaisesRegex(ValueError,"cannot be reconstructed"):
                jpt(base,panel,samples)
if __name__=="__main__":unittest.main()
