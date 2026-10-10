"""1KG 504 EAS exact ref/alt and dosage flipping, EAS/JPT subset sensitivity."""
import csv,json,sys,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from extract_nonaldh2_alcohol_1kg_eas_ld import extract_chrom,get_plan,REGIONS,build_ld
from audit_nonaldh2_alcohol_jpt_ld_sensitivity import run as sensitivity

def make_sources(root,samples,style="balanced"):
    root.mkdir(parents=True,exist_ok=True)
    all_vars={c:[f"{chrom}:{pos}:{ref}:{alt}" for chrom,pos,ref,alt in v] for c,v in REGIONS.items()}
    for c,vs in all_vars.items():
        with (root/f"chr{c}.7snp_eas.dosages.tsv").open("w",newline="") as f:
            w=csv.writer(f,delimiter="\t")
            w.writerow(["EAS_sample_ID"]+vs)
            for i,name in enumerate(samples):
                # Each variable varies in EAS and JPT subset.
                w.writerow([name]+[(i+j)%3 for j in range(len(vs))])
    return all_vars
class LDRealSourceTests(unittest.TestCase):
    def test_chrom4_reversed_REF_ALT_converts_dosage(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);out=root/"out";out.mkdir()
            samples=[f"ID{x:03d}" for x in range(504)]
            file=root/"samples.txt";file.write_text("\n".join(samples)+"\n")
            target=out/"chr4.7snp_eas.grch37.vcf.gz";target.write_bytes(b"valid stub")
            (out/(target.name+".tbi")).write_bytes(b"index stub")
            dosages1=["0|0"]*252+["1|1"]*252
            dosages2=["0|0"]*252+["1|1"]*252
            stdout="\n".join(["\t".join(["4","39413780","G","A"]+dosages1),
                 "\t".join(["4","100239319","T","C"]+dosages2)])+"\n"
            with patch("extract_nonaldh2_alcohol_1kg_eas_ld.subprocess.run",
                 side_effect=[SimpleNamespace(returncode=0,stdout="\n".join(samples)+"\n",stderr=""),
                   SimpleNamespace(returncode=0,stdout=stdout,stderr="")]):
                q=extract_chrom("4",out,file,True)
            self.assertEqual(q["swapped_REF_ALT_source_count"],1)
            self.assertEqual(q["source_genotypes_exact_verified"],2)
            with (out/"chr4.7snp_eas.dosages.tsv").open() as f:
                r=list(csv.DictReader(f,delimiter="\t"))
            self.assertEqual(float(r[0]["4:39413780:A:G"]),2)
            self.assertEqual(float(r[0]["4:100239319:T:C"]),0)
            self.assertEqual(float(r[-1]["4:39413780:A:G"]),0)
    def test_ld7_and_jpt104_without_causal_promotion(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);names=[f"J{i:03d}" for i in range(504)]
            make_sources(root,names)
            outcome=build_ld(root,{c:{"test":c} for c in REGIONS})
            self.assertEqual(outcome["pairwise_r2_results"],21)
            panel=root/"panel.tsv"
            with panel.open("w",newline="") as f:
                w=csv.writer(f,delimiter="\t")
                w.writerow(["sample","pop","super_pop","gender"])
                for i,s in enumerate(names):
                    w.writerow([s,"JPT" if i<104 else "CHB","EAS","female"])
            audit=sensitivity(root,panel)
            self.assertEqual(audit["sample_count_JPT"],104)
            self.assertEqual(audit["pair_count"],21)
            self.assertTrue(audit["JPT_r2_sample_limited_no_causal_MR_gate"])
if __name__=="__main__":unittest.main()
