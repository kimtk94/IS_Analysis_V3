import csv
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from importlib import import_module
build = import_module("scripts.is.audit_is_calhm2_pairwise_ld").build

class TestCALHM2PairwiseLD(unittest.TestCase):
    def make_inputs(self,root):
        qtl=root/"qtl.tsv"
        geno=root/"geno.raw"
        panel=root/"population.tsv"
        keys=["10:105410001:A:G","10:105410002:C:T","10:105459834:T:C"]
        with qtl.open("w",newline="") as f:
            fields=["locus","symbol","SNP_GRCh37",
                "BBJ_exact_Fasta37","GIGA_EAS_AIS_exact_Fasta37",
                "palindromic_review","QTD000161_nominal_p"]
            writer=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
            writer.writeheader()
            for k in keys:
                writer.writerow(dict(locus="BBJ_IS_L002",symbol="CALHM2",
                    SNP_GRCh37=k,BBJ_exact_Fasta37="True",
                    GIGA_EAS_AIS_exact_Fasta37="True",
                    palindromic_review="False",QTD000161_nominal_p="0.00001"))
        ids=[f"sample{i}" for i in range(504)]
        with panel.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=["sample","pop"],delimiter="\t")
            w.writeheader()
            for i,k in enumerate(ids):
                w.writerow({"sample":k,"pop":"JPT" if i<104 else "CHB"})
        rng=np.random.default_rng(42)
        a=rng.binomial(2,.4,504)
        b=rng.binomial(2,.5,504)
        c=rng.binomial(2,.3,504)
        frame=pd.DataFrame({"FID":0,"IID":ids,"PAT":0,"MAT":0,"SEX":0,"PHENOTYPE":-9,
            "10:105459834:T:C_C":a,
            "10:105599770:T:G_G":b,
            keys[0]+"_G":a,
            keys[1]+"_T":a,
            "UNUSED_SYNTHETIC_C_C":c})
        frame.to_csv(geno,sep="\t",index=False)
        return qtl,geno,panel
    def test_outputs_and_reference_bounds(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);q,g,p=self.make_inputs(root)
            result=build(q,g,p,root/"out")
            self.assertEqual(len(result),18) # 4 unique variants choose 2 x 3 populations
            pairs=[r for r in result if r["population"]=="JPT104"]
            self.assertEqual(len(pairs),6)
            self.assertTrue(all(r["n"]==104 for r in pairs))
            self.assertTrue(all(0<=r["r2"]<=1.000001 for r in result))
            self.assertTrue((root/"out"/"CALHM2_PAIRWISE_LD_MANIFEST.json").exists())
            with self.assertRaises(FileExistsError):
                build(q,g,p,root/"out")
    def test_snp_count_drift_fails_closed(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);q,g,p=self.make_inputs(root)
            records=pd.read_csv(q,sep="\t").iloc[:2]
            records.to_csv(q,sep="\t",index=False)
            with self.assertRaises(ValueError):
                build(q,g,p,root/"out")
