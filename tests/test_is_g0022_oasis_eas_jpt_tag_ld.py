"""G0022 OASIS tag-LD fail-closed offline regression."""
import csv
import numpy as np
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_g0022_oasis_eas_jpt_tag_ld import (
    ANCHOR,PROXY,pairwise_ld,panel_mapping,read_ld_genotypes,
    read_qtl_pairs,analyze_group
)

class TestOasisTagLd(unittest.TestCase):
    def test_pairwise_high_tag(self):
        x=np.array([0,0,1,1,2,2]*20,dtype=float)
        y=x.copy()
        r=pairwise_ld(x,y,min_n=80)
        self.assertEqual(r["n_pair"],120)
        self.assertAlmostEqual(r["r2"],1,places=10)
        self.assertGreater(r["r"],.99)

    def test_null_genotype_cannot_be_tag(self):
        x=np.array([0,1,2]*40,dtype=float)
        y=x.copy();y[:60]=float("nan")
        r=pairwise_ld(x,y,min_n=80)
        self.assertIsNone(r["r2"])
        self.assertEqual(r["n_pair"],60)

    def test_jpt_must_be_exact_subgroup(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"panel.tsv"
            with p.open("w") as f:
                f.write("sample\tpop\tsuper_pop\tgender\n")
                for i in range(504):
                    f.write(f"S{i:03d}\t{'JPT' if i<104 else 'CHB'}\tEAS\tmale\n")
            with self.assertRaisesRegex(ValueError,"composition"):
                panel_mapping(p,[f"S{i:03d}_S{i:03d}" for i in range(504)])

    def test_source_csv_block_if_causal_promotion(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"regional.tsv"
            row={"gene":"ALDH2","cell":"Mono-L1",
                 "gwas_variant_grch37":PROXY,
                 "oasis_rsid":"rs11066015",
                 "oasis_p":"0.000913","gwas_p":"9.515e-18",
                 "source_ref_alt_id_match":"EXACT_AFTER_PRIMARY_UCSC_CHAIN",
                 "valid_for_full_locus_coloc":"True"}
            with p.open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=list(row),delimiter="\t")
                w.writeheader();w.writerow(row)
            with self.assertRaisesRegex(ValueError,"closed coloc gate"):
                read_qtl_pairs(p)

    def test_ld_lead_descriptive(self):
        a=np.array([0,1,2]*50+ [0,1,2]*18,dtype=float)
        j=np.arange(104)
        v=a.copy()
        gen={ANCHOR:a,PROXY:v}
        sample={"gene":"ALDH2","cell":"Mono-L1","gwas_variant_grch37":PROXY,
           "oasis_rsid":"rs11066015","oasis_p":"0.000913","gwas_p":"9.515e-18"}
        report=analyze_group("ALDH2","Mono-L1",[sample],gen,j)
        self.assertEqual(report["n_EAS_rs671_r2_ge_08"],1)
        self.assertEqual(report["n_EAS_rs671_r2_ge_08_QTL_nominal_p_lt_05"],1)
        self.assertEqual(report["claim"],"PANEL_TAGGING_DESCRIPTIVE_ONLY")
        self.assertAlmostEqual(report["strongest_oasis_QTL_ref_LD"]["JPT_r2"],1.0)

if __name__=="__main__":
    unittest.main()
