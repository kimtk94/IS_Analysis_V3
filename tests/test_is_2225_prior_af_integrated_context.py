"""Broad candidate retention and exact stable-gene prior/AF join guardrails."""
import importlib.util
import sys
import unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"scripts/is/build_is_2225_prior_af_integrated_context.py"
spec=importlib.util.spec_from_file_location("prior_af_gene",SCRIPT)
mod=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=mod
spec.loader.exec_module(mod)

def fixture():
    genes=[];links=[];grid=[];af=[];ns=[];replay=[]
    for i in range(2225):
        gid=f"ENSG{i:011d}"
        n=2 if i<200 else 1
        genes.append({
          "gene_id_stable":gid,"n_positional_regions":str(n),
          "archived_EAS_GTEx_ABF_assays":"646" if i==0 else "0",
          "archived_EAS_GTEx_max_PP_H4":"0.1" if i==0 else "",
        })
        for j in range(n):
            links.append({
              "gene_id_stable":gid,
              "ancestry":"EAS_AND_JAPANESE" if i==0 and j==0 else "EUR",
              "legacy_original_bbj_locus":"BBJ_IS_L001" if i==0 and j==0 else "",
              "archived_ABF_assays":"646" if i==0 and j==0 else "0",
            })
    for i in range(646):
        k={"locus":"BBJ_IS_L001","dataset_key":f"GTEx_V8__Tissue_{i:03d}",
           "gene_base":genes[0]["gene_id_stable"]}
        replay.append({**k,"status":"PASS","delta_max":"0","new_PP_H4":"0.1"})
        af.append({**k,"fraction_MAF_difference_GT_0p1":"0.4","n_palindromic":"2"})
        ns.append({**k,"max_abs_H4_shift_vs_first":"0.001","archived_H4":"0.1"})
        for p in mod.PRIORS:
            grid.append({**k,"p12":str(p),"PP_H4":"0.1","PP_H3":"0.2"})
    return genes,links,grid,af,ns,replay

class TestGenePriorAFOverlay(unittest.TestCase):
    def test_full_positional_universe_and_646_assays_preserved(self):
        merged,details=mod.integrate(*fixture())
        self.assertEqual(len(merged),2225)
        self.assertEqual(len(details),646)
        self.assertEqual(sum(x["new_DIRECT_prior_grid_assays"] for x in merged),646)
        self.assertEqual(sum(x["science_gate"]=="MOLECULAR_QTL_NOT_TESTED" for x in merged),2224)
        self.assertTrue(all(x["causal_status"]=="NOT_ESTABLISHED" for x in merged))

    def test_eur_only_locus_must_not_inherit_BBJ_qtl(self):
        g,l,grid,af,ns,replay=fixture()
        l[0]["ancestry"]="EUR"
        with self.assertRaises(ValueError):
            mod.integrate(g,l,grid,af,ns,replay)

    def test_missing_one_prior_fails(self):
        g,l,grid,af,ns,replay=fixture()
        with self.assertRaises(ValueError):
            mod.integrate(g,l,grid[:-1],af,ns,replay)

    def test_discordant_reproduced_h4_fails(self):
        g,l,grid,af,ns,replay=fixture()
        replay[0]["new_PP_H4"]="0.9"
        with self.assertRaises(ValueError):
            mod.integrate(g,l,grid,af,ns,replay)

if __name__=="__main__":
    unittest.main()
