"""Full 4.05-Mb G0022 pipeline: planning and provenance safety checks."""
import argparse,csv,gzip,json,sys,tempfile,unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from prepare_g0022_full_locus_signed_ld import build as prepare,LEADS,source
from audit_g0022_full_locus_gene_qtl import annotate
from audit_g0022_full_locus_publication_gate import audit

def tsv(path,records):
    path.parent.mkdir(exist_ok=True,parents=True)
    with path.open("w") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t")
        w.writeheader();w.writerows(records)
class FullG0022Test(unittest.TestCase):
    def test_safe_plan_no_genotype_export(self):
        with tempfile.TemporaryDirectory() as t:
            base=Path(t);root=base/"root";root.mkdir();out=base/"must_not_create"
            variants=[]
            for i in range(3000):
                variants.append({"group_id":"IS_XDATA_G0022","dataset":"GCST90104545",
                    "qc_status":"MATCH_ALT_EFFECT","variant_id":f"12:{110000000+i}:A:G",
                    "ref":"A","alt":"G","alt_effect_beta":".1","se":".03"})
            for vid in LEADS:
                variants.append({"group_id":"IS_XDATA_G0022","dataset":"GCST90104545",
                    "qc_status":"MATCH_ALT_EFFECT","variant_id":vid,
                    "ref":vid.split(":")[2],"alt":vid.split(":")[3],
                    "alt_effect_beta":".2","se":".03"})
            with gzip.open(root/"IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz","wt") as f:
                w=csv.DictWriter(f,fieldnames=list(variants[0]),delimiter="\t")
                w.writeheader();w.writerows(variants)
            args=argparse.Namespace(root=root,out=out,reference=base/"notreal",
                                    execute=False,max_variants=6500)
            result=prepare(args)
            self.assertEqual(result["source_harmonized_variants"],3003)
            self.assertEqual(result["status"],"PLAN_ONLY")
            self.assertFalse(out.exists())
            self.assertEqual(len(source(root/"IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz")),3003)
    def test_prior_qtl_provenance_and_publication_gate(self):
        with tempfile.TemporaryDirectory() as t:
            base=Path(t);root=base/"full";root.mkdir();pilot=base/"pilot";pilot.mkdir()
            (root/"G0022_FULL_AIS_SUSIE_SUMMARY.json").write_text(json.dumps({
                "status":"FULL_LOCUS_EXPLORATORY_ONLY","converged":True,
                "n_credible_sets":1,"credible_set_variant_ids":{"L1":["12:120:A:G"]},
                "credible_set_sizes":1,"max_pip_snp":"12:120:A:G",
                "max_pip":.4,"n_snps":4000,"approximate_n_eff":70474,
                "n_eff_is_per_snp":False,"final_causal_signal_validated":False,
                "final_causal_gene_validated":False}))
            tsv(root/"G0022_FULL_AIS_PIP.tsv",[dict(variant_id="12:120:A:G",pos=120,pip=.4)])
            genes=[dict(group_id="IS_XDATA_G0022",gene_id=f"ENSG{i}.1",gene_symbol=s,
                   gene_start=str(i*100),gene_end=str(i*100+150),biotype="protein_coding")
                   for i,s in enumerate(("ALDH2","PTPN11","IFT81","ATP2A2"),1)]
            gf=base/"genes.tsv";tsv(gf,genes)
            coloc=[dict(gene_base="ENSG1",locus="BBJ_IS_L003",**{
                "PP.H4":"0.1","tissue_label":"Artery - Aorta"})]
            cf=base/"coloc.tsv";tsv(cf,coloc)
            r=annotate(root,gf,cf)
            self.assertEqual(r["new_AIS_QTL_coloc_completed"],0)
            self.assertEqual(r["legacy_BBJ_coloc_h4_above_0_8"],0)
            (root/"FULL_LOCUS_INPUT_QC.json").write_text(json.dumps({
                "status":"FULL_G0022_SIGNED_LD_INPUT_QC_PASS",
                "genotype_snp_count":4000,"genotype_n":504,
                "source_harmonized_variants":4500}))
            sens=[dict(config=c,converged="TRUE",purity_filtered_cs="1",
               cs_sizes="1",max_pip_variant="12:120:A:G",max_pip="0.4",
               n_assumed="70474") for c in (
                "LOWER_N20000_L8","APPROX_NEFF_L3","APPROX_NEFF_L8_LD_SHRINK1PCT")]
            tsv(root/"G0022_FULL_AIS_SENSITIVITY.tsv",sens)
            tsv(pilot/"G0022_SUSIE_PILOT_EVIDENCE_GATES.tsv",
                [dict(trait="AIS") for i in range(3)]+[dict(trait="AS") for i in range(2)])
            gate=audit(root,pilot)
            self.assertTrue(gate["all_sensitivity_converged_and_same_CS_top_snp"])
            self.assertFalse(gate["causal_signal_publication_ready"])
            self.assertFalse(gate["causal_gene_publication_ready"])
            sens[1]["converged"]="FALSE"
            tsv(root/"G0022_FULL_AIS_SENSITIVITY.tsv",sens)
            gate2=audit(root,pilot)
            self.assertFalse(gate2["all_sensitivity_converged_and_same_CS_top_snp"])
            self.assertFalse(gate2["causal_signal_publication_ready"])
if __name__=="__main__":
    unittest.main()
