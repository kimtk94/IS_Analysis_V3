"""Contract: molecular evidence is anchored to a GWAS locus and Ensembl gene."""
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"scripts/is/build_locus_gene_evidence_v2.py"

def make(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        w=csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

class TestLocusAwareEvidence(unittest.TestCase):
    def test_no_symbol_only_or_cross_locus_transfer(self):
        with tempfile.TemporaryDirectory() as t:
            t=Path(t)
            root=t/"out"
            root.mkdir()
            base=t/"base"
            make(root/"IS_ALL_GENE_WINDOW_UNIVERSE.tsv",[
                dict(group_id="G1",gene_id="ENSG111.4",gene_symbol="FGF5",biotype="protein_coding",gene_start=100,gene_end=200,lead_distance_bp=0),
                dict(group_id="G2",gene_id="ENSG222.1",gene_symbol="FGF5",biotype="protein_coding",gene_start=200,gene_end=300,lead_distance_bp=40)])
            make(root/"CROSS_DATASET_INTERVAL_GROUPS.tsv",[
                dict(group_id="G1",chr="4",start=100,end=200,has_gws=1,lead_p="1e-10",phenotypes="AIS"),
                dict(group_id="G2",chr="5",start=200,end=300,has_gws=0,lead_p="2e-7",phenotypes="CES")])
            make(root/"LEGACY_ANCHOR_CANDIDATES.tsv",[dict(gene="FGF5"),dict(gene="NEURL1")])
            make(base/"stage1_gwas_qc/japan/bbj/BBJ_IS_PROVISIONAL_LOCI.tsv",[
                dict(locus_id="BBJ_IS_L001",chr="4",start=120,end=140),
                dict(locus_id="BBJ_IS_L002",chr="6",start=100,end=200)])
            make(base/"stage5_functional/phase9c_convergence/COLOC_ABF_MASTER_ANNOTATED_V2.tsv",[
                dict(locus="BBJ_IS_L001",gene_base="ENSG111",gene_symbol="FGF5",status="PASS",**{"PP.H4":"0.91","PP.H3":"0.08"},tissue_label="Brain",dataset_key="GTEx"),
                dict(locus="BBJ_IS_L002",gene_base="ENSG222",gene_symbol="FGF5",status="PASS",**{"PP.H4":"0.98","PP.H3":"0.01"},tissue_label="Artery",dataset_key="GTEx")])
            make(base/"stage5_functional/phase9f_e_human_r3_4_import/HUMAN_TARGET_GENE_SUMMARY.tsv",[
                dict(gene="FGF5",present="TRUE",top_celltype="EC"),
                dict(gene="NEURL1",present="TRUE",top_celltype="Astrocyte")])
            p=subprocess.run([sys.executable,str(SCRIPT),"--root",str(root),"--base",str(base)],capture_output=True,text=True)
            self.assertEqual(p.returncode,0,p.stderr)
            report=json.loads((root/"IS_LOCUS_GENE_EVIDENCE_V2_SUMMARY.json").read_text())
            self.assertEqual(report["locus_gene_matched_coloc"],1)
            with (root/"IS_LOCUS_GENE_EVIDENCE_V2.tsv").open() as f:
                rows=list(csv.DictReader(f,delimiter="\t"))
            r1=next(x for x in rows if x["group_id"]=="G1")
            r2=next(x for x in rows if x["group_id"]=="G2")
            rn=next(x for x in rows if x["gene_symbol"]=="NEURL1")
            self.assertEqual(r1["coloc_max_pp_h4"],"0.91")
            self.assertEqual(r2["coloc_status"],"NOT_TESTED_IN_THIS_GROUP")
            self.assertEqual(r2["coloc_max_pp_h4"],"")
            self.assertEqual(rn["molecular_scope"],"ANCHOR_OUTSIDE_CURRENT_WINDOWS")
            self.assertEqual(rn["human_atlas_status"],"PRESENT")
if __name__=="__main__":
    unittest.main()
