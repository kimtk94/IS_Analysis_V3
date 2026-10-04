import csv
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"bridge_ckd_master_localization.py"
WRAPPER=ROOT/"server"/"ckd_master_materialize_existing.sh"

def write(path,header,rows):
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(header);w.writerows(rows)

def test_localization_bridge(tmp_path):
    s3a=tmp_path/"s3a.tsv";s3b=tmp_path/"s3b.tsv";out=tmp_path/"out.tsv"
    write(s3a,["gene_symbol","hirohama_maintext_kidney_pqtl","kidney_eqtl_meta686_hits"],[
      ["UMOD",0,2],["GSTA3",0,0]
    ])
    write(s3b,["gene_symbol","top_cell_type","top_nCPM","top_compartment","top_second_ratio","tau_specificity","positive_celltype_n"],[
      ["UMOD","loop of Henle",2450,"renal_epithelial",20,0.95,3],
      ["GSTA3","loop of Henle",0.2,"renal_epithelial",2,0.5,2],
    ])
    p=subprocess.run([
      sys.executable,str(SCRIPT),"--stage3a",str(s3a),"--stage3b",str(s3b),"--output",str(out)
    ],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    with out.open(encoding="utf-8") as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    by={r["gene_symbol"]:r for r in rows}
    assert by["UMOD"]["localization_evidence_class"]=="TISSUE_CELL"
    assert by["GSTA3"]["localization_evidence_class"]=="SINGLE_LAYER"
    assert by["UMOD"]["spatial_source_n"]=="0"

def test_materialize_wrapper_bash_syntax():
    p=subprocess.run(["bash","-n",str(WRAPPER)],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stderr
