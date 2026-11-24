import csv
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"run_master_phenotype_matrix.py"

def write(path,header,rows):
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(header);w.writerows(rows)

def test_phenotype_matrix(tmp_path):
    a=tmp_path/"a.tsv";b=tmp_path/"b.tsv";m=tmp_path/"manifest.tsv"
    long=tmp_path/"long.tsv";wide=tmp_path/"wide.tsv"
    h=["gene_symbol","phenotype","beta","se","p"]
    write(a,h,[["G1","IS",0.2,0.05,0.001],["G2","IS",-0.3,0.08,0.01]])
    write(b,h,[["G1","CES",0.25,0.07,0.002],["G2","CES",-0.2,0.09,0.03]])
    write(m,["mr_file","label","group","coloc_file"],[
      [str(a),"IS","stroke",""],[str(b),"CES","stroke_subtype",""]
    ])
    p=subprocess.run([sys.executable,str(SCRIPT),"--manifest",str(m),"--long-output",str(long),"--wide-output",str(wide)],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    with wide.open(encoding="utf-8") as f:rows=list(csv.DictReader(f,delimiter="\t"))
    by={r["gene_symbol"]:r for r in rows}
    assert by["G1"]["direction_consistent"]=="1"
    assert by["G2"]["direction_consistent"]=="1"
    assert by["G1"]["IS__evidence"]=="MR"
