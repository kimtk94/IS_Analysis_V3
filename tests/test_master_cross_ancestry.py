import csv
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"run_master_cross_ancestry.py"

def write(path,header,rows):
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(header);w.writerows(rows)

def test_cross_ancestry_tiers(tmp_path):
    d=tmp_path/"disc.tsv"; r=tmp_path/"repl.tsv"; c=tmp_path/"coloc.tsv"; out=tmp_path/"out.tsv"
    h=["protein_id","gene_symbol","phenotype","beta","se","p","fdr_bh","n_instruments"]
    write(d,h,[
      ["P1","G1","IS",0.5,0.1,1e-6,1e-5,3],
      ["P2","G2","IS",-0.4,0.1,1e-5,1e-4,2],
      ["P3","G3","IS",0.3,0.1,1e-3,0.01,2],
    ])
    write(r,h,[
      ["P1","G1","IS",0.4,0.15,0.01,0.02,2],
      ["P2","G2","IS",-0.2,0.12,0.04,0.05,2],
      ["P3","G3","IS",-0.1,0.12,0.40,0.5,2],
    ])
    write(c,["locus","PP.H3","PP.H4","p12"],[
      ["G1",0.05,0.90,1e-5],
      ["G2",0.10,0.60,1e-5],
      ["G3",0.20,0.20,1e-5],
    ])
    p=subprocess.run([
      sys.executable,str(SCRIPT),
      "--discovery-mr",str(d),"--replication-mr",str(r),
      "--replication-coloc",str(c),
      "--output",str(out),
      "--discovery-ancestry","EUR","--replication-ancestry","EAS"
    ],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    with out.open(encoding="utf-8") as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    by={x["gene_symbol"]:x for x in rows}
    assert by["G1"]["replication_tier"]=="A"
    assert by["G2"]["replication_tier"]=="B"
    assert by["G3"]["replication_tier"]=="D"
    assert by["G1"]["direction"]=="concordant"
