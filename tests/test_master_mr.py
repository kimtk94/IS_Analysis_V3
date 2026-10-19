import csv
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"run_master_mr.py"

def write(path,header,rows):
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(header);w.writerows(rows)

def test_generic_mr_smoke(tmp_path):
    exp=tmp_path/"exp.tsv"; out=tmp_path/"out.tsv"
    harm=tmp_path/"harm.tsv"; mr=tmp_path/"mr.tsv"
    write(exp,
      ["protein_id","gene_symbol","rsid","effect_allele","other_allele","beta","se","eaf","n"],
      [
        ["P1","GENE1","rs1","A","G",0.20,0.02,0.30,50000],
        ["P1","GENE1","rs2","C","T",0.15,0.02,0.25,50000],
        ["P2","GENE2","rs3","G","A",0.25,0.03,0.20,50000],
      ])
    write(out,
      ["rsid","effect_allele","other_allele","beta","se","eaf","n","effect_scale"],
      [
        ["rs1","A","G",0.10,0.03,0.30,100000,"log_odds"],
        ["rs2","T","C",-0.08,0.03,0.75,100000,"log_odds"],
        ["rs3","G","A",-0.05,0.04,0.20,100000,"log_odds"],
      ])
    p=subprocess.run([
      sys.executable,str(SCRIPT),
      "--exposure",str(exp),"--outcome",str(out),
      "--harmonized-output",str(harm),"--mr-output",str(mr),
      "--ancestry","EUR","--phenotype","IS"
    ],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    with mr.open(encoding="utf-8") as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    assert len(rows)==2
    p1=next(r for r in rows if r["protein_id"]=="P1")
    p2=next(r for r in rows if r["protein_id"]=="P2")
    assert p1["method"]=="IVW"
    assert p2["method"]=="Wald ratio"
    assert float(p1["beta"])>0
    assert float(p2["beta"])<0
