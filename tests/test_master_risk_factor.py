import csv
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"run_master_risk_factor.py"

def write(path,header,rows):
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(header);w.writerows(rows)

def test_risk_factor_mechanism(tmp_path):
    d=tmp_path/"disease.tsv";r=tmp_path/"risk.tsv";s=tmp_path/"sign.tsv";out=tmp_path/"out.tsv"
    h=["gene_symbol","phenotype","beta","se","p"]
    write(d,h,[["FURIN","IS",0.3,0.1,0.001]])
    write(r,h,[["FURIN","SBP",0.4,0.1,0.002],["FURIN","LDL_C",-0.1,0.1,0.4]])
    write(s,["phenotype","disease_risk_sign"],[["SBP",1],["LDL_C",1]])
    p=subprocess.run([sys.executable,str(SCRIPT),"--disease-mr",str(d),"--risk-mr",str(r),"--disease-phenotype","IS","--risk-sign-map",str(s),"--output",str(out)],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    with out.open(encoding="utf-8") as f:rows=list(csv.DictReader(f,delimiter="\t"))
    by={x["risk_factor"]:x for x in rows}
    assert by["SBP"]["mechanism_support"]=="1"
    assert by["LDL_C"]["mechanism_support"]=="0"
