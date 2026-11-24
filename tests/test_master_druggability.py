import csv,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"run_master_druggability.py"
def write(p,h,rows):
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(h);w.writerows(rows)
def test_druggability_direction(tmp_path):
    drugs=tmp_path/"d.tsv";causal=tmp_path/"c.tsv";out=tmp_path/"out.tsv"
    write(drugs,["gene_symbol","drug_name","mechanism","approved","source"],[
      ["F11","DrugA","inhibitor",1,"manual"],["PROCR","DrugB","agonist",0,"manual"]
    ])
    write(causal,["gene_symbol","beta"],[["F11",0.3],["PROCR",-0.2]])
    p=subprocess.run([sys.executable,str(SCRIPT),"--drug-table",str(drugs),"--causal-direction",str(causal),"--output",str(out)],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    with out.open(encoding="utf-8") as f:rows=list(csv.DictReader(f,delimiter="\t"))
    by={r["gene_symbol"]:r for r in rows}
    assert by["F11"]["translation_tier"]=="A"
    assert by["PROCR"]["translation_tier"]=="B"
