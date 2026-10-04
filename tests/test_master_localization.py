import csv,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"run_master_localization.py"

def write(p,h,rows):
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(h);w.writerows(rows)

def test_localization_integration(tmp_path):
    cand=tmp_path/"cand.tsv";t=tmp_path/"t.tsv";c=tmp_path/"c.tsv";s=tmp_path/"s.tsv"
    tm=tmp_path/"tm.tsv";cm=tmp_path/"cm.tsv";sm=tmp_path/"sm.tsv";out=tmp_path/"out.tsv"
    write(cand,["gene_symbol"],[["UMOD"],["CPVL"]])
    write(t,["gene_symbol","tissue","expression","tau_specificity"],[["UMOD","kidney",100,0.95],["CPVL","kidney",20,0.7]])
    write(c,["gene_symbol","cell_type","expression"],[["UMOD","loop of Henle",2500],["CPVL","cDC",1500]])
    write(s,["gene_symbol","spatial_region","expression","enrichment"],[["UMOD","TAL",200,3.0]])
    write(tm,["source_file","source_name","tissue_or_region"],[[str(t),"HPA","kidney"]])
    write(cm,["source_file","source_name","tissue_or_region"],[[str(c),"HPA_cell","kidney"]])
    write(sm,["source_file","source_name","tissue_or_region"],[[str(s),"spatial","kidney"]])
    p=subprocess.run([sys.executable,str(SCRIPT),"--candidates",str(cand),"--tissue-manifest",str(tm),"--cell-manifest",str(cm),"--spatial-manifest",str(sm),"--output",str(out)],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    with out.open(encoding="utf-8") as f:rows=list(csv.DictReader(f,delimiter="\t"))
    by={r["gene_symbol"]:r for r in rows}
    assert by["UMOD"]["localization_evidence_class"]=="TISSUE_CELL_SPATIAL"
    assert by["CPVL"]["localization_evidence_class"]=="TISSUE_CELL"
    assert by["UMOD"]["top_cell_type"]=="loop of Henle"
