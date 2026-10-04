import csv,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"run_master_phewas.py"
def write(p,h,rows):
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(h);w.writerows(rows)
def test_phewas(tmp_path):
    src=tmp_path/"p.tsv";man=tmp_path/"m.tsv";long=tmp_path/"long.tsv";summary=tmp_path/"sum.tsv"
    write(src,["gene_symbol","phenotype","beta","se","p","direction_is_adverse"],[
      ["F11","VTE",0.4,0.1,1e-8,1],["F11","LDL",-0.2,0.1,0.5,0]
    ])
    write(man,["source_file","source_name","phenotype_group"],[[str(src),"FinnGen","disease"]])
    p=subprocess.run([sys.executable,str(SCRIPT),"--manifest",str(man),"--long-output",str(long),"--summary-output",str(summary)],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    with summary.open(encoding="utf-8") as f:r=list(csv.DictReader(f,delimiter="\t"))[0]
    assert r["gene_symbol"]=="F11"
    assert r["safety_attention"]=="1"
