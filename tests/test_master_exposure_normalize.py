import csv,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"run_master_exposure_normalize.py"
def write(p,h,rows):
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(h);w.writerows(rows)
def test_exposure_normalize(tmp_path):
    src=tmp_path/"src.tsv";mp=tmp_path/"map.tsv";out=tmp_path/"norm.tsv"
    write(src,["protein","gene","snp","ea","oa","b","se","eaf"],[
      ["P1","g1","rs1","a","g",0.2,0.02,0.3]
    ])
    write(mp,["canonical_column","source_column"],[
      ["protein_id","protein"],["gene_symbol","gene"],["rsid","snp"],
      ["effect_allele","ea"],["other_allele","oa"],["beta","b"],["se","se"],["eaf","eaf"]
    ])
    p=subprocess.run([sys.executable,str(SCRIPT),"--input",str(src),"--column-map",str(mp),"--output",str(out),"--ancestry","EUR","--genome-build","GRCh37","--platform","Olink"],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    with out.open(encoding="utf-8") as f:r=list(csv.DictReader(f,delimiter="\t"))[0]
    assert r["gene_symbol"]=="G1"
    assert r["effect_allele"]=="A"
    assert r["ancestry"]=="EUR"
