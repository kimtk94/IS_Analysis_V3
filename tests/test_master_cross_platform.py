import csv
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"run_master_cross_platform.py"

def write(path,header,rows):
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(header);w.writerows(rows)

def test_cross_platform_with_mapping(tmp_path):
    d=tmp_path/"olink.tsv"; r=tmp_path/"soma.tsv"
    m=tmp_path/"map.tsv"; c=tmp_path/"coloc.tsv"; out=tmp_path/"out.tsv"
    h=["assay_id","gene_symbol","phenotype","beta","se","p","fdr_bh","n_instruments"]
    write(d,h,[
      ["OID1","F11","IS",0.5,0.1,1e-6,1e-5,2],
      ["OID2","MMUT","IS",-0.4,0.1,1e-5,1e-4,2],
      ["OID3","CEP85","IS",0.3,0.1,1e-3,0.01,2],
    ])
    write(r,h,[
      ["SID11","F11","IS",0.4,0.15,0.01,0.02,2],
      ["SID22","MMUT","IS",-0.2,0.12,0.04,0.05,2],
      ["SID33","CEP85","IS",-0.1,0.12,0.40,0.5,2],
    ])
    write(m,["discovery_id","replication_id","gene_symbol"],[
      ["OID1","SID11","F11"],
      ["OID2","SID22","MMUT"],
      ["OID3","SID33","CEP85"],
    ])
    write(c,["locus","PP.H3","PP.H4","p12"],[
      ["F11",0.05,0.90,1e-5],
      ["MMUT",0.10,0.60,1e-5],
      ["CEP85",0.20,0.20,1e-5],
    ])
    p=subprocess.run([
      sys.executable,str(SCRIPT),
      "--discovery-mr",str(d),"--replication-mr",str(r),
      "--mapping",str(m),
      "--discovery-id-column","assay_id",
      "--replication-id-column","assay_id",
      "--replication-coloc",str(c),
      "--discovery-platform","Olink",
      "--replication-platform","SomaScan",
      "--output",str(out)
    ],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    with out.open(encoding="utf-8") as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    by={x["gene_symbol"]:x for x in rows}
    assert by["F11"]["platform_replication_tier"]=="A"
    assert by["MMUT"]["platform_replication_tier"]=="B"
    assert by["CEP85"]["platform_replication_tier"]=="D"
    assert by["F11"]["direction"]=="concordant"
