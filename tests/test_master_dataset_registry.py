import csv,subprocess,sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"run_master_dataset_registry.py"
def write(p,h,rows):
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(h);w.writerows(rows)
def test_dataset_registry(tmp_path):
    data=tmp_path/"data.tsv";man=tmp_path/"manifest.tsv";out=tmp_path/"registry.tsv";summary=tmp_path/"summary.json"
    data.write_text("gene\tbeta\nG1\t0.2\n",encoding="utf-8")
    write(man,["dataset_id","role","path","ancestry","genome_build","phenotype","platform","source_name"],[
      ["d1","exposure",str(data),"EUR","GRCh37","","Olink","test"]
    ])
    p=subprocess.run([sys.executable,str(SCRIPT),"--manifest",str(man),"--output",str(out),"--summary",str(summary)],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    meta=json.loads(summary.read_text())
    assert meta["datasets"]==1
    assert meta["missing_files"]==[]
