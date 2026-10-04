import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
MAP=ROOT/"config"/"ckd"/"master_stage_map.tsv"
AUDIT=ROOT/"scripts"/"audit_ckd_master_existing_outputs.py"
WRAPPER=ROOT/"server"/"ckd_master_end_to_end_audit.sh"

def test_ckd_stage_map_is_complete():
    with MAP.open(encoding="utf-8") as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    assert len(rows)==20
    assert [int(r["stage"]) for r in rows]==list(range(20))
    assert all(r["current_status"] in {"EXISTING","PARTIAL","PLANNED"} for r in rows)

def test_ckd_existing_output_audit_smoke(tmp_path):
    existing=tmp_path/"existing.tsv"
    existing.write_text("gene_symbol\tp\nG1\t0.01\n",encoding="utf-8")
    missing=tmp_path/"missing.tsv"
    m=tmp_path/"map.tsv"
    with m.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(["stage","stage_name","current_status","path_pattern","canonical_use","note"])
        w.writerow([0,"dataset_audit","EXISTING",str(existing),"x",""])
        w.writerow([1,"exposure","PLANNED",str(missing),"x",""])
    out=tmp_path/"readiness.tsv"
    summary=tmp_path/"summary.json"
    p=subprocess.run([
      sys.executable,str(AUDIT),"--map",str(m),
      "--output",str(out),"--summary",str(summary)
    ],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    meta=json.loads(summary.read_text())
    assert meta["artifacts_found"]==1
    assert meta["existing_stage_artifacts_found"]==[0]
    assert meta["planned_stages"]==[1]

def test_ckd_master_wrapper_bash_syntax():
    p=subprocess.run(["bash","-n",str(WRAPPER)],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stderr
