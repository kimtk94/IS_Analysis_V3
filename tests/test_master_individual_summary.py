import csv
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"summarize_master_individual_validation.py"

def write(path,header,rows):
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(header);w.writerows(rows)

def test_individual_validation_summary(tmp_path):
    r1=tmp_path/"G1.tsv"; man=tmp_path/"manifest.tsv"; out=tmp_path/"summary.tsv"
    write(r1,
      ["model","endpoint","method","status","n","events","beta","se","statistic","p","predictor","effect_scale","effect_ratio","ci95_low","ci95_high","ratio_ci95_low","ratio_ci95_high"],
      [
        ["baseline_egfr","baseline_f0_egfr","linear_regression","ok",100,"",0.2,0.08,2.5,0.012,"dosage","linear_beta","",0.043,0.357,"",""],
        ["incident_ckd_logistic","incident_ckd","logistic_regression","ok",95,10,0.4,0.2,2.0,0.045,"dosage","log_odds_or_log_hazard",1.49,0.008,0.792,1.01,2.21],
      ])
    write(man,["gene_symbol","predictor","result_file"],[["G1","dosage",str(r1)]])
    p=subprocess.run([sys.executable,str(SCRIPT),"--manifest",str(man),"--output",str(out)],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    with out.open(encoding="utf-8") as f: rows=list(csv.DictReader(f,delimiter="\t"))
    assert len(rows)==2
    assert rows[0]["gene_symbol"]=="G1"
    assert rows[0]["nominal_p05"]=="1"
    assert rows[1]["method"]=="logistic_regression"
