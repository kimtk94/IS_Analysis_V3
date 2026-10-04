import csv
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"run_master_mr_robustness.py"

def test_mr_robustness_smoke(tmp_path):
    inp=tmp_path/"harmonized.tsv"
    out=tmp_path/"robust.tsv"
    fields=[
      "protein_id","gene_symbol","ancestry","phenotype","rsid",
      "beta_exposure","se_exposure","beta_outcome","se_outcome",
      "exposure_n","outcome_n"
    ]
    rows=[
      ["P1","GENE1","EUR","IS","rs1",0.20,0.02,0.10,0.03,50000,100000],
      ["P1","GENE1","EUR","IS","rs2",0.15,0.02,0.08,0.03,50000,100000],
      ["P1","GENE1","EUR","IS","rs3",0.25,0.03,0.12,0.04,50000,100000],
      ["P1","GENE1","EUR","IS","rs4",0.18,0.02,0.09,0.03,50000,100000],
    ]
    with inp.open("w",newline="",encoding="utf-8") as fh:
        w=csv.writer(fh,delimiter="\t",lineterminator="\n")
        w.writerow(fields); w.writerows(rows)
    p=subprocess.run(
      [sys.executable,str(SCRIPT),"--input",str(inp),"--output",str(out)],
      cwd=ROOT,capture_output=True,text=True
    )
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    with out.open(encoding="utf-8") as fh:
        got=list(csv.DictReader(fh,delimiter="\t"))
    assert len(got)==1
    r=got[0]
    assert r["protein_id"]=="P1"
    assert int(r["n_instruments"])==4
    assert float(r["ivw_beta"])>0
    assert 0 <= float(r["cochran_q_p"]) <= 1
    assert r["egger_intercept"] != ""
    assert r["weighted_median_beta"] != ""
    assert int(r["steiger_n"])==4
