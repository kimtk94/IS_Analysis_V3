import csv
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"run_master_transcriptomics.py"

def write(path,header,rows):
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(header);w.writerows(rows)

def test_transcript_evidence_tiers(tmp_path):
    pmr=tmp_path/"pmr.tsv"; smr=tmp_path/"smr.tsv"; coloc=tmp_path/"coloc.tsv"; out=tmp_path/"out.tsv"
    write(pmr,["gene_symbol","phenotype","beta","se","p"],[
      ["G1","CKD",0.3,0.1,0.001],
      ["G2","CKD",-0.2,0.1,0.02],
      ["G3","CKD",0.2,0.1,0.03],
    ])
    write(smr,["Gene","ProbeID","b_SMR","se_SMR","p_SMR","p_HEIDI","nsnp_HEIDI"],[
      ["G1","P1",0.4,0.1,0.001,0.50,20],
      ["G2","P2",-0.3,0.1,0.002,0.20,18],
      ["G3","P3",-0.2,0.1,0.003,0.01,15],
    ])
    write(coloc,["locus","PP.H3","PP.H4","p12"],[
      ["G1",0.05,0.90,1e-5],
      ["G2",0.20,0.50,1e-5],
      ["G3",0.70,0.10,1e-5],
    ])
    p=subprocess.run([
      sys.executable,str(SCRIPT),
      "--protein-mr",str(pmr),"--smr",str(smr),
      "--phenotype","CKD","--eqtl-coloc",str(coloc),
      "--output",str(out)
    ],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    with out.open(encoding="utf-8") as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    by={r["gene_symbol"]:r for r in rows}
    assert by["G1"]["transcript_evidence_tier"]=="T1"
    assert by["G2"]["transcript_evidence_tier"]=="T2"
    assert by["G3"]["transcript_evidence_tier"]=="CONFLICT"
    assert by["G1"]["heidi_status"]=="no_evidence_of_heterogeneity"
