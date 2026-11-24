import csv,json,subprocess,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"run_master_evidence_integration.py"
def write(p,h,rows):
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(h);w.writerows(rows)
def test_final_evidence_tier(tmp_path):
    cand=tmp_path/"cand.tsv";mr=tmp_path/"mr.tsv";col=tmp_path/"col.tsv";anc=tmp_path/"anc.tsv";loc=tmp_path/"loc.tsv";ind=tmp_path/"ind.tsv";man=tmp_path/"m.tsv";out=tmp_path/"out.tsv"
    write(cand,["gene_symbol"],[["G1"],["G2"]])
    write(mr,["gene_symbol","p","fdr_bh"],[["G1",1e-5,2e-5],["G2",0.01,0.02]])
    write(col,["gene_symbol","PP.H4"],[["G1",0.95],["G2",0.90]])
    write(anc,["gene_symbol","replication_tier"],[["G1","A"],["G2","D"]])
    write(loc,["gene_symbol","localization_evidence_class"],[["G1","TISSUE_CELL"],["G2","TISSUE_CELL"]])
    write(ind,["gene_symbol","status","p"],[["G1","ok",0.03],["G2","ok",0.2]])
    write(man,["stage","file","key_column"],[
      ["mr",str(mr),"gene_symbol"],["coloc",str(col),"gene_symbol"],["ancestry",str(anc),"gene_symbol"],
      ["localization",str(loc),"gene_symbol"],["individual",str(ind),"gene_symbol"]
    ])
    p=subprocess.run([sys.executable,str(SCRIPT),"--candidates",str(cand),"--manifest",str(man),"--output",str(out)],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    with out.open(encoding="utf-8") as f:rows=list(csv.DictReader(f,delimiter="\t"))
    by={r["gene_symbol"]:r for r in rows}
    assert by["G1"]["final_tier"]=="Tier 1"
    assert by["G2"]["final_tier"]=="Tier 3"
    assert "ANCESTRY_DIRECTION" in by["G2"]["conflict_flags"]

CASES=json.loads((ROOT/"tests/fixtures/master_evidence_integration/gates.json").read_text())

def check_scientific_evidence_gates(tmp_path,case):
    candidate=tmp_path/"candidate.tsv";manifest=tmp_path/"manifest.tsv";output=tmp_path/"output.tsv"
    write(candidate,["gene_symbol"],[["SYNTHETIC"]])
    stages={stage:[dict(row)] for stage,row in CASES["baseline"].items()}
    for stage,rows in case.get("overrides",{}).items():
        stages[stage]=rows
    entries=[]
    for stage,rows in stages.items():
        if not rows:continue
        path=tmp_path/(stage+".tsv")
        fields=["gene_symbol"]+sorted({key for row in rows for key in row})
        write(path,fields,[["SYNTHETIC"]+[row.get(key,"") for key in fields[1:]] for row in rows])
        entries.append([stage,str(path),"gene_symbol"])
    write(manifest,["stage","file","key_column"],entries)
    completed=subprocess.run([sys.executable,str(SCRIPT),"--candidates",str(candidate),
        "--manifest",str(manifest),"--output",str(output)],cwd=ROOT,capture_output=True,text=True,timeout=10)
    assert completed.returncode==0,completed.stdout+completed.stderr
    with output.open(encoding="utf-8") as stream:
        row=next(csv.DictReader(stream,delimiter="\t"))
    assert row["final_tier"]==case["expected_tier"]
    assert row["conflict_flags"]==case.get("expected_conflict","")
    for key,value in case.get("expected_scores",{}).items():
        assert row["score_"+key]==str(value)

class ScientificEvidenceGateTests(unittest.TestCase):
    def test_fixture_evidence_gates(self):
        for case in CASES["cases"]:
            with self.subTest(case=case["name"]), tempfile.TemporaryDirectory() as tmp:
                check_scientific_evidence_gates(Path(tmp),case)

if __name__=="__main__":
    unittest.main()
