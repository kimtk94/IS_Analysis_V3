#!/usr/bin/env python3
"""Standardize common GIGASTROKE column aliases to a canonical TSV."""
import argparse, csv, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workflow.common import existing_file
ALIASES={"chr":["chr","chromosome"],"pos":["pos","position","bp"],"effect_allele":["effect_allele","ea","alt"],"other_allele":["other_allele","oa","ref"],"beta":["beta","effect"],"se":["se","stderr"],"p_value":["p_value","p","pval"]}
def main():
    p=argparse.ArgumentParser(); p.add_argument("--input",type=existing_file,required=True); p.add_argument("--output",type=Path,required=True); a=p.parse_args()
    with a.input.open(newline="") as f: rows=list(csv.DictReader(f,delimiter="\t"))
    if not rows: raise SystemExit("input contains no rows")
    lookup={k.lower():k for k in rows[0]}; selected={canon:next((lookup[x] for x in aliases if x in lookup),None) for canon,aliases in ALIASES.items()}
    missing=[k for k,v in selected.items() if v is None]
    if missing: raise SystemExit("missing columns: "+", ".join(missing))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(ALIASES),delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows({k:r[v] for k,v in selected.items()} for r in rows)
if __name__ == "__main__": main()
