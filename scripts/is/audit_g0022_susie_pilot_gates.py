#!/usr/bin/env python3
"""Integrate LD-consistency and SuSiE pilot sensitivity; never promote to causal.

Requires 5 real windows, two hypothetical sample-size assumptions each.
"""
import argparse,csv,json
from collections import defaultdict
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/susie_g0022_pilot")
def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def build(root):
    results=read(root/"G0022_SUSIE_RSS_PILOT_SUMMARY.tsv")
    diag=read(root/"G0022_GWAS_LD_CONSISTENCY.tsv")
    inputs=json.loads((root/"G0022_SUSIE_PILOT_INPUT_MANIFEST.json").read_text())
    if len(inputs)!=5 or len(results)!=10 or len(diag)!=10:
        raise ValueError("Expected five windows and 10 two-N scenarios")
    if any(x["status"]!="PILOT_INPUT_QC_PASS" for x in inputs):
        raise ValueError("Not all input regions passed initial QC")
    by_r=defaultdict(dict);by_d=defaultdict(dict)
    for x in results:
        key=x["signal"];n=int(x["hypothetical_n"])
        if n in by_r[key]:raise ValueError("Duplicate SuSiE scenario")
        by_r[key][n]=x
    for x in diag:
        key=x["signal"];n=int(x["hypothetical_n"])
        if n in by_d[key]:raise ValueError("Duplicate LD diagnostic")
        by_d[key][n]=x
    out=[]
    for i in inputs:
        signal=i["signal"]
        runs=by_r.get(signal,{})
        checks=by_d.get(signal,{})
        if set(runs)!={20000,256274} or set(checks)!={20000,256274}:
            raise ValueError("Missing pilot N assumption in "+signal)
        max_s=max(float(d["ld_mismatch_s"]) for d in checks.values())
        converged=all(x["converged"]=="TRUE" for x in runs.values())
        cs_stable=len({x["cs_count"] for x in runs.values()})==1
        top_stable=len({x["top_pip_variant"] for x in runs.values()})==1
        if max_s>0.2:
            tag="BLOCKED_VERY_HIGH_LD_MISMATCH"
        elif max_s>0.05:
            tag="BLOCKED_LD_MISMATCH"
        elif not converged:
            tag="BLOCKED_NONCONVERGENCE"
        else:
            tag="EXPLORATORY_SIGNAL_ONLY_NEEDS_EFFECTIVE_N_AND_FULL_LOCI"
        out.append({"signal":signal,"study":i["study"],"trait":i["trait"],
            "n_variants":i["variant_count"],"ref_samples":i["reference_samples"],
            "n_scenarios":2,"both_converged":int(converged),
            "max_ld_mismatch_s":round(max_s,6),
            "credible_set_counts":";".join(runs[n]["cs_count"] for n in sorted(runs)),
            "top_variant_stable_between_n":int(top_stable),
            "cs_count_stable_between_n":int(cs_stable),
            "scientific_gate":tag,
            "claim":"NOT_FINE_MAPPING_VALIDATED"})
    with (root/"G0022_SUSIE_PILOT_EVIDENCE_GATES.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter="\t")
        w.writeheader();w.writerows(out)
    report={"windows":len(out),"total_scenarios":len(results),
        "both_n_converged":sum(x["both_converged"] for x in out),
        "high_ld_mismatch":sum(x["scientific_gate"].startswith("BLOCKED_") and "LD_MISMATCH" in x["scientific_gate"] for x in out),
        "exploratory_low_mismatch_converged":sum(x["scientific_gate"].startswith("EXPLORATORY_") for x in out),
        "final_fine_mapping_validated":0,
        "status":"PILOT_EVIDENCE_QC_ONLY",
        "warning":"No case/control effective N; window truncation; ancestry and GWAS LD mismatch; GIGASTROKE AS/AIS overlapping samples."}
    (root/"G0022_SUSIE_PILOT_EVIDENCE_GATES_SUMMARY.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    return report
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,default=ROOT)
    build(ap.parse_args().root)
