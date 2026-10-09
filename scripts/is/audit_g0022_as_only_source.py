#!/usr/bin/env python3
"""Audit AS-only G0022 pilot markers against original AIS GRCh37 canonical rows.

'AS-only' means missing from the specific AIS reference-matched pilot; distinguish
true absent canonical positions from allele mismatch and upstream QC exclusion.
"""
import csv,gzip,math,json,argparse
from pathlib import Path
from collections import defaultdict,Counter

ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/susie_g0022_pilot")
AIS_CANONICAL=Path("/srv/is-analysis/data/is/processed/gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz")
MATCHED=ROOT.parent/"IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz"
def data(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def source_index(src,lookup):
    hits=defaultdict(list)
    with gzip.open(src,"rt") as f:
        for row in csv.DictReader(f,delimiter="\t"):
            if row["chr"]!="12":continue
            try:p=int(row["pos"])
            except (ValueError,TypeError):continue
            if p in lookup:
                hits[p].append(row)
    return hits
def run(root,ais_path,matched_path):
    prior=set()
    with gzip.open(matched_path,"rt") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            if r["dataset"]=="GCST90104545" and r["group_id"]=="IS_XDATA_G0022" and r["qc_status"] in ("MATCH_ALT_EFFECT","MATCH_REF_EFFECT"):
                prior.add(r["variant_id"])
    chosen=[]
    for path in sorted(root.glob("GCST90104544_AS_signal*/variants.tsv")):
        for r in data(path):
            if r["variant_id"] not in prior:
                chosen.append((path.parent.name,r))
    if not chosen:raise ValueError("No AS-only markers found")
    positions={int(r["pos"]) for _,r in chosen}
    lookup=source_index(ais_path,positions)
    reports=[]
    for name,r in chosen:
        vid=r["variant_id"];pos=int(r["pos"])
        alleles=set(vid.split(":")[2:4])
        at=lookup.get(pos,[])
        matched=[x for x in at if set([x["effect_allele"],x["other_allele"]])==alleles]
        if not at:state="NO_AIS_CANONICAL_POSITION"
        elif not matched:state="POSITION_PRESENT_DIFFERENT_ALLELES"
        else:state="ALLELE_PAIR_PRESENT_BUT_FAILED_JOIN_OR_QC"
        report={"window":name,"variant_id":vid,"pos":pos,
          "AS_z":r["gwas_z"],"AS_p":r["gwas_p"],
          "AS_reference_alt_af":r["reference_alt_af"],
          "AS_gwas_alt_eaf":r["gwas_alt_eaf"],
          "AIS_canonical_rows_at_position":len(at),
          "AIS_canonical_matching_allele_pairs":len(matched),
          "AIS_canonical_first_matching_p":matched[0]["p"] if matched else "",
          "audit_status":state,
          "interpretation":"DO_NOT_AUTO_DISCARD_VARIANT_OR_INFER_TECHNICAL_ERROR"}
        reports.append(report)
    with (root/"G0022_AS_ONLY_SNP_CANONICAL_SOURCE_AUDIT.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(reports[0]),delimiter="\t")
        w.writeheader();w.writerows(reports)
    counts=Counter(x["audit_status"] for x in reports)
    bywin=Counter(x["window"] for x in reports)
    summary={"as_only_marker_rows":len(reports),"unique_variant_ids":len({x["variant_id"] for x in reports}),
      "by_window":dict(bywin),"source_presence_status":dict(counts),
      "max_abs_z_as_only":max(abs(float(x["AS_z"])) for x in reports),
      "as_only_gws_p_le_5e8":sum(float(x["AS_p"])<=5e-8 for x in reports),
      "interpretation":"Presence/QC audit only, not proof of error and not basis for changing canonical GWAS."}
    (root/"G0022_AS_ONLY_SNP_CANONICAL_SOURCE_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    return summary
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,default=ROOT)
    ap.add_argument("--ais",type=Path,default=AIS_CANONICAL)
    ap.add_argument("--matched",type=Path,default=MATCHED)
    a=ap.parse_args()
    run(a.root,a.ais,a.matched)
