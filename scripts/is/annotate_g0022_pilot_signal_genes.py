#!/usr/bin/env python3
"""Annotate G0022 pilot top-PIP variants by nearby gene bodies, no causal assignment."""
import argparse,csv,json
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/susie_g0022_pilot")
BASE=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
def read(path):
    with path.open() as h:return list(csv.DictReader(h,delimiter="\t"))
def build(root,base):
    gate=read(root/"G0022_SUSIE_PILOT_EVIDENCE_GATES.tsv")
    runs=read(root/"G0022_SUSIE_RSS_PILOT_SUMMARY.tsv")
    byrun={r["signal"]:r for r in runs if r["hypothetical_n"]=="256274"}
    gene=read(base/"IS_ALL_GENE_WINDOW_UNIVERSE.tsv")
    candidates=[r for r in gene if r["group_id"]=="IS_XDATA_G0022"]
    if not candidates:raise ValueError("Missing G0022 gene universe")
    out=[]
    for line in gate:
        signal=line["signal"]
        fit=byrun.get(signal)
        if not fit:raise ValueError("Missing hypothetical 256274-sample pilot fit")
        variant=fit["top_pip_variant"]
        tokens=variant.split(":")
        if len(tokens)!=4 or tokens[0]!="12":
            raise ValueError("Invalid G0022 top SNP identifier")
        pos=int(tokens[1])
        local=[]
        for g in candidates:
            low=int(g["gene_start"]);high=int(g["gene_end"])
            distance=max(low-pos,pos-high,0)
            if distance>500000:continue
            local.append((distance,g))
        local.sort(key=lambda x:(x[0],x[1]["gene_id"]))
        for rank,(distance,g) in enumerate(local[:5],1):
            out.append({
                "signal":signal,"lead_top_pip_snp":variant,
                "pilot_top_pip":fit["top_pip"],
                "gene_proximity_rank":rank,
                "gene_id":g["gene_id"],"gene_symbol":g["gene_symbol"],
                "biotype":g["biotype"],"gene_start":g["gene_start"],
                "gene_end":g["gene_end"],"gene_body_distance_bp":distance,
                "overlaps_gene_body":int(distance==0),
                "evidence_gate":line["scientific_gate"],
                "annotation":"GENCODE_V19_GRCH37_POSITIONAL_ONLY_NOT_CAUSAL"})
    if not out:raise ValueError("No nearby GENCODE genes")
    with (root/"G0022_PILOT_TOP_SNP_NEAREST_GENES.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=list(out[0]))
        w.writeheader();w.writerows(out)
    summary={"pilot_signals":len(gate),"gene_link_rows":len(out),
       "unique_nearby_gene_ids":len({r["gene_id"] for r in out}),
       "signals_with_nearby_annotated_gene":len({r["signal"] for r in out}),
       "ancestry":"EAS_GWAS_PLUS_1KG_EAS_LD",
       "status":"PROXIMITY_ONLY_NO_CAUSAL_GENE_CLAIM"}
    (root/"G0022_PILOT_NEAREST_GENES_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    return summary
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--base",type=Path,default=BASE)
    a=p.parse_args()
    build(a.root,a.base)
