#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv
from collections import defaultdict
from pathlib import Path

def num(x):
    try: return float(x)
    except Exception: return 0.0

def compartment(cell):
    x=(cell or "").lower()
    groups={
      "immune":["macrophage","monocyte","neutrophil","t-cell","t cell","b-cell","b cell","cdc","dendritic","mast","nk cell","plasma cell"],
      "vascular":["endothelial","vascular smooth muscle","pericyte"],
      "renal_epithelial":["proximal tubule","distal tubule","loop of henle","collecting duct","connecting tubule","papillary tip","podocyte","urothelial","epithelial"],
      "stromal":["fibroblast","mesangial","stromal"],
    }
    for group, terms in groups.items():
        if any(t in x for t in terms): return group
    return "other"

def read_tsv(p):
    with Path(p).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(p, rows):
    p=Path(p); p.parent.mkdir(parents=True, exist_ok=True)
    fields=list(rows[0]) if rows else []
    with p.open("w", newline="", encoding="utf-8") as f:
        if not fields: return
        w=csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

def summarize(rows):
    gct=defaultdict(lambda: defaultdict(float))
    for r in rows:
        gene=r["gene_symbol"].upper()
        ct=r["cell_type"]
        val=num(r.get("expression", r.get("nCPM",0)))
        gct[gene][ct]=max(gct[gene][ct], val)
    all_ct=sorted({ct for d in gct.values() for ct in d})
    out=[]
    for gene in sorted(gct):
        vals=sorted([(ct,gct[gene].get(ct,0.0)) for ct in all_ct], key=lambda x:x[1], reverse=True)
        top_ct,top=vals[0]; second=vals[1][1] if len(vals)>1 else 0.0
        total=sum(v for _,v in vals)
        ratio=top/second if second>0 else (float("inf") if top>0 else 0.0)
        tau=sum(1-v/top for _,v in vals)/(len(vals)-1) if top>0 and len(vals)>1 else 0.0
        comp_max=defaultdict(float); comp_cell={}
        for ct,v in vals:
            c=compartment(ct)
            if v>comp_max[c]: comp_max[c]=v; comp_cell[c]=ct
        epi=comp_max["renal_epithelial"]; imm=comp_max["immune"]
        epi_imm=epi/imm if imm>0 else (float("inf") if epi>0 else 0.0)
        pattern="very_low_expression" if top<1 else ("focused" if ratio>=5 else ("moderately_focused" if ratio>=2 else "broad"))
        top_comp=max(comp_max,key=comp_max.get) if comp_max else ""
        out.append({
          "gene_symbol":gene,"top_cell_type":top_ct,"top_nCPM":top,"second_nCPM":second,
          "top_second_ratio":ratio,"top_expression_share":top/total if total else 0.0,
          "tau_specificity":tau,"positive_celltype_n":sum(v>0 for _,v in vals),
          "top_compartment":top_comp,"renal_epithelial_max":epi,
          "renal_epithelial_top_cell":comp_cell.get("renal_epithelial",""),
          "immune_max":imm,"immune_top_cell":comp_cell.get("immune",""),
          "vascular_max":comp_max["vascular"],"stromal_max":comp_max["stromal"],
          "epithelial_to_immune_ratio":epi_imm,"localization_pattern":pattern,
        })
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--cluster-expression",type=Path,required=True)
    ap.add_argument("--stage3a-evidence",type=Path,required=True)
    ap.add_argument("--output-root",type=Path,required=True)
    a=ap.parse_args()
    summary=summarize(read_tsv(a.cluster_expression))
    write_tsv(a.output_root/"STAGE3B_KIDNEY_SPECIFICITY.tsv",summary)
    smap={x["gene_symbol"]:x for x in summary}
    merged=[]
    for r in read_tsv(a.stage3a_evidence):
        x=dict(r); x.update(smap.get(r["gene_symbol"].upper(),{})); merged.append(x)
    write_tsv(a.output_root/"STAGE3B_INTEGRATED_EVIDENCE.tsv",merged)
    for x in summary:
        print(x["gene_symbol"],x["top_cell_type"],x["top_nCPM"],x["top_second_ratio"],x["tau_specificity"],x["localization_pattern"],sep="\t")
    print("CKD_STAGE3B_REFINED_PASS")

if __name__=="__main__": main()
