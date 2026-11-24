#!/usr/bin/env python3
"""Bridge existing CKD Stage3A/3B evidence into MASTER localization summary.

No new biological inference is performed:
- tissue evidence is "present" only when Stage3A contains a positive kidney
  pQTL/eQTL/supplement hit for the candidate;
- cell evidence is "present" only when Stage3B has a non-empty top cell type
  with positive reported expression;
- spatial evidence is always absent in this bridge.

Therefore this bridge can emit TISSUE_CELL or SINGLE_LAYER, but never
TISSUE_CELL_SPATIAL.
"""
from __future__ import annotations
import argparse,csv,math
from pathlib import Path

def read_tsv(p):
    with Path(p).open("r",encoding="utf-8",newline="") as f:
        yield from csv.DictReader(f,delimiter="\t")

def fnum(x):
    try:
        v=float(str(x).strip())
    except Exception:
        return None
    return v if math.isfinite(v) else None

def positive(r, keys):
    for k in keys:
        v=fnum(r.get(k))
        if v is not None and v>0:
            return True
    return False

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--stage3a",type=Path,required=True)
    ap.add_argument("--stage3b",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()

    s3a={r["gene_symbol"].upper():r for r in read_tsv(a.stage3a)}
    s3b={r["gene_symbol"].upper():r for r in read_tsv(a.stage3b)}
    genes=sorted(set(s3a)|set(s3b))
    rows=[]
    tissue_keys=[
      "hirohama_maintext_kidney_pqtl",
      "kidney_pqtl_significant_supp_hits",
      "kidney_eqtl_same_study_supp_hits",
      "kidney_qtl_egfr_coloc_supp_hits",
      "kidney_eqtl_meta686_hits",
      "kidney_eqtl_tubule356_hits",
      "kidney_eqtl_glomerulus303_hits",
    ]
    for g in genes:
        ta=s3a.get(g,{})
        cb=s3b.get(g,{})
        tissue=positive(ta,tissue_keys)
        top_cell=str(cb.get("top_cell_type","")).strip()
        top_expr=fnum(cb.get("top_nCPM"))
        cell=bool(top_cell) and top_expr is not None and top_expr>0
        if tissue and cell:
            cls="TISSUE_CELL"
        elif tissue or cell:
            cls="SINGLE_LAYER"
        else:
            cls="NO_LOCALIZATION"
        rows.append({
          "gene_symbol":g,
          "tissue_source_n":1 if tissue else 0,
          "top_tissue":"kidney" if tissue else "",
          "top_tissue_expression":"",
          "top_tissue_share":"",
          "max_tissue_specificity":"",
          "cell_source_n":1 if cell else 0,
          "top_cell_type":top_cell if cell else "",
          "top_cell_compartment":cb.get("top_compartment","") if cell else "",
          "top_cell_expression":top_expr if cell else "",
          "top_second_ratio":cb.get("top_second_ratio","") if cell else "",
          "max_cell_specificity":cb.get("tau_specificity","") if cell else "",
          "positive_celltype_n":cb.get("positive_celltype_n","") if cell else 0,
          "spatial_source_n":0,
          "top_spatial_region":"",
          "top_spatial_expression":"",
          "top_spatial_enrichment":"",
          "top_spatial_spot_fraction":"",
          "localization_evidence_class":cls,
          "bridge_note":"Existing CKD Stage3A/3B only; no spatial evidence inferred."
        })
    fields=list(rows[0]) if rows else []
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
        w.writeheader();w.writerows(rows)
    print(f"CKD_MASTER_LOCALIZATION_BRIDGE_PASS genes={len(rows)} output={a.output}")

if __name__=="__main__":
    main()
