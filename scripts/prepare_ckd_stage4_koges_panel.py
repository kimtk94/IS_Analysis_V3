#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, gzip, json
from collections import defaultdict
from pathlib import Path

def open_text(p):
    p=Path(p)
    return gzip.open(p,"rt",encoding="utf-8",errors="replace",newline="") if p.suffix==".gz" else p.open("r",encoding="utf-8",errors="replace",newline="")

def read_tsv(p):
    with open_text(p) as f: return list(csv.DictReader(f,delimiter="\t"))

def write_tsv(p,rows):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    fields=list(rows[0]) if rows else []
    with p.open("w",newline="",encoding="utf-8") as f:
        if not fields: return
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
        w.writeheader(); w.writerows(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--stage1-summary",type=Path,required=True)
    ap.add_argument("--instrument-qc",type=Path,required=True)
    ap.add_argument("--stage3b-integrated",type=Path,required=True)
    ap.add_argument("--output-root",type=Path,required=True)
    a=ap.parse_args()

    s1=read_tsv(a.stage1_summary)
    inst=read_tsv(a.instrument_qc)
    s3=read_tsv(a.stage3b_integrated)
    candidates={r["gene_symbol"].upper():r for r in s3}

    # Stage1 summary retains the pre-specified strongest-F EUR eGFR anchor.
    s1c=[r for r in s1 if r.get("gene_symbol","").upper() in candidates]
    by_protein=defaultdict(list)
    for r in inst: by_protein[r.get("protein_id","")].append(r)

    anchors=[]; allcis=[]
    for r in s1c:
        gene=r["gene_symbol"].upper(); protein=r["protein_id"]; anchor=r.get("eur_egfr_rsid","")
        rows=[x for x in by_protein.get(protein,[]) if x.get("instrument_qc")=="ok"]
        for x in rows:
            role="anchor" if x.get("rsid")==anchor else "sensitivity"
            base={
              "gene_symbol":gene,"protein_id":protein,"panel_role":role,
              "rsid":x.get("rsid",""),"variant_id":x.get("variant_id",""),
              "chrom_hg19":x.get("chrom_hg19",""),"pos_hg19":x.get("pos_hg19",""),
              "ref_allele":x.get("ref_allele",""),"effect_allele_alt":x.get("alt_allele",""),
              "beta_cond_alt":x.get("beta_cond",""),"se_cond":x.get("se_cond",""),
              "alt_freq_eur":x.get("alt_freq",""),"pip":x.get("pip",""),
              "f_stat":x.get("f_stat",""),
              "stage2_class":candidates[gene].get("stage2_class",""),
              "localization_pattern":candidates[gene].get("localization_pattern",""),
              "top_compartment":candidates[gene].get("top_compartment",""),
              "top_cell_type":candidates[gene].get("top_cell_type",""),
            }
            allcis.append(base)
            if role=="anchor": anchors.append(base)

    genes_expected=sorted(candidates)
    genes_anchor=sorted({r["gene_symbol"] for r in anchors})
    missing=sorted(set(genes_expected)-set(genes_anchor))

    a.output_root.mkdir(parents=True,exist_ok=True)
    write_tsv(a.output_root/"STAGE4_KOGES_ANCHOR_PANEL.tsv",sorted(anchors,key=lambda x:x["gene_symbol"]))
    write_tsv(a.output_root/"STAGE4_KOGES_ALL_CIS_PANEL.tsv",sorted(allcis,key=lambda x:(x["gene_symbol"],0 if x["panel_role"]=="anchor" else 1,-float(x["f_stat"] or 0))))
    meta={
      "stage":"CKD Stage 4 KoGES individual-level validation preparation",
      "candidate_genes":genes_expected,
      "anchor_genes":genes_anchor,
      "missing_anchor_genes":missing,
      "primary_rule":"Use the pre-specified strongest-F EUR cis-pQTL anchor from Stage 1; do not select a new KoGES SNP by outcome association.",
      "effect_allele":"ALT from UKB-PPP ST16 / Stage1 instrument_qc",
      "sensitivity_rule":"Additional strong cis instruments are exported separately; multi-SNP analyses require KoGES ancestry-matched LD handling.",
      "restricted_data_note":"KoGES individual-level genotype/phenotype data must remain local and must not be committed or synced to public/shared storage."
    }
    (a.output_root/"STAGE4_KOGES_PANEL_PROVENANCE.json").write_text(json.dumps(meta,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(meta,indent=2))
    print("CKD_STAGE4_KOGES_PANEL_PASS")

if __name__=="__main__": main()
