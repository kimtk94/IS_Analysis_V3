#!/usr/bin/env python3
"""Select EUR-confirmed candidates before East-Asian replication.

This stage does not label a protein as fully replicated. It only summarizes
EUR evidence using the locked project rules:
- MR FDR < 0.05,
- favorable direction in metabolic domains,
- coloc PP.H4 >= 0.80 at the default p12,
- multi-domain priority when >=2 metabolic domains are supported.

EAS evidence is added later and must not be inferred here.
"""

from __future__ import annotations

import csv
import json
import math
import os
from collections import defaultdict
from pathlib import Path

ROOT = Path(os.environ.get("IS_ANALYSIS_ROOT", "/srv/is-analysis"))
MR = ROOT / "results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D5_PRIMARY_MR_SUMMARY.tsv"
COLOC = ROOT / "results/metabolic_resilience/stage3_coloc/STAGE3E3_COLOC_SUMMARY.tsv"
OUTDIR = ROOT / "results/metabolic_resilience/stage4_eas"
OUT = OUTDIR / "STAGE4C_EUR_CANDIDATE_GATE.tsv"
SUMMARY = OUTDIR / "STAGE4C_EUR_CANDIDATE_GATE.json"

OUTDIR.mkdir(parents=True, exist_ok=True)

MR_FDR = 0.05
COLOC_H4 = 0.80
COLOC_RATIO = 5.0


def fnum(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None


def first(row, *names):
    for name in names:
        if name in row and row[name] not in ("", None, "NA"):
            return row[name]
    return None


if not MR.exists() or not COLOC.exists():
    payload = {
        "status": "HOLD_MISSING_STAGE3_INPUT",
        "mr_exists": MR.exists(),
        "coloc_exists": COLOC.exists(),
    }
    SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    raise SystemExit(2)

with MR.open("r", encoding="utf-8", newline="") as f:
    mr_rows = list(csv.DictReader(f, delimiter="\t"))

with COLOC.open("r", encoding="utf-8", newline="") as f:
    coloc_rows = list(csv.DictReader(f, delimiter="\t"))

# Strong coloc at default p12 is already the intended Stage3E3 summary view,
# but keep direct numeric checks so the gate is auditable.
coloc_by_gene_trait = {}
for r in coloc_rows:
    gene = first(r, "gene", "gene_symbol")
    trait = first(r, "trait")
    if not gene or not trait:
        continue
    h4 = fnum(first(r, "PP.H4", "PP_H4", "PP.H4.abf"))
    h3 = fnum(first(r, "PP.H3", "PP_H3", "PP.H3.abf"))
    if h4 is None:
        continue
    ratio = None if h3 in (None, 0) else h4 / h3
    coloc_by_gene_trait[(gene, trait)] = {
        "h4": h4,
        "h3": h3,
        "h4_h3_ratio": ratio,
        "coloc_strong": int(h4 >= COLOC_H4),
        "ratio_gt_5": int(ratio is not None and ratio > COLOC_RATIO),
    }

by_gene = defaultdict(list)

for r in mr_rows:
    if r.get("mode") not in ("", None, "primary"):
        continue

    gene = first(r, "gene", "gene_symbol")
    trait = first(r, "trait")
    domain = first(r, "domain") or ""
    if not gene or not trait:
        continue

    p = fnum(first(r, "p", "primary_p"))
    fdr = fnum(first(r, "FDR_mode", "primary_fdr", "fdr"))
    beta = fnum(first(r, "beta", "primary_beta"))
    fav = first(r, "direction_favorable")
    fav = None if fav is None else str(fav).strip() in {"1", "TRUE", "True", "true"}

    coloc = coloc_by_gene_trait.get((gene, trait), {})
    mr_fdr_sig = fdr is not None and fdr < MR_FDR
    metabolic = domain != "disease_validation"

    by_gene[gene].append({
        "trait": trait,
        "domain": domain,
        "beta": beta,
        "p": p,
        "fdr": fdr,
        "favorable": fav,
        "metabolic": metabolic,
        "mr_fdr_sig": mr_fdr_sig,
        "coloc_h4": coloc.get("h4"),
        "coloc_h3": coloc.get("h3"),
        "coloc_ratio": coloc.get("h4_h3_ratio"),
        "coloc_strong": coloc.get("coloc_strong", 0) == 1,
        "coloc_ratio_gt_5": coloc.get("ratio_gt_5", 0) == 1,
    })

rows = []
for gene in sorted(by_gene):
    x = by_gene[gene]
    metabolic = [r for r in x if r["metabolic"]]
    strong_pairs = [
        r for r in metabolic
        if r["mr_fdr_sig"] and r["favorable"] is True and r["coloc_strong"]
    ]
    domains = sorted({r["domain"] for r in strong_pairs if r["domain"]})
    opposite_fdr = [
        r for r in metabolic
        if r["mr_fdr_sig"] and r["favorable"] is False
    ]

    if len(domains) >= 2 and not opposite_fdr:
        gate = "EUR_MULTI_DOMAIN_PRIORITY"
    elif strong_pairs and not opposite_fdr:
        gate = "EUR_STRONG_SINGLE_OR_SAME_DOMAIN"
    elif any(r["mr_fdr_sig"] for r in metabolic):
        gate = "REVIEW_MR_WITHOUT_STRONG_COLOC_OR_DIRECTION"
    else:
        gate = "NO_EUR_STRONG_EVIDENCE"

    rows.append({
        "gene": gene,
        "metabolic_traits_available": len(metabolic),
        "mr_fdr_lt_0p05": sum(r["mr_fdr_sig"] for r in metabolic),
        "strong_mr_coloc_favorable_pairs": len(strong_pairs),
        "strong_domains": len(domains),
        "strong_domain_names": ";".join(domains),
        "strong_trait_names": ";".join(sorted(r["trait"] for r in strong_pairs)),
        "opposite_direction_fdr_traits": ";".join(
            sorted(r["trait"] for r in opposite_fdr)
        ),
        "coloc_ratio_gt5_pairs": sum(r["coloc_ratio_gt_5"] for r in strong_pairs),
        "stage4_gate": gate,
    })

fields = [
    "gene",
    "metabolic_traits_available",
    "mr_fdr_lt_0p05",
    "strong_mr_coloc_favorable_pairs",
    "strong_domains",
    "strong_domain_names",
    "strong_trait_names",
    "opposite_direction_fdr_traits",
    "coloc_ratio_gt5_pairs",
    "stage4_gate",
]
with OUT.open("w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)

payload = {
    "status": "PASS",
    "genes": len(rows),
    "eur_multi_domain_priority": sum(r["stage4_gate"] == "EUR_MULTI_DOMAIN_PRIORITY" for r in rows),
    "eur_strong_single_or_same_domain": sum(r["stage4_gate"] == "EUR_STRONG_SINGLE_OR_SAME_DOMAIN" for r in rows),
    "mr_fdr_threshold": MR_FDR,
    "coloc_h4_threshold": COLOC_H4,
    "coloc_h4_h3_ratio_flag": COLOC_RATIO,
    "important": "EAS replication is not inferred by this stage.",
}
SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
print(json.dumps(payload, indent=2))
