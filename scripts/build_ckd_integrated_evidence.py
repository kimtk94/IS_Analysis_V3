#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Any


ROOT = Path("/srv/is-analysis")
RES = ROOT / "results/ckd"
OUT = RES / "integrated_evidence"

GENES = [
    "ACP1",
    "CPVL",
    "F12",
    "GSTA1",
    "GSTA3",
    "HLA-E",
    "INHBC",
    "SDCCAG8",
    "UMOD",
]


def read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        print(f"[WARN] missing: {path}")
        return []

    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: List[Dict[str, Any]], fields: List[str]):
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=fields,
            extrasaction="ignore",
        )
        w.writeheader()

        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})


def norm_gene(x: str) -> str:
    return (x or "").strip().upper()


def fnum(x):
    try:
        y = float(x)
        if math.isfinite(y):
            return y
    except Exception:
        pass
    return None


def find_first(paths):
    for p in paths:
        if p.exists():
            return p
    return None


def by_gene(rows):
    out = {}
    for r in rows:
        g = norm_gene(r.get("gene_symbol", ""))
        if g:
            out[g] = r
    return out


def by_gene_pheno(rows):
    out = {}
    for r in rows:
        g = norm_gene(r.get("gene_symbol", ""))
        p = (r.get("phenotype", "") or "").strip()
        if g and p:
            out[(g, p)] = r
    return out


# ============================================================
# Source discovery
# ============================================================

stage3b_path = RES / "stage3b_celltype/STAGE3B_INTEGRATED_EVIDENCE.tsv"
local_path = RES / "stage3b_celltype/STAGE3B_KIDNEY_LOCALIZATION.tsv"
spec_path = RES / "stage3b_celltype/STAGE3B_KIDNEY_SPECIFICITY.tsv"
koges_path = RES / "stage4_koges/STAGE4_KOGES_ANCHOR_PANEL.tsv"

eas_anchor_path = RES / "stage4_eas/anchor_mr/EAS_ANCHOR_WALD_MR.tsv"
eas_integrated_path = RES / "stage4_eas/anchor_mr/EAS_ANCHOR_INTEGRATED_EVIDENCE.tsv"

egfr_coloc_path = (
    RES / "stage4_eas/regional/eGFR/coloc_results/EAS_COLOC_DEFAULT.tsv"
)
bun_coloc_path = (
    RES / "stage4_eas/regional/BUN/coloc_results/EAS_COLOC_DEFAULT.tsv"
)

dual_path = (
    RES
    / "stage4_eas/dual_ld_susie/susie_results/DUAL_LD_SUSIE_DEFAULT.tsv"
)

f12_status_path = (
    RES / "stage4_eas/f12_final_diagnostic/F12_ANALYSIS_STATUS.tsv"
)


# ============================================================
# Read sources
# ============================================================

stage3b = by_gene(read_tsv(stage3b_path))
localization = by_gene(read_tsv(local_path))
specificity = by_gene(read_tsv(spec_path))
koges = by_gene(read_tsv(koges_path))

eas_anchor_rows = read_tsv(eas_anchor_path)
eas_integrated_rows = read_tsv(eas_integrated_path)

eas_anchor = by_gene_pheno(eas_anchor_rows)
eas_integrated = by_gene(eas_integrated_rows)

egfr_coloc = by_gene(read_tsv(egfr_coloc_path))
bun_coloc = by_gene(read_tsv(bun_coloc_path))

dual = by_gene_pheno(read_tsv(dual_path))

f12_status_rows = read_tsv(f12_status_path)


# ============================================================
# Generic anchor-field extraction
# ============================================================

def pick(row, candidates, default=""):
    if not row:
        return default

    for k in candidates:
        if k in row and row[k] not in ("", None):
            return row[k]

    return default


def anchor_for(gene, pheno):
    r = eas_anchor.get((gene, pheno), {})

    return {
        "beta": pick(
            r,
            [
                "wald_beta",
                "beta_wald",
                "beta_mr",
                "mr_beta",
                "beta",
            ],
        ),
        "se": pick(
            r,
            [
                "wald_se",
                "se_wald",
                "se_mr",
                "mr_se",
                "se",
            ],
        ),
        "p": pick(
            r,
            [
                "wald_p",
                "p_wald",
                "p_mr",
                "mr_p",
                "p",
                "pvalue",
            ],
        ),
        "outcome_beta": pick(
            r,
            [
                "beta_outcome",
                "outcome_beta",
            ],
        ),
        "outcome_se": pick(
            r,
            [
                "se_outcome",
                "outcome_se",
            ],
        ),
        "outcome_p": pick(
            r,
            [
                "p_outcome",
                "outcome_p",
            ],
        ),
    }


# ============================================================
# Evidence interpretation
# Descriptive categories only — no score/rank.
# ============================================================

def classify_eas_coloc(row):
    if not row:
        return "missing"

    h3 = fnum(row.get("PP.H3"))
    h4 = fnum(row.get("PP.H4"))

    if h3 is None or h4 is None:
        return "uninformative"

    if h4 >= 0.80:
        return "shared_signal_supported"

    if h3 >= 0.80 and h4 < 0.20:
        return "distinct_signal_supported"

    if h4 >= 0.50:
        return "shared_signal_suggestive"

    return "uninformative"


def classify_dual(row):
    if not row:
        return "missing"

    status = (row.get("comparison_status") or "").strip()

    if status:
        return status

    n_out = fnum(row.get("outcome_credible_sets"))
    max_h4 = fnum(row.get("max_PP.H4"))

    if n_out == 0:
        return "uninformative_no_outcome_cs"

    if max_h4 is not None and max_h4 >= 0.80:
        return "shared_signal_supported"

    if max_h4 is not None and max_h4 >= 0.50:
        return "shared_signal_suggestive"

    return "unresolved"


def overall_pattern(row):
    eur_class = row.get("eur_stage2_class", "")

    e1 = row.get("eas_egfr_coloc_status", "")
    e2 = row.get("eas_bun_coloc_status", "")

    d1 = row.get("dual_egfr_status", "")
    d2 = row.get("dual_bun_status", "")

    if row["gene_symbol"] == "F12":
        return "EUR_EVIDENCE_EAS_DISTINCT_FINE_MAPPING_TECHNICAL_UNRESOLVED"

    eas_shared = any(
        x == "shared_signal_supported"
        for x in [e1, e2, d1, d2]
    )

    eas_distinct = any(
        x == "distinct_signal_supported"
        for x in [e1, e2, d1, d2]
    )

    eas_uninformative = all(
        x in {
            "",
            "missing",
            "uninformative",
            "unresolved",
            "uninformative_no_outcome_cs",
        }
        for x in [e1, e2, d1, d2]
    )

    # EUR shared evidence must be explicit.
    # Do NOT use substring matching because "nonshared"
    # also contains the string "shared".
    eur_class_norm = eur_class.strip().lower()

    eur_explicit_shared_classes = {
        "core_shared_signal",
    }

    eur_h4 = fnum(row.get("eur_abf_h4"))
    eur_susie_h4 = fnum(row.get("eur_susie_max_h4"))

    eur_shared = (
        eur_class_norm in eur_explicit_shared_classes
        or (
            eur_h4 is not None
            and eur_susie_h4 is not None
            and eur_h4 >= 0.80
            and eur_susie_h4 >= 0.80
        )
    )

    if eur_shared and eas_shared:
        return "EUR_SHARED_EAS_SHARED"

    if eur_shared and eas_distinct:
        return "EUR_SHARED_EAS_DISTINCT"

    if eur_shared and eas_uninformative:
        return "EUR_SHARED_EAS_UNINFORMATIVE"

    if eas_shared:
        return "EAS_SHARED_SIGNAL_SUPPORTED"

    if eas_distinct:
        return "EAS_DISTINCT_SIGNAL_SUPPORTED"

    return "MIXED_OR_UNRESOLVED"


# ============================================================
# Build wide matrix
# ============================================================

matrix = []

for gene in GENES:

    s3 = stage3b.get(gene, {})
    loc = localization.get(gene, {})
    sp = specificity.get(gene, {})
    kg = koges.get(gene, {})

    eg = egfr_coloc.get(gene, {})
    bu = bun_coloc.get(gene, {})

    de = dual.get((gene, "eGFR"), {})
    db = dual.get((gene, "BUN"), {})

    am_e = anchor_for(gene, "eGFR")
    am_b = anchor_for(gene, "BUN")

    row = {
        # identity
        "gene_symbol": gene,

        # EUR discovery
        "eur_stage2_class":
            s3.get("stage2_class", ""),

        "eur_comparison_status":
            s3.get("stage2_comparison_status", ""),

        "eur_abf_h4":
            s3.get("stage2_abf_h4", ""),

        "eur_susie_max_h4":
            s3.get("stage2_susie_max_h4", ""),

        # kidney tissue
        "kidney_pqtl_maintext":
            s3.get("hirohama_maintext_kidney_pqtl", ""),

        "kidney_eqtl_meta686_hits":
            s3.get("kidney_eqtl_meta686_hits", ""),

        "kidney_eqtl_tubule356_hits":
            s3.get("kidney_eqtl_tubule356_hits", ""),

        "kidney_eqtl_glomerulus303_hits":
            s3.get("kidney_eqtl_glomerulus303_hits", ""),

        # cell localization
        "top_cell_type":
            loc.get(
                "top_kidney_cell_type",
                s3.get("top_cell_type", ""),
            ),

        "top_nCPM":
            loc.get(
                "top_nCPM",
                s3.get("top_nCPM", ""),
            ),

        "top_compartment":
            sp.get(
                "top_compartment",
                s3.get("top_compartment", ""),
            ),

        "tau_specificity":
            sp.get(
                "tau_specificity",
                s3.get("tau_specificity", ""),
            ),

        "localization_pattern":
            sp.get(
                "localization_pattern",
                s3.get("localization_pattern", ""),
            ),

        # EAS anchor MR
        "eas_egfr_wald_beta": am_e["beta"],
        "eas_egfr_wald_se": am_e["se"],
        "eas_egfr_wald_p": am_e["p"],

        "eas_bun_wald_beta": am_b["beta"],
        "eas_bun_wald_se": am_b["se"],
        "eas_bun_wald_p": am_b["p"],

        # EAS ABF coloc
        "eas_egfr_abf_h3": eg.get("PP.H3", ""),
        "eas_egfr_abf_h4": eg.get("PP.H4", ""),
        "eas_egfr_coloc_status": classify_eas_coloc(eg),

        "eas_bun_abf_h3": bu.get("PP.H3", ""),
        "eas_bun_abf_h4": bu.get("PP.H4", ""),
        "eas_bun_coloc_status": classify_eas_coloc(bu),

        # dual LD SuSiE
        "dual_egfr_pqtl_cs":
            de.get("pqtl_credible_sets", ""),

        "dual_egfr_outcome_cs":
            de.get("outcome_credible_sets", ""),

        "dual_egfr_max_h4":
            de.get("max_PP.H4", ""),

        "dual_egfr_status":
            classify_dual(de),

        "dual_bun_pqtl_cs":
            db.get("pqtl_credible_sets", ""),

        "dual_bun_outcome_cs":
            db.get("outcome_credible_sets", ""),

        "dual_bun_max_h4":
            db.get("max_PP.H4", ""),

        "dual_bun_status":
            classify_dual(db),

        # KoGES anchor
        "koges_rsid":
            kg.get("rsid", ""),

        "koges_variant_id":
            kg.get("variant_id", ""),

        "koges_effect_allele":
            kg.get("effect_allele_alt", ""),

        "koges_beta_pqtl":
            kg.get("beta_cond_alt", ""),

        "koges_se_pqtl":
            kg.get("se_cond", ""),

        "koges_eaf_eur":
            kg.get("alt_freq_eur", ""),

        "koges_pip":
            kg.get("pip", ""),

        "koges_f_stat":
            kg.get("f_stat", ""),

        # explicit technical status
        "technical_status":
            (
                "F12_PQTL_SUSIE_TECHNICALLY_UNRESOLVED"
                if gene == "F12"
                else "NO_FROZEN_TECHNICAL_EXCEPTION"
            ),
    }

    row["evidence_pattern"] = overall_pattern(row)

    matrix.append(row)


FIELDS = [
    "gene_symbol",

    "eur_stage2_class",
    "eur_comparison_status",
    "eur_abf_h4",
    "eur_susie_max_h4",

    "kidney_pqtl_maintext",
    "kidney_eqtl_meta686_hits",
    "kidney_eqtl_tubule356_hits",
    "kidney_eqtl_glomerulus303_hits",

    "top_cell_type",
    "top_nCPM",
    "top_compartment",
    "tau_specificity",
    "localization_pattern",

    "eas_egfr_wald_beta",
    "eas_egfr_wald_se",
    "eas_egfr_wald_p",

    "eas_bun_wald_beta",
    "eas_bun_wald_se",
    "eas_bun_wald_p",

    "eas_egfr_abf_h3",
    "eas_egfr_abf_h4",
    "eas_egfr_coloc_status",

    "eas_bun_abf_h3",
    "eas_bun_abf_h4",
    "eas_bun_coloc_status",

    "dual_egfr_pqtl_cs",
    "dual_egfr_outcome_cs",
    "dual_egfr_max_h4",
    "dual_egfr_status",

    "dual_bun_pqtl_cs",
    "dual_bun_outcome_cs",
    "dual_bun_max_h4",
    "dual_bun_status",

    "koges_rsid",
    "koges_variant_id",
    "koges_effect_allele",
    "koges_beta_pqtl",
    "koges_se_pqtl",
    "koges_eaf_eur",
    "koges_pip",
    "koges_f_stat",

    "technical_status",
    "evidence_pattern",
]

write_tsv(
    OUT / "CKD_9GENE_EVIDENCE_MATRIX.tsv",
    matrix,
    FIELDS,
)


# ============================================================
# Long-format evidence
# ============================================================

long_rows = []

def add(gene, domain, evidence, value, status=""):
    long_rows.append({
        "gene_symbol": gene,
        "domain": domain,
        "evidence": evidence,
        "value": value,
        "status": status,
    })


for r in matrix:

    g = r["gene_symbol"]

    add(
        g,
        "EUR",
        "ABF_coloc_PP.H4",
        r["eur_abf_h4"],
        r["eur_comparison_status"],
    )

    add(
        g,
        "EUR",
        "SuSiE_max_PP.H4",
        r["eur_susie_max_h4"],
        r["eur_stage2_class"],
    )

    add(
        g,
        "EAS_eGFR",
        "ABF_coloc_PP.H4",
        r["eas_egfr_abf_h4"],
        r["eas_egfr_coloc_status"],
    )

    add(
        g,
        "EAS_BUN",
        "ABF_coloc_PP.H4",
        r["eas_bun_abf_h4"],
        r["eas_bun_coloc_status"],
    )

    add(
        g,
        "DualLD_eGFR",
        "SuSiE_max_PP.H4",
        r["dual_egfr_max_h4"],
        r["dual_egfr_status"],
    )

    add(
        g,
        "DualLD_BUN",
        "SuSiE_max_PP.H4",
        r["dual_bun_max_h4"],
        r["dual_bun_status"],
    )

    add(
        g,
        "Kidney",
        "top_cell_type",
        r["top_cell_type"],
        r["localization_pattern"],
    )

    add(
        g,
        "KoGES",
        "anchor_rsid",
        r["koges_rsid"],
        "pending_individual_level_validation",
    )


write_tsv(
    OUT / "CKD_9GENE_EVIDENCE_LONG.tsv",
    long_rows,
    [
        "gene_symbol",
        "domain",
        "evidence",
        "value",
        "status",
    ],
)


# ============================================================
# Interpretation table
# ============================================================

interpretation = []

for r in matrix:

    interpretation.append({
        "gene_symbol": r["gene_symbol"],
        "eur_pattern": r["eur_stage2_class"],
        "eas_egfr_pattern": r["eas_egfr_coloc_status"],
        "eas_bun_pattern": r["eas_bun_coloc_status"],
        "dual_egfr_pattern": r["dual_egfr_status"],
        "dual_bun_pattern": r["dual_bun_status"],
        "kidney_localization": r["top_cell_type"],
        "technical_status": r["technical_status"],
        "integrated_pattern": r["evidence_pattern"],
        "koges_next_step": (
            "individual_level_longitudinal_validation"
            if r["koges_rsid"]
            else "anchor_availability_check"
        ),
    })


write_tsv(
    OUT / "CKD_9GENE_INTERPRETATION.tsv",
    interpretation,
    [
        "gene_symbol",
        "eur_pattern",
        "eas_egfr_pattern",
        "eas_bun_pattern",
        "dual_egfr_pattern",
        "dual_bun_pattern",
        "kidney_localization",
        "technical_status",
        "integrated_pattern",
        "koges_next_step",
    ],
)


# ============================================================
# QA
# ============================================================

qa = {
    "expected_genes": GENES,
    "n_expected": len(GENES),
    "n_matrix": len(matrix),
    "n_long": len(long_rows),
    "missing_stage3b": [
        g for g in GENES if g not in stage3b
    ],
    "missing_koges_anchor": [
        g for g in GENES if g not in koges
    ],
    "missing_eas_egfr_coloc": [
        g for g in GENES if g not in egfr_coloc
    ],
    "missing_eas_bun_coloc": [
        g for g in GENES if g not in bun_coloc
    ],
    "missing_dual_egfr": [
        g for g in GENES if (g, "eGFR") not in dual
    ],
    "missing_dual_bun": [
        g for g in GENES if (g, "BUN") not in dual
    ],
}

qa["pass"] = (
    qa["n_matrix"] == 9
    and not qa["missing_stage3b"]
    and not qa["missing_koges_anchor"]
    and not qa["missing_eas_egfr_coloc"]
    and not qa["missing_eas_bun_coloc"]
)

with (OUT / "CKD_9GENE_EVIDENCE_QA.json").open("w") as f:
    json.dump(qa, f, indent=2)


# ============================================================
# Human-readable summary
# ============================================================

with (OUT / "CKD_9GENE_EVIDENCE_SUMMARY.txt").open("w") as f:

    for r in matrix:
        f.write(
            f"{r['gene_symbol']}\t"
            f"{r['evidence_pattern']}\t"
            f"EUR_H4={r['eur_abf_h4']}\t"
            f"EAS_eGFR_H4={r['eas_egfr_abf_h4']}\t"
            f"EAS_BUN_H4={r['eas_bun_abf_h4']}\t"
            f"cell={r['top_cell_type']}\t"
            f"KoGES={r['koges_rsid']}\n"
        )


print("============================================================")
print("CKD 9-GENE INTEGRATED EVIDENCE")
print("============================================================")

for r in matrix:
    print(
        r["gene_symbol"],
        r["evidence_pattern"],
        "EUR_H4=" + str(r["eur_abf_h4"]),
        "EAS_eGFR_H4=" + str(r["eas_egfr_abf_h4"]),
        "EAS_BUN_H4=" + str(r["eas_bun_abf_h4"]),
        "KoGES=" + str(r["koges_rsid"]),
        sep="\t",
    )

print()
print("QA:")
print(json.dumps(qa, indent=2))

print()
print("OUTPUT =", OUT)
