#!/usr/bin/env python3
"""Stratified same-source reproducibility of TM DEX pseudobulk edgeR effects.

Uses REAL previously fitted public aggregate gene-level outputs; no donor-level
counts, original genotypes, or repeat raw 10x matrix downloads. These two
sensitivity contrasts share source donors/libraries and are NOT independent
biological replications. Binning on baseline magnitude is descriptive.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

CONTRASTS = [
    ("BeamA classification score ≥0.75 vs ≥0.50 (57 paired donors)", "score"),
    ("BeamA vs scPred TM_Fibro (same 55 paired donors)", "label"),
]
EFFECT_CUTS = [(0, .25, "|LFC| < 0.25"), (.25, 1, "0.25 ≤ |LFC| < 1"),
               (1, 3, "1 ≤ |LFC| < 3"), (3, math.inf, "|LFC| ≥ 3")]
CORE = ("MYOC", "TGFBI", "MMP2")
PANEL = ("MYOC", "TGFBI", "MMP2", "ACTA2", "FKBP5", "ERFE")


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_table(path):
    entries = {}
    with gzip.open(path, "rt", encoding="utf8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not {"gene_id", "gene_symbol", "logFC", "logCPM", "FDR"}.issubset(reader.fieldnames or []):
            raise ValueError("MISSING_EDGER_GENE_FIELDS")
        for row in reader:
            gene = row["gene_id"].split(".")[0]
            if not gene.startswith("ENSG") or not gene[4:].isdigit() or len(gene) != 15:
                raise ValueError("BAD_GENE_ID")
            if gene in entries:
                raise ValueError("DUPLICATE_GENE")
            vals = {k: float(row[k]) for k in ("logFC", "logCPM", "FDR")}
            if not all(math.isfinite(v) for v in vals.values()) or not (0 <= vals["FDR"] <= 1):
                raise ValueError("NONFINITE_GENE_STAT_OR_FDR")
            vals["gene_symbol"] = row["gene_symbol"]
            entries[gene] = vals
    if not entries:
        raise ValueError("EMPTY_GENE_TABLE")
    return entries


def compare(a, b, key):
    shared = sorted(set(a) & set(b))
    if not shared:
        raise ValueError("NO_MATCHED_GENE_IDS")
    x = np.array([a[g]["logFC"] for g in shared], dtype=float)
    y = np.array([b[g]["logFC"] for g in shared], dtype=float)
    abundance = np.array([a[g]["logCPM"] for g in shared], dtype=float)
    fa = np.array([a[g]["FDR"] for g in shared], dtype=float)
    fb = np.array([b[g]["FDR"] for g in shared], dtype=float)
    signs = (x * y) > 0
    double_fdr = (fa < .05) & (fb < .05)
    absdelta = np.abs(x - y)
    corr = float(np.corrcoef(x, y)[0, 1])
    if not math.isfinite(corr):
        raise ValueError("UNDEFINED_GENE_CORRELATION")
    result = []
    for lo, hi, label in EFFECT_CUTS:
        mask = (np.abs(x) >= lo) & (np.abs(x) < hi)
        n = int(mask.sum())
        result.append({
            "comparison": key, "stratum_type": "baseline_abs_log2FC",
            "stratum": label, "n_shared_tested": n,
            "sign_concordant_n": int(signs[mask].sum()),
            "sign_concordance_percent": round(100 * signs[mask].mean(), 3) if n else None,
            "both_fdr_lt_0_05_n": int(double_fdr[mask].sum()),
            "median_absolute_delta_log2FC": round(float(np.median(absdelta[mask])), 6) if n else None,
            "baseline_bin_definition_not_independent": True,
        })
    # Explicitly derived from baseline logCPM distribution on the shared test universe.
    q = np.quantile(abundance, [0, .25, .5, .75, 1])
    for idx in range(4):
        mask = ((abundance >= q[idx]) &
                ((abundance <= q[idx + 1]) if idx == 3 else (abundance < q[idx + 1])))
        n = int(mask.sum())
        result.append({
            "comparison": key, "stratum_type": "baseline_logCPM_quartile",
            "stratum": f"Q{idx + 1}", "n_shared_tested": n,
            "sign_concordant_n": int(signs[mask].sum()),
            "sign_concordance_percent": round(100 * signs[mask].mean(), 3) if n else None,
            "both_fdr_lt_0_05_n": int(double_fdr[mask].sum()),
            "median_absolute_delta_log2FC": round(float(np.median(absdelta[mask])), 6) if n else None,
            "baseline_logCPM_lower": float(q[idx]),
            "baseline_logCPM_upper": float(q[idx + 1]),
            "baseline_bin_definition_not_independent": True,
        })
    markers = []
    for symbol in PANEL:
        matching = [g for g in shared if a[g]["gene_symbol"] == symbol and b[g]["gene_symbol"] == symbol]
        if len(matching) > 1:
            raise ValueError("MULTIPLE_ENSEMBL_IDS_FOR_PANEL_SYMBOL_" + symbol)
        if not matching:
            markers.append({"comparison": key, "gene_symbol": symbol, "status": "NOT_TESTED_IN_BOTH"})
        else:
            g = matching[0]
            markers.append({
                "comparison": key, "gene_symbol": symbol, "gene_id": g,
                "baseline_log2FC": a[g]["logFC"], "sensitivity_log2FC": b[g]["logFC"],
                "baseline_FDR": a[g]["FDR"], "sensitivity_FDR": b[g]["FDR"],
                "same_direction": (a[g]["logFC"] * b[g]["logFC"] > 0),
                "both_FDR_lt_0_05": (a[g]["FDR"] < .05 and b[g]["FDR"] < .05),
                "status": "BOTH_TESTED",
            })
    return {
        "summary": {
            "comparison": key, "baseline_tested": len(a), "sensitivity_tested": len(b),
            "shared_tested": len(shared), "same_direction_n": int(signs.sum()),
            "same_direction_fraction": round(float(signs.mean()), 8),
            "both_FDR_lt_0_05_n": int(double_fdr.sum()),
            "pearson_log2FC": round(corr, 8),
            "median_abs_delta_log2FC": round(float(np.median(absdelta)), 8),
            "no_independent_biological_replication": True,
            "library_treatment_aliasing_unresolved": True,
        },
        "bins": result,
        "panel": markers,
    }


def write_tsv(path, rows):
    fields = sorted(set().union(*(x.keys() for x in rows)))
    with path.open("w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def validate_frozen(a, b, c, d, r1, r2):
    if (len(a), len(b), len(c), len(d)) != (16050, 16728, 16051, 14003):
        raise ValueError("FROZEN_TEST_UNIVERSE_DRIFT")
    s1, s2 = r1["summary"], r2["summary"]
    if s1["shared_tested"] != 16050 or abs(s1["pearson_log2FC"] - 0.99911597) > .0002:
        raise ValueError("57DONOR_THRESHOLD_COMPARISON_SOURCE_DRIFT")
    if s2["shared_tested"] != 13975 or abs(s2["pearson_log2FC"] - .94990493) > .0002:
        raise ValueError("55DONOR_ANNOTATION_COMPARISON_SOURCE_DRIFT")
    if s2["same_direction_n"] != 12634 or s2["both_FDR_lt_0_05_n"] != 10255:
        raise ValueError("55DONOR_SIGN_OR_FDR_COUNT_DRIFT")
    high = lambda dic: {g for g, row in dic.items()
                        if row["FDR"] < .05 and abs(row["logFC"]) >= 3}
    if (len(high(a)), len(high(b)), len(high(a) & high(b))) != (83, 88, 81):
        raise ValueError("57DONOR_HIGH_EFFECT_INTERNAL_DRIFT")
    for result in (r1, r2):
        for gene in CORE:
            hit = next((r for r in result["panel"] if r["gene_symbol"] == gene), None)
            if not hit or hit["status"] != "BOTH_TESTED" or not (
                hit["same_direction"] and hit["both_FDR_lt_0_05"]):
                raise ValueError("CORE_THREE_GENE_SIGN_FDR_DRIFT_" + gene)
    erfe = next(r for r in r2["panel"] if r["gene_symbol"] == "ERFE")
    if erfe["status"] != "NOT_TESTED_IN_BOTH":
        raise ValueError("ERFE_MUST_BE_MARKED_UNTESTED_NOT_NULL_FDR")


def figure(rows, output):
    """Draw caption-safe plots: per-bin n values belong in dedicated tables.

    The archived 2026-11-02 initial plot placed n labels at y≈100 in both
    models and visibly overlapped at the top. New layout has no in-plot count
    labels; transparent source counts are fixed in tabular areas beneath plots.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.8))
    fig.subplots_adjust(left=.085, right=.985, top=.84, bottom=.36, wspace=.24)
    for ax, kind, xlabel, labels in [
        (axes[0], "baseline_abs_log2FC", "Baseline |log2 fold-change|",
         [x[2] for x in EFFECT_CUTS]),
        (axes[1], "baseline_logCPM_quartile", "Baseline logCPM quartile",
         [f"Q{i}" for i in range(1, 5)]),
    ]:
        counts=[]
        for key, name, fmt in [
            ("score", "57 paired: BeamA 0.75 vs 0.50", "o-"),
            ("label", "55 paired: BeamA vs TM_Fibro", "s--"),
        ]:
            data=[x for x in rows if x["comparison"] == key and x["stratum_type"] == kind]
            order={label: i for i, label in enumerate(labels)}
            data.sort(key=lambda x: order[x["stratum"]])
            if len(data)!=4: raise ValueError("EXPECTED_FOUR_STRATA")
            y=[d["sign_concordance_percent"] for d in data]
            if any(v is None for v in y): raise ValueError("STRATUM_EMPTY_SIGN")
            ax.plot(range(4), y, fmt, markersize=6, linewidth=2, label=name)
            counts.append([f"{d['n_shared_tested']:,}" for d in data])
        ax.set_ylim(0, 105)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_ylabel("Same-sign log2FC (%)")
        ax.set_xticks(range(4))
        ax.set_xticklabels(labels, fontsize=8, rotation=10 if kind=="baseline_abs_log2FC" else 0)
        ax.set_xlabel(xlabel, labelpad=9)
        ax.grid(axis="y", alpha=.2)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        # Dedicated n table: avoids annotation collisions when sign agreement
        # is 98–100% in two series. Each plotted point has an explicit n.
        table=ax.table(cellText=counts, rowLabels=["57 donor n", "55 donor n"],
                       cellLoc="center", rowLoc="center", bbox=[.015,-.59,.97,.27])
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        for (r,c),cell in table.get_celld().items():
            cell.set_linewidth(.35)
            if c==-1: cell.set_facecolor("#f2f4f7")
    axes[0].legend(frameon=False, loc="lower right", fontsize=8)
    fig.suptitle("GSE313736 — within-cohort sign sensitivity by baseline effect/expression",
                 y=.965, fontsize=12)
    fig.text(.5, -.035,
             "N values are shared tested genes in each baseline-defined bin. "
             "Both comparisons reuse original 10x libraries and are NOT independent biological replication.",
             ha="center", va="center", fontsize=9, color="#4b5563")
    fig.savefig(output, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)

def selftest():
    a = {}
    b = {}
    for i, (x, y, lcpm) in enumerate([
        (.1, -.1, 1), (.6, .5, 2), (1.5, 1.0, 3), (4.0, 3.0, 4),
        (.0, .0, 0),
    ], 1):
        gene = f"ENSG{i:011d}"
        a[gene] = {"logFC": x, "logCPM": lcpm, "FDR": .01, "gene_symbol": "GENE" + str(i)}
        b[gene] = {"logFC": y, "logCPM": lcpm, "FDR": .01, "gene_symbol": "GENE" + str(i)}
    o = compare(a, b, "selftest")
    assert o["summary"]["shared_tested"] == 5
    assert o["summary"]["same_direction_n"] == 3
    assert sum(d["n_shared_tested"] for d in o["bins"] if d["stratum_type"] == "baseline_abs_log2FC") == 5
    assert sum(d["n_shared_tested"] for d in o["bins"] if d["stratum_type"] == "baseline_logCPM_quartile") == 5
    with tempfile.TemporaryDirectory() as t:
        p = Path(t) / "test.tsv.gz"
        with gzip.open(p, "wt", encoding="utf8") as f:
            f.write("gene_id\tgene_symbol\tlogFC\tlogCPM\tFDR\n")
            f.write("ENSG00000000001\tMYOC\t4\t3\t0.001\n")
        assert len(read_table(p)) == 1
    print("OCULAR_STRATIFIED_SELFTEST_PASS", flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--base57", type=Path)
    p.add_argument("--score050", type=Path)
    p.add_argument("--beam55", type=Path)
    p.add_argument("--fibro55", type=Path)
    p.add_argument("--out", type=Path)
    p.add_argument("--self-test", action="store_true")
    a = p.parse_args()
    if a.self_test:
        selftest()
        return
    if not all((a.base57, a.score050, a.beam55, a.fibro55, a.out)):
        p.error("--base57, --score050, --beam55, --fibro55 and --out required")
    data = [read_table(f) for f in (a.base57, a.score050, a.beam55, a.fibro55)]
    r1 = compare(data[0], data[1], "score")
    r2 = compare(data[2], data[3], "label")
    validate_frozen(*data, r1, r2)
    bins = r1["bins"] + r2["bins"]
    panel = r1["panel"] + r2["panel"]
    a.out.mkdir(parents=True, exist_ok=True)
    write_tsv(a.out / "ocular_effect_stratum_sign_and_fdr.tsv", bins)
    write_tsv(a.out / "ocular_sensitivity_panel_effects.tsv", panel)
    figure(bins, a.out / "ocular_effect_stratified_sign_concordance.png")
    report = {
        "status": "PASS_REAL_AGGREGATE_STRATIFIED_SENSITIVITY",
        "dataset": "GSE313736",
        "source_ci_runs": {"score075_frozen": 41085982185,
                           "score050": 41199532259,
                           "dual55_original": 41206289766},
        "model_units": {"score": "same 57 paired biological donors",
                        "label": "same 55 paired biological donors"},
        "source_file_sha256": dict(zip(("score075", "score050", "beam55", "fibro55"),
                                       [sha256_file(f) for f in (
                                           a.base57, a.score050, a.beam55, a.fibro55)])),
        "models": [r1["summary"], r2["summary"]],
        "effect_strata": bins,
        "markers": panel,
        "original_authors_80_gene_identity": "HOLD_UNAVAILABLE",
        "confidence_or_cell_annotation_sensitivity_is_not_independent_replication": True,
        "response_is_confounded_with_original_source_10x_library": True,
        "p_values_not_validated_by_this_descriptive_stratification": True,
        "quantile_bins_data_derived_and_not_pre_registered": True,
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    (a.out / "ocular_stratified_sensitivity_manifest.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf8")
    md = [
        "# GSE313736: stratified actual edgeR gene-effect robustness",
        "Descriptive within-cohort analysis only; no independent biological replication or batch deconfounding.",
        "",
        "| Analysis | Shared tested | Pearson log2FC | Same sign | FDR<0.05 both |",
        "|---|---:|---:|---:|---:|",
    ]
    for r in (r1, r2):
        m = r["summary"]
        md.append(f"| {m['comparison']} | {m['shared_tested']:,} | {m['pearson_log2FC']:.5f} | "
                  f"{m['same_direction_n']:,} | {m['both_FDR_lt_0_05_n']:,} |")
    md.extend(["", "## Effect-size and abundance strata", "",
               "| Comparison | Stratum | n | Same sign (%) | FDR both n | Median |Δlog2FC| |",
               "|---|---|---:|---:|---:|---:|"])
    for z in bins:
        md.append(f"| {z['comparison']} | {z['stratum']} | {z['n_shared_tested']:,} | "
                  f"{z['sign_concordance_percent']} | {z['both_fdr_lt_0_05_n']:,} | "
                  f"{z['median_absolute_delta_log2FC']} |")
    md.extend(["", "Caution: all bins are defined using the **baseline** model; the grouping is post-hoc and not independently validated.",
               "Near-zero effects have intrinsically unstable signs, and Pearson correlation can hide this effect.",
               "The DEX and Vehicle arms originated in separate treatment-specific 10x libraries.",
               "No patient steroid ΔIOP, original authors' DEG80 IDs, protein activity or causal validation used."])
    (a.out / "ocular_stratified_sensitivity_report.md").write_text("\n".join(md) + "\n", encoding="utf8")
    print("OCULAR_REAL_STRATIFIED_PASS",
          json.dumps({"score": r1["summary"], "label": r2["summary"]}), flush=True)
    for row in bins:
        print("OCULAR_STRATUM", json.dumps({k: row[k] for k in (
            "comparison", "stratum", "n_shared_tested", "sign_concordance_percent",
            "both_fdr_lt_0_05_n")}), flush=True)


if __name__ == "__main__":
    main()
