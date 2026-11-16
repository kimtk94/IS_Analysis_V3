#!/usr/bin/env python3
"""Reconcile the eight GSE331114 QC artifacts and report donor-level descriptive CPM.

Inputs: exactly eight already completed GitHub Actions artifacts from run
43010007530. Outputs: sample-wise normalized expression and gene-level
descriptive summaries. NO inferential DEG, eQTL, patient risk prediction,
single-cell-as-replicate tests, or adjusted cell-type comparison.

The 10x MEX samples are CD45-enriched human ischemic stroke brain cells.
Total-library normalization across all cells cannot control composition/batch.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from html import escape
import json
import math
from pathlib import Path
import statistics
import sys

SOURCE_RUN_ID = 43010007530
EXPECTED = {
    ("AIS", 0): "GSM9739251",
    ("AIS", 1): "GSM9739252",
    ("AIS", 2): "GSM9739253",
    ("AIS", 3): "GSM9739254",
    ("Control", 0): "GSM9739255",
    ("Control", 1): "GSM9739256",
    ("Control", 2): "GSM9739257",
    ("Control", 3): "GSM9739258",
}
IS_CANDIDATES = (
    "ABO", "ALDH2", "CD40", "F11", "FGF5", "FURIN", "GRK5",
    "IL6R", "KIAA0319", "MMP12", "PROCR", "SCARA5",
    "SWAP70", "TFPI", "TMPRSS5",
)
CELL_CONTEXT = (
    "PTPRC", "AIF1", "TYROBP", "CD3D", "CD14", "MS4A1",
    "CX3CR1", "TREM2", "PECAM1", "GFAP", "RBFOX3",
)
TARGETS = IS_CANDIDATES + CELL_CONTEXT


def write_table(path, rows, headers):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def finite_nonnegative(value, label):
    try:
        x = float(value)
    except (ValueError, TypeError) as ex:
        raise ValueError(f"{label} is not numeric") from ex
    if not math.isfinite(x) or x < 0:
        raise ValueError(f"{label} must be finite and non-negative")
    return x


def reconcile(root: Path):
    """Validate all sample/job identities and fields before any comparison."""
    samples = []
    gene_rows = []
    evidence = []
    for (group, index), gsm in EXPECTED.items():
        artifact_name = f"is-gse331114-{group}-{index}"
        folder = root / artifact_name
        p = folder / "gene_counts.json"
        if not p.is_file():
            raise ValueError(f"Missing expected CI artifact {p}")
        record = json.loads(p.read_text(encoding="utf-8"))
        if (record.get("status") != "SAMPLE_MEX_GENE_QC_PASS"
                or record.get("group") != group
                or record.get("group_position_0based") != index
                or record.get("accession") != gsm):
            raise ValueError(f"Sample provenance/GEO identity or QC mismatch: {artifact_name}")
        if record.get("matched_genotype") is not False or record.get("clinical_gwas_patient_matching") is not False:
            raise ValueError("Invalid claims of genotype matching")
        if record.get("donor_unique_identity") != "NOT_INDEPENDENTLY_VERIFIED":
            raise ValueError("Premature donor-independence claim")
        if record.get("genes") != 36601:
            raise ValueError(f"Feature universe mismatch: {gsm}")
        barcodes = int(record["cell_barcodes"])
        count_total = finite_nonnegative(record["total_umi_count"], f"{gsm} total UMI")
        if not (100 <= barcodes <= 500000) or count_total <= 0:
            raise ValueError(f"Implausible barcodes/UMI library {gsm}")
        nrows = record.get("target_gene_rows", [])
        if len(nrows) != len(TARGETS):
            raise ValueError(f"Unexpected target panel size in {gsm}")
        by_gene = {}
        for r in nrows:
            gene = r.get("gene_symbol")
            if gene not in TARGETS or gene in by_gene:
                raise ValueError(f"Unexpected/duplicate gene {gene} in {gsm}")
            if r.get("sample_GSM") != gsm or r.get("sample_group") != group:
                raise ValueError(f"Gene-level GSM/group metadata mismatch for {gsm}")
            by_gene[gene] = r
        if set(by_gene) != set(TARGETS):
            raise ValueError(f"Missing genes for {gsm}")
        sample = {
            "sample_GSM": gsm,
            "group": group,
            "group_position_0based": index,
            "cell_barcodes": barcodes,
            "genes": record["genes"],
            "total_umis": round(count_total, 3),
            "umis_per_cell": round(count_total / barcodes, 2),
        }
        samples.append(sample)
        source_records = []
        for gene in TARGETS:
            item = by_gene[gene]
            counts = finite_nonnegative(item.get("umis_sum"), f"{gsm} {gene} UMI")
            detection = finite_nonnegative(item.get("fraction_cells_nonzero"), f"{gsm} {gene} detection")
            if detection > 1 or counts > count_total:
                raise ValueError(f"Invalid count or detection: {gsm} {gene}")
            present = item.get("gene_present_in_features")
            if present is not True and present is not False:
                raise ValueError("Feature presence status must be boolean")
            if counts > 0 and present is not True:
                raise ValueError("Count observed but gene absent from feature table")
            cpm = counts * 1000000.0 / count_total
            row = {
                "sample_GSM": gsm, "group": group,
                "gene_symbol": gene,
                "is_prespecified_IS_candidate": gene in IS_CANDIDATES,
                "gene_present_in_features": present,
                "umis_sum": int(counts) if counts.is_integer() else counts,
                "library_total_umis": round(count_total, 3),
                "normalized_CPM": round(cpm, 8),
                "log2_CPM_plus_1": round(math.log2(cpm + 1.0), 8),
                "fraction_cells_nonzero": round(detection, 8),
                "cell_barcodes": barcodes,
                "normalization_note": "All-cell library size only; not cell-type adjusted",
            }
            gene_rows.append(row)
            source_records.append(row)
        evidence.append({"artifact": artifact_name, "GSM": gsm, "source":"GSE331114",
                         "source_run_id": SOURCE_RUN_ID,
                         "per_sample_QC": record["status"], "genes_present": len(by_gene)})
    if len(samples) != 8 or len(gene_rows) != 8 * len(TARGETS):
        raise ValueError("Not all complete sample data were integrated")
    if len(set(x["sample_GSM"] for x in samples)) != 8:
        raise ValueError("Duplicate GSM across samples")
    return samples, gene_rows, evidence


def summarize_gene_rows(gene_rows):
    # Pairing is not presumed; only a small-n sample-wise descriptive comparison.
    result = []
    by = {}
    for r in gene_rows:
        by[(r["gene_symbol"], r["group"], r["sample_GSM"])] = r
    for gene in TARGETS:
        ais = [r for r in gene_rows if r["gene_symbol"] == gene and r["group"] == "AIS"]
        ctrl = [r for r in gene_rows if r["gene_symbol"] == gene and r["group"] == "Control"]
        if len(ais) != 4 or len(ctrl) != 4:
            raise ValueError(f"Requires 4+4 complete samples for {gene}")
        if not all(r["gene_present_in_features"] for r in ais + ctrl):
            # Missing features should NOT be converted into true zero expression.
            status = "FEATURE_MISSING_REVIEW"
        elif any(r["umis_sum"] == 0 for r in ais + ctrl):
            status = "ZERO_COUNT_SENSITIVITY_REQUIRED"
        else:
            status = "DESCRIPTIVE_ONLY"
        ais_cpm = [r["normalized_CPM"] for r in ais]
        ctrl_cpm = [r["normalized_CPM"] for r in ctrl]
        ais_log = [r["log2_CPM_plus_1"] for r in ais]
        ctrl_log = [r["log2_CPM_plus_1"] for r in ctrl]
        ais_detection = [r["fraction_cells_nonzero"] for r in ais]
        ctrl_detection = [r["fraction_cells_nonzero"] for r in ctrl]
        # Sample-level median log CPM difference; NOT a fitted DE effect.
        delta = statistics.median(ais_log) - statistics.median(ctrl_log)
        result.append({
            "gene_symbol": gene,
            "candidate_or_marker": "IS_candidate" if gene in IS_CANDIDATES else "cell_context",
            "n_AIS_GSM": 4, "n_Control_GSM": 4,
            "AIS_median_CPM": round(statistics.median(ais_cpm), 6),
            "Control_median_CPM": round(statistics.median(ctrl_cpm), 6),
            "AIS_median_log2_CPM_plus_1": round(statistics.median(ais_log), 6),
            "Control_median_log2_CPM_plus_1": round(statistics.median(ctrl_log), 6),
            "descriptive_delta_median_log2_CPM_plus_1": round(delta, 6),
            "AIS_min_log2_CPM_plus_1": round(min(ais_log), 6),
            "AIS_max_log2_CPM_plus_1": round(max(ais_log), 6),
            "Control_min_log2_CPM_plus_1": round(min(ctrl_log), 6),
            "Control_max_log2_CPM_plus_1": round(max(ctrl_log), 6),
            "AIS_median_detection_fraction": round(statistics.median(ais_detection), 6),
            "Control_median_detection_fraction": round(statistics.median(ctrl_detection), 6),
            "status": status,
            "p_value": "NOT_TESTED",
            "fdr": "NOT_TESTED",
            "celltype_adjustment": "NOT_PERFORMED",
            "donor_independence": "NOT_INDEPENDENTLY_VERIFIED",
        })
    return result


def svg_plot(path, summary):
    """Simple descriptive dot plot (all GSM sample summaries are in TSV)."""
    rows = [r for r in summary if r["candidate_or_marker"] == "IS_candidate"]
    width, height = 1300, 72 * len(rows) + 130
    left, right = 180, 115
    x0, x1 = left, width - right
    vals = [r[k] for r in rows for k in
            ("AIS_median_log2_CPM_plus_1", "Control_median_log2_CPM_plus_1")]
    xmax = max(1, math.ceil(max(vals) * 1.1))
    def xp(x):
        return x0 + (x1-x0)*(x / xmax)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#222} .t{font-size:18px;font-weight:600} .s{font-size:13px} .g{font-size:15px}</style>',
        '<text x="25" y="28" class="t">GSE331114 — descriptive median log2(CPM + 1)</text>',
        '<text x="25" y="48" class="s">AIS n=4; non-stroke n=4. CD45+ brain immune samples; NO DEG test or cell-type adjustment.</text>',
        '<circle cx="28" cy="75" r="6" fill="#111"/><text x="42" y="80" class="s">AIS</text>',
        '<circle cx="115" cy="75" r="6" fill="#999"/><text x="129" y="80" class="s">Control</text>',
    ]
    for index, r in enumerate(rows):
        y = 116 + index * 72
        if index % 2 == 0:
            parts.append(f'<rect x="10" y="{y-27}" width="{width-20}" height="66" fill="#f5f5f5"/>')
        for tick in range(xmax+1):
            x = xp(tick)
            if tick % max(1, math.ceil(xmax/10)) == 0:
                parts.append(f'<line x1="{x:.1f}" y1="{y-20}" x2="{x:.1f}" y2="{y+31}" stroke="#ddd"/>')
        xa = xp(r["AIS_median_log2_CPM_plus_1"])
        xc = xp(r["Control_median_log2_CPM_plus_1"])
        parts.extend([
            f'<text class="g" x="18" y="{y+6}">{escape(r["gene_symbol"])}</text>',
            f'<line x1="{min(xa,xc):.1f}" y1="{y+1}" x2="{max(xa,xc):.1f}" y2="{y+1}" stroke="#777" stroke-width="2"/>',
            f'<circle cx="{xc:.1f}" cy="{y+1}" r="7" fill="#999"/>',
            f'<circle cx="{xa:.1f}" cy="{y+1}" r="7" fill="#111"/>',
        ])
    parts.append(f'<text class="s" x="{x0}" y="{height-18}">Median log2(CPM + 1) per sample (library-size normalized; cell-composition confounding not corrected)</text>')
    parts.append('</svg>')
    path.write_text("\n".join(parts), encoding="utf-8")


def run(root, out):
    samples, genes, evidence = reconcile(root)
    results = summarize_gene_rows(genes)
    out.mkdir(parents=True, exist_ok=True)
    write_table(out / "samples_qc.tsv", samples,
                ["sample_GSM","group","group_position_0based","cell_barcodes","genes","total_umis","umis_per_cell"])
    write_table(out / "per_sample_candidate_genes_cpm.tsv", genes,
                ["sample_GSM","group","gene_symbol","is_prespecified_IS_candidate",
                 "gene_present_in_features","umis_sum","library_total_umis","normalized_CPM",
                 "log2_CPM_plus_1","fraction_cells_nonzero","cell_barcodes","normalization_note"])
    write_table(out / "gene_group_descriptive.tsv", results,
                list(results[0].keys()))
    svg_plot(out / "descriptive_gene_medians.svg", results)
    audit = {
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "source_run": SOURCE_RUN_ID,
        "source_run_url": f"https://github.com/kimtk94/IS_Analysis_V3/actions/runs/{SOURCE_RUN_ID}",
        "source_evidence": evidence,
        "verified_8of8_artifact_identity_and_source_qc": True,
        "n_GSM": len(samples), "n_AIS": 4, "n_control": 4,
        "n_tested_gene_targets": len(TARGETS),
        "n_target_rows": len(genes),
        "statistical_testing": "NOT_PERFORMED",
        "pseudobulk_cell_type": "NOT_PERFORMED",
        "donor_independence": "NOT_INDEPENDENTLY_VERIFIED",
        "interpretation": "Bulk across all CD45+ cells per GSM, library CPM, descriptive medians only. Clinical groups not causally exchangeable.",
        "genotype_matched": False,
        "ROSMAP_eQTL_substituted": False,
        "assay_data": "human acute ischemic stroke brain CD45-enriched expression; NOT eQTL",
        "status": "EIGHT_SAMPLE_DESCRIPTIVE_AGGREGATION_PASS",
    }
    (out / "evidence.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    md = [
        "# GSE331114 AIS/control: complete 8-sample descriptive gene summary",
        "",
        "**8/8 original 10x sample QC artifacts independently reconciled; no P values and no DEG claims.**",
        "",
        "| Gene | AIS median CPM | Control median CPM | Delta median log2(CPM+1) | QC status |",
        "|---|---:|---:|---:|---|",
    ]
    for r in results:
        if r["candidate_or_marker"] != "IS_candidate":
            continue
        md.append(f"| {r['gene_symbol']} | {r['AIS_median_CPM']:.3f} | {r['Control_median_CPM']:.3f} | "
                  f"{r['descriptive_delta_median_log2_CPM_plus_1']:+.3f} | {r['status']} |")
    md += [
        "",
        "### Interpretation limit",
        "- Sample median changes are **descriptive**. Small n=4+4, donor uniqueness and clinical matching not independently verified.",
        "- Sample CPM denominators are total raw UMI across CD45+ enriched cells; cell-type proportions may shift between AIS and control.",
        "- Existing preselected genes were chosen for historical IS research, not genome-wide prioritized from unbiased new GWAS loci.",
        "- A change after stroke is not proof of inherited ischemic stroke susceptibility, colocalization, or causal mediation.",
        "- Outputs include sample table, gene medians, raw per-sample public aggregate counts, evidence JSON and SVG plot.",
    ]
    summary = "\n".join(md)+"\n"
    (out / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)
    if path := __import__("os").environ.get("GITHUB_STEP_SUMMARY"):
        with open(path, "a", encoding="utf-8") as f:
            f.write(summary)
    return 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", default="downloaded_gse331114")
    p.add_argument("--out-dir", default="artifacts/is_gse331114_descriptive")
    args = p.parse_args()
    try:
        return run(Path(args.input), Path(args.out_dir))
    except (OSError, ValueError, KeyError, TypeError) as e:
        print(f"FAIL_CLOSED: {type(e).__name__}: {e}", file=sys.stderr)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
