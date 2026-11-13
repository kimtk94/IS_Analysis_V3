#!/usr/bin/env python3
"""Aggregate eight GSE331114 CellTypist pilot model predictions by GEO sample.

No labels have been independently validated. The sole purpose is to detect
cross-library instability and record whether reference-based annotation merits
separate validation. Cell barcodes never become independent biological n.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import statistics

MODEL="Immune_All_High.pkl"
EXPECTED={(group,i):f"GSM{9739251+i+(4 if group=='Control' else 0)}"
          for group in ("AIS","Control") for i in range(4)}
COARSE_LABELS=("Monocytes","Macrophages","T cells","B cells","DC","ILC","NK cells","Mono-mac")


def validate_and_collect(root):
    rows=[]
    common_sha=None
    common_gene_count=None
    for (group,index),gsm in EXPECTED.items():
        directory=root/f"is-celltypist-reference-pilot-{group}-{index}"
        path=directory/"reference_pilot.json"
        if not path.is_file():
            raise ValueError("Missing model evidence for 8/8 GSM: "+str(path))
        d=json.loads(path.read_text(encoding="utf-8"))
        if (d.get("GSM")!=gsm or d.get("group")!=group or
                d.get("group_index_0based")!=index or d.get("source")!="NCBI_GEO_GSE331114"
                or d.get("status")!="REFERENCE_PILOT_PROVISIONAL_COMPLETE"
                or d.get("model")!=MODEL or d.get("original_source_SHA256_match") is not True
                or d.get("reference_predictions_are_provisional") is not True
                or d.get("source_cell_labels_confirmed") is not False
                or d.get("brain_microglia_versus_peripheral_myeloid_validated") is not False
                or d.get("independent_donor_identity_confirmed") is not False
                or d.get("pseudobulk_DE")!="NOT_RUN"):
            raise ValueError("Data provenance or biological inference mismatch for "+gsm)
        if d.get("reference_genes_matched_fraction",0)<0.6:
            raise ValueError("Reference features too incomplete for "+gsm)
        sha=d.get("model_sha256")
        if not sha or len(sha)!=64:
            raise ValueError("No valid model hash for "+gsm)
        if common_sha is not None and sha!=common_sha:
            raise ValueError("Different model hashes across 8 samples")
        common_sha=sha
        overlap=d["reference_genes_matched"]
        if common_gene_count is not None and overlap!=common_gene_count:
            raise ValueError("Model feature overlap differs across samples")
        common_gene_count=overlap
        n=int(d["n_predictions"])
        if n!=int(d["qc_barcodes_passed"]) or n<100:
            raise ValueError("Model prediction count differs from QC for "+gsm)
        labels=d["provisional_labels_top_20"]
        if len({r["model_label"] for r in labels})!=len(labels):
            raise ValueError("Duplicate model class in top 20")
        counts={r["model_label"]:int(r["n_cells"]) for r in labels}
        total=sum(counts.values())
        if total>n or any(v<0 for v in counts.values()):
            raise ValueError("Impossible label count for "+gsm)
        row={
            "GSM":gsm, "group":group, "n_cells_QC_and_predicted":n,
            "original_barcodes":d["original_barcodes"],
            "model_gene_overlap_fraction":d["reference_genes_matched_fraction"],
            "top_class_score_median_uncalibrated":d["top_score_median"],
            "fraction_top_score_ge_0_5":d["fraction_score_ge_threshold"],
            "other_labels_not_in_top20_count":n-total,
        }
        for label in COARSE_LABELS:
            key=label.lower().replace(" ","_").replace("-","_")
            row[key+"_model_fraction"]=round(counts.get(label,0)/n,6)
        rows.append(row)
    if len(rows)!=8 or len({r["GSM"] for r in rows})!=8:
        raise ValueError("Duplicate or absent GEO sample")
    return rows,common_sha,common_gene_count


def summarize(rows):
    results=[]
    label_keys=[k for k in rows[0] if k.endswith("_model_fraction")]
    for key in label_keys:
        a=[r[key] for r in rows if r["group"]=="AIS"]
        c=[r[key] for r in rows if r["group"]=="Control"]
        if len(a)!=4 or len(c)!=4:
            raise ValueError("Must have four libraries per group")
        results.append({
            "model_label":key.replace("_model_fraction",""),
            "AIS_median_provisional_fraction":round(statistics.median(a),6),
            "Control_median_provisional_fraction":round(statistics.median(c),6),
            "AIS_range":f"{min(a):.3f}..{max(a):.3f}",
            "Control_range":f"{min(c):.3f}..{max(c):.3f}",
            "difference_of_group_medians":round(statistics.median(a)-statistics.median(c),6),
            "p_value":"NOT_TESTED",
            "validated_celltype":"NO",
        })
    return results


def tsv(path, rows):
    with open(path,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def run(source,out):
    samples,model_sha,genes=validate_and_collect(source)
    summary=summarize(samples)
    out.mkdir(parents=True,exist_ok=True)
    tsv(out/"per_GSM_reference_predictions.tsv",samples)
    tsv(out/"AIS_control_provisional_group_summaries.tsv",summary)
    evidence={
        "completed_at_utc":datetime.now(timezone.utc).isoformat(),
        "dataset":"GSE331114",
        "biological_samples":8,"AIS_sample_records":4,"Control_sample_records":4,
        "same_reference_model_used":True,"model":MODEL,
        "model_sha256":model_sha,
        "matched_reference_genes_each_sample":genes,
        "all_eight_original_files_verified_against_prior_QC":True,
        "scientific_status":"PROVISIONAL_REFERENCE_COHORT_QC_ONLY",
        "cell_types_validated":False,
        "microglia_specific_reference_validation":False,
        "independent_donors_confirmed":False,
        "celltype_pseudobulk_gene_counts_ready":False,
        "celltype_DE":"NOT_RUN",
        "eQTL_GWAS_causal_inference":"NOT_RUN",
        "warnings":[
            "High CellTypist scores are not calibrated cell-type accuracies",
            "Blood/immune reference may classify brain microglia as monocytes or macrophages",
            "CD45 enrichment, clinical sampling and batch differences can confound cell proportions",
            "Numerical results are sample-median model label fractions, not prevalence or differential abundance tests",
            "Do not treat 8 GSMs as verified independent donors until metadata is confirmed",
        ],
    }
    (out/"evidence.json").write_text(json.dumps(evidence,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    lines=[
        "# GSE331114 — 8 GEO sample CellTypist reference feasibility synthesis",
        "",
        "**Model labels are preliminary. No formal cell-type identity, donor-pseudobulk or DE statistics.**",
        "",
        "| Model provisional label | AIS median | Control median | Difference |",
        "|---|---:|---:|---:|",
    ]
    for x in summary:
        lines.append(f"| {x['model_label']} | "
                     f"{x['AIS_median_provisional_fraction']:.3f} | "
                     f"{x['Control_median_provisional_fraction']:.3f} | "
                     f"{x['difference_of_group_medians']:+.3f} |")
    lines += [
        "",
        f"Reference gene overlap {genes} per sample; identical model SHA-256 \`{model_sha}\`.",
        "",
        "**Decision: HOLD formal pseudobulk DEG pending independent brain immune/microglia label validation,",
        "doublet assessment, donor independence, clinical/batch metadata, and celltype gene-wide count inputs.**",
    ]
    md="\n".join(lines)+"\n"
    (out/"summary.md").write_text(md,encoding="utf-8")
    print(md)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"],"a",encoding="utf-8") as f:
            f.write(md)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",default="downloaded_reference_predictions")
    p.add_argument("--out-dir",default="artifacts/is_celltypist_8sample_synthesis")
    a=p.parse_args()
    run(Path(a.source),Path(a.out_dir))


if __name__=="__main__":
    main()
