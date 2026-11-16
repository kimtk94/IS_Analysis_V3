#!/usr/bin/env python3
"""Reference-based annotation FEASIBILITY pilot on GSE331114 public 10x MEX.

Runs the official CellTypist human Immune_All_High reference on a single
GEO GSM (one AIS or one Control) with library/mitochondrial QC and documents
source + model hashes and gene overlap. Reference labels are PROVISIONAL;
in particular, microglia vs recruited myeloid cells is not validated.

Writes only GSM-level aggregates, never cell barcodes, patient identifiers,
per-cell probabilities, or raw/restricted inputs into GitHub artifacts.
No pseudobulk differential expression, ancestry inference or eQTL.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import tempfile

from audit_stroke_versus_rosmap_sources import (
    download_metadata, geo_soft_url, parse_soft, derive_group,
)
from audit_stroke_gse331114_mex_pilot import (
    group_files, download, CAPS,
)
from sall1_microglia_orthogonal_check import ann_data_binary_markers

MODEL_NAME = "Immune_All_High.pkl"
MODEL_REFERENCE = "CellTypist human pan-tissue immune high-hierarchy, not stroke brain trained"
MIN_MODEL_OVERLAP_FRACTION = 0.60
PASS_N_GENES_MIN = 200
PASS_N_GENES_MAX = 7500
PASS_MT_PCT_LT = 20.0
MIN_PILOT_CELLS = 100
SCORE_CUTOFF_DESCRIPTIVE = 0.5


def select_gsm(group, index, samples):
    subset = sorted((x for x in samples.values() if derive_group(x) == group),
                    key=lambda x: x["accession"])
    if len(subset) != 4 or index not in range(4):
        raise ValueError("Expected four biological GEO libraries per group")
    return subset[index]


def summarize_predictions(labels, scores, expected_n, cutoff=SCORE_CUTOFF_DESCRIPTIVE):
    import numpy as np
    if len(labels) != expected_n or len(scores) != expected_n:
        raise ValueError("Prediction rows do not match QC passing cells")
    if expected_n < MIN_PILOT_CELLS:
        raise ValueError("Not enough QC passing cells for reference pilot")
    arr = np.asarray(scores, dtype=float)
    if not np.all(np.isfinite(arr)) or np.any(arr < 0) or np.any(arr > 1):
        raise ValueError("CellTypist score values outside [0,1]")
    counted = Counter(str(x) for x in labels)
    if not counted or sum(counted.values()) != expected_n:
        raise ValueError("Label count mismatch")
    resolved = sum(float(score) >= cutoff for score in scores)
    result = {
        "celltypist_raw_top_class_score_is_not_calibrated_probability": True,
        "score_threshold_descriptive_only": cutoff,
        "n_predictions": expected_n,
        "n_with_top_score_at_least_threshold": resolved,
        "fraction_score_ge_threshold": round(resolved/expected_n, 5),
        "top_score_median": round(float(np.median(arr)), 5),
        "top_score_min": round(float(np.min(arr)), 5),
        "top_score_max": round(float(np.max(arr)), 5),
        "provisional_labels_top_20": [
            {"model_label":label,"n_cells":n,
             "fraction":round(n/expected_n, 5)}
            for label,n in counted.most_common(20)
        ],
    }
    return result


def run(group, index, source_root, out):
    import numpy as np
    import scanpy as sc
    import celltypist
    from celltypist import models

    raw, _ = download_metadata(geo_soft_url("GSE331114"))
    info, samples = parse_soft(raw)
    if info.get("accession") != "GSE331114":
        raise ValueError("GEO source mismatch")
    sample = select_gsm(group, index, samples)
    name = sample["accession"]
    source_qc = source_root/f"is-gse331114-{group}-{index}"/"gene_counts.json"
    prior = json.loads(source_qc.read_text(encoding="utf-8"))
    if (prior.get("accession") != name or prior.get("group") != group or
            prior.get("status") != "SAMPLE_MEX_GENE_QC_PASS"):
        raise ValueError("Prior original source GSM or QC mismatch")
    files = group_files(sample)
    if set(files) != {"matrix","features","barcodes"}:
        raise ValueError("Missing public 10x MEX triplet")
    with tempfile.TemporaryDirectory(prefix="is_reference_pilot_") as temp:
        root=Path(temp)
        downloaded={}
        for component in ("features","barcodes","matrix"):
            dest=root/f"{component}.tsv.gz" if component != "matrix" else root/"matrix.mtx.gz"
            digest=download(files[component], dest, CAPS[component])["sha256_local"]
            if digest != prior.get(component+"_sha256_not_source_verified"):
                raise ValueError("Public source SHA256 differs from previous verified run")
            downloaded[component]=dest
        adata=sc.read_10x_mtx(root, var_names="gene_symbols",make_unique=True,gex_only=True)
        if adata.n_obs != prior["cell_barcodes"] or adata.n_vars != prior["genes"]:
            raise ValueError("MEX cell / gene dimensions changed")
        total_by_cell=np.asarray(adata.X.astype(np.float64).sum(axis=1)).ravel()
        if abs(float(total_by_cell.sum())-float(prior["total_umi_count"])) > 0.1:
            raise ValueError("MEX total UMI differs from source QC")
        n_genes=np.asarray((adata.X>0).sum(axis=1)).ravel()
        mt=adata.var_names.str.startswith("MT-")
        if int(mt.sum())<10:
            raise ValueError("No mitochondrial features in source")
        mt_totals=np.asarray(adata[:,mt].X.sum(axis=1)).ravel()
        mt_pct=100.0*mt_totals/np.maximum(total_by_cell,1)
        mask=((n_genes>=PASS_N_GENES_MIN)&(n_genes<=PASS_N_GENES_MAX)&
              (mt_pct<PASS_MT_PCT_LT)&(total_by_cell>0))
        n_pass=int(mask.sum())
        if n_pass < MIN_PILOT_CELLS:
            raise ValueError("Too few cells pass simple QC")
        adata=adata[mask,:].copy()
        adata.X=adata.X.astype(np.float32)
        # CellTypist docs expect ln(1 + normalized 10,000 counts per cell).
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        models.download_models(model=MODEL_NAME)
        model=models.Model.load(model=MODEL_NAME)
        model_path=Path(models.get_model_path(MODEL_NAME))
        model_sha256=hashlib.sha256(model_path.read_bytes()).hexdigest()
        ref_features=set(str(x) for x in model.features)
        overlap=len(ref_features.intersection(set(adata.var_names)))
        fraction=overlap/max(1,len(ref_features))
        if fraction<MIN_MODEL_OVERLAP_FRACTION:
            raise ValueError("Reference gene overlap too low: "+str(round(fraction,4)))
        predictions=celltypist.annotate(adata,model=model,majority_voting=False,mode="best match")
        labels=predictions.predicted_labels["predicted_labels"].astype(str).tolist()
        # Sigmoid per-class scores are NOT independently validated probabilities.
        probs=predictions.probability_matrix
        best_scores=probs.to_numpy().max(axis=1).tolist()
        summary=summarize_predictions(labels,best_scores,n_pass)
        sall1_evidence=ann_data_binary_markers(adata,labels)

        evidence={
            "source":"NCBI_GEO_GSE331114",
            "GSM":name,"group":group,"group_index_0based":index,
            "original_source_run_id":43010007530,
            "original_source_SHA256_match":True,
            "original_barcodes":int(prior["cell_barcodes"]),
            "qc_barcodes_passed":n_pass,
            "qc_pass_fraction":round(n_pass/int(prior["cell_barcodes"]),5),
            "QC_limits":{"n_genes_min":PASS_N_GENES_MIN,
                         "n_genes_max":PASS_N_GENES_MAX,"mt_percent_lt":PASS_MT_PCT_LT},
            "model":MODEL_NAME,
            "model_reference":MODEL_REFERENCE,
            "model_sha256":model_sha256,
            "celltypist_version":celltypist.__version__,
            "scanpy_version":sc.__version__,
            "reference_gene_count":len(ref_features),
            "reference_genes_matched":overlap,
            "reference_genes_matched_fraction":round(fraction,5),
            "source_cell_labels_confirmed":False,
            "brain_microglia_versus_peripheral_myeloid_validated":False,
            "reference_predictions_are_provisional":True,
            "independent_donor_identity_confirmed":False,
            "pseudobulk_DE":"NOT_RUN",
            "genotype_eQTL":"NOT_AVAILABLE",
            "status":"REFERENCE_PILOT_PROVISIONAL_COMPLETE",
            "SALL1_per_GSM_orthogonal_RNA_markers":sall1_evidence,
            "clinical_covariate_values_collected":False,
            "GEO_sample_characteristic_names_only":sorted({
                x.split(':',1)[0].strip().lower()
                for x in sample.get('characteristics',[]) if ':' in x
            }),
            **summary,
        }
    out.mkdir(parents=True,exist_ok=True)
    (out/"reference_pilot.json").write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    lines=[f"# GSE331114 {name} {group}: CellTypist reference pilot","",
           f"**{evidence['status']}**, {n_pass}/{prior['cell_barcodes']} QC barcodes annotated (provisional only).",
           "",
           f"- Reference: {MODEL_NAME}, SHA-256: \`{model_sha256}\`",
           f"- Matched reference model genes: {overlap}/{len(ref_features)} ({fraction:.1%})",
           f"- Median top-class sigmoid score: {summary['top_score_median']:.3f} (NOT calibrated accuracy)",
           f"- Fraction top score ≥{SCORE_CUTOFF_DESCRIPTIVE}: {summary['fraction_score_ge_threshold']:.3f}",
           "",
           "| Provisional reference label | QC barcodes | Fraction |","|---|---:|---:|"]
    for row in evidence["provisional_labels_top_20"]:
        lines.append(f"| {row['model_label']} | {row['n_cells']} | {row['fraction']:.3f} |")
    lines += ["", "### SALL1 orthogonal marker RNA co-detection (NOT cell-type validation)", "",
              "| Model provisional label | n | SALL1 detected | SALL1 + ≥1 canonical | Strict SALL1/TMEM119/P2RY12 | ≥2 peripheral myeloid |",
              "|---|---:|---:|---:|---:|---:|"]
    for d in sall1_evidence["model_categories"]:
        if d["counts"] is None:
            continue
        c=d["counts"]
        lines.append(f"| {d['model_provisional_label']} | {d['n_barcodes']} | "
                     f"{c['SALL1_any']} | {c['SALL1_plus_one_canonical']} | "
                     f"{c['SALL1_TMEM119_P2RY12_triple']} | "
                     f"{c['at_least_two_FCN1_CD14_S100A8']} |")
    lines += ["", "These are binary RNA detection counts, NOT the published paper's SALL1-high classifier.",
              "Absence of detectable SALL1 can reflect single-cell RNA dropout, and a model label",
              "does NOT establish resident microglia or infiltrating monocyte lineage."]
    lines+=["","No cell-level predictions exported; donor reference mapping/validation, doublets and microglia-specific annotation remain HOLD.",
            "This pilot does NOT permit cell-type differential expression, independent clinical replication or GWAS mediation claims."]
    md="\n".join(lines)+"\n"
    (out/"summary.md").write_text(md,encoding="utf-8")
    print(md)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"],"a",encoding="utf-8") as f:
            f.write(md)
    return evidence


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--group",required=True,choices=("AIS","Control"))
    p.add_argument("--index",required=True,type=int,choices=range(4))
    p.add_argument("--source-root",default="downloaded_gse331114")
    p.add_argument("--out-dir",default="artifacts/is_celltypist_pilot")
    args=p.parse_args()
    run(args.group,args.index,Path(args.source_root),Path(args.out_dir))


if __name__=="__main__":
    main()
