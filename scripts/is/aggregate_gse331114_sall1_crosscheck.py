#!/usr/bin/env python3
"""Reconcile real 8-sample SALL1 marker crosschecks from one CellTypist CI run.

The sample is the unit of descriptive comparison, NOT cells as independent
replicates. NO diagnosis, genotype, donor-specific lineage, DE, or formal
cell-abundance inference. Every source marker result must be independent
provenance-checked against original GEO identifiers and identical model SHA.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import statistics
import sys

EXPECTED={(g,i):f"GSM{9739251+i+(4 if g=='Control' else 0)}"
          for g in ("AIS","Control") for i in range(4)}
MARKER_FIELDS=(
    "SALL1_any",
    "SALL1_plus_one_canonical",
    "SALL1_TMEM119_P2RY12_triple",
    "TMEM119_P2RY12_pair_regardless_SALL1",
    "at_least_two_FCN1_CD14_S100A8",
    "both_SALL1_canonical_and_peripheral",
)
KEY_LABELS=("Monocytes","Macrophages","Mono-mac","ALL_PREDICTED_MY_MARKERS")


def collect(root):
    data=[]
    refhash=None
    for (group,index),gsm in EXPECTED.items():
        path=root/f"is-celltypist-reference-pilot-{group}-{index}"/"reference_pilot.json"
        if not path.is_file():
            raise ValueError("Missing source CellTypist/SALL1 GEO evidence: "+str(path))
        source=json.loads(path.read_text(encoding="utf-8"))
        if (source.get("GSM")!=gsm or source.get("group")!=group
                or source.get("group_index_0based")!=index
                or source.get("status")!="REFERENCE_PILOT_PROVISIONAL_COMPLETE"
                or source.get("source")!="NCBI_GEO_GSE331114"
                or source.get("original_source_SHA256_match") is not True
                or source.get("source_cell_labels_confirmed") is not False
                or source.get("independent_donor_identity_confirmed") is not False
                or source.get("pseudobulk_DE")!="NOT_RUN"):
            raise ValueError("Sample lineage/provenance failure: "+gsm)
        sh=source.get("model_sha256")
        if not sh or len(sh)!=64 or (refhash is not None and sh!=refhash):
            raise ValueError("Reference model is not identical or hash absent: "+gsm)
        refhash=sh
        orth=source.get("SALL1_per_GSM_orthogonal_RNA_markers",{})
        if (orth.get("scope")!="PUBLIC_GEO_PER_GSM_ORTHOGONAL_BINARY_RNA_DETECTION"
                or orth.get("paper_SALL1_high_classifier_reproduced") is not False
                or orth.get("celltypes_validated") is not False
                or orth.get("myeloid_resident_versus_infiltrating_resolved") is not False):
            raise ValueError("No valid orthogonal SALL1 source or overclaimed celltype validity: "+gsm)
        if orth.get("n_total_barcodes")!=source.get("qc_barcodes_passed"):
            raise ValueError("SALL1/CellTypist barcode denominators disagree: "+gsm)
        groups={r.get("model_provisional_label"):r for r in orth.get("model_categories",[])}
        if len(groups)!=len(orth.get("model_categories",[])) or any(x not in groups for x in KEY_LABELS):
            raise ValueError("Missing/duplicate model label crosscheck: "+gsm)
        for label in KEY_LABELS:
            item=groups[label]
            n=int(item["n_barcodes"])
            counts=item.get("counts")
            rates=item.get("rates")
            if n<0 or n>source["qc_barcodes_passed"]:
                raise ValueError("Invalid selected label count")
            if n>0:
                if not isinstance(counts,dict) or not isinstance(rates,dict):
                    raise ValueError("No counts/rates for model label: "+label)
                if set(counts)!=set(MARKER_FIELDS) or set(rates)!=set(MARKER_FIELDS):
                    raise ValueError("Wrong six-marker signature universe")
                for key in MARKER_FIELDS:
                    num=counts[key]
                    if not isinstance(num,int) or not 0<=num<=n:
                        raise ValueError("Marker count out of range")
                    if abs(rates[key]-num/n)>1.5e-6:
                        raise ValueError("Marker rate/count mismatch")
            else:
                if counts is not None or rates is not None:
                    raise ValueError("Zero-cells should have no marker fractions")
            row={"GSM":gsm,"group":group,"group_index":index,
                 "provisional_reference_label":label,
                 "qc_cells_all":source["qc_barcodes_passed"],
                 "model_label_cells":n,
                 "model_hash":sh,
                 "SALL1_threshold":"detectable_RNA_only_NOT_published_high"}
            for key in MARKER_FIELDS:
                row[key+"_count"]=counts[key] if counts is not None else ""
                row[key+"_fraction"]=rates[key] if rates is not None else ""
            data.append(row)
        a=groups["ALL_PREDICTED_MY_MARKERS"]["n_barcodes"]
        b=sum(groups[label]["n_barcodes"] for label in ("Monocytes","Macrophages","Mono-mac"))
        if a!=b:
            raise ValueError("Combined provisional myeloid counts inconsistent")
    if len(data)!=8*len(KEY_LABELS):
        raise ValueError("Incomplete 8 sample by model label matrix")
    return data,refhash


def descriptive(data):
    out=[]
    for label in KEY_LABELS:
        for marker in MARKER_FIELDS:
            ai=[r[marker+"_fraction"] for r in data
                if r["group"]=="AIS" and r["provisional_reference_label"]==label
                and r["model_label_cells"]>0]
            ctrl=[r[marker+"_fraction"] for r in data
                  if r["group"]=="Control" and r["provisional_reference_label"]==label
                  and r["model_label_cells"]>0]
            status="DESCRIPTIVE_ONLY" if len(ai)==4 and len(ctrl)==4 else "INCOMPLETE_LABEL_SUPPORT"
            row={
                "reference_label_is_provisional":label,
                "binary_RNA_marker_signature":marker,
                "n_AIS_GSM":len(ai),
                "n_Control_GSM":len(ctrl),
                "AIS_sample_median_fraction":round(statistics.median(ai),6) if ai else "",
                "Control_sample_median_fraction":round(statistics.median(ctrl),6) if ctrl else "",
                "delta_sample_median":round(statistics.median(ai)-statistics.median(ctrl),6) if ai and ctrl else "",
                "status":status,
                "p_value":"NOT_TESTED",
                "FDR":"NOT_TESTED",
                "microglia_identity":"NOT_VALIDATED",
            }
            out.append(row)
    return out


def tsv(path,rows):
    with path.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run(root,out):
    data,refhash=collect(root)
    summary=descriptive(data)
    out.mkdir(parents=True,exist_ok=True)
    tsv(out/"GSM_label_SALL1_marker_detection.tsv",data)
    tsv(out/"descriptive_SALL1_group_medians.tsv",summary)
    allmy=[r for r in summary if r["reference_label_is_provisional"]=="ALL_PREDICTED_MY_MARKERS"]
    report={
        "finished_utc":datetime.now(timezone.utc).isoformat(),
        "source":"GSE331114",
        "all_eight_original_GEO_libraries_included":True,
        "n_AIS_GSM":4,"n_Control_GSM":4,
        "model_sha256":refhash,
        "original_10x_source_hash_reverified_each_sample":True,
        "actual_cell_type_validation":"HOLD",
        "microglia_lineage_distinguished":"HOLD",
        "pseudobulk_gene_counts":"NOT_CREATED",
        "pseudobulk_DEG":"NOT_RUN",
        "clinical_covariates_or_donor_independence_verified":False,
        "status":"SALL1_8SAMPLE_DESCRIPTIVE_CROSSCHECK_PASS",
        "scope":"Binary marker detection in provisional reference myeloid calls, not paper threshold or independent benchmark",
        "striking_limit":"SALL1 RNA absence may reflect transcript dropout; specificity/sensitivity not measured in this cohort",
    }
    (out/"evidence.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    lines=[
        "# Human IS CD45+ immune cells: 8-sample SALL1 cross-check",
        "",
        "**All 8 original GEO sample results reconciled. NOT cell-type validation and NOT published SALL1-high classifier.**",
        "",
        "| Predicted myeloid (n=4 GSM/group) marker evidence | AIS median | Control median | Delta |",
        "|---|---:|---:|---:|",
    ]
    for r in allmy:
        a=r["AIS_sample_median_fraction"]
        c=r["Control_sample_median_fraction"]
        d=r["delta_sample_median"]
        lines.append(f"| {r['binary_RNA_marker_signature']} | "
                     f"{a:.4f} | {c:.4f} | {d:+.4f} |"
                     if r["status"]=="DESCRIPTIVE_ONLY" else
                     f"| {r['binary_RNA_marker_signature']} | INCOMPLETE | INCOMPLETE | HOLD |")
    lines.extend(["","SALL1 positivity alone is not proven human microglia. "
                  "No statistical tests or donor-specific genotyping were performed.",
                  "Cell-type pseudobulk DEG remains HOLD until reference validation, donor metadata and batch checks."])
    md="\n".join(lines)+"\n"
    (out/"summary.md").write_text(md,encoding="utf-8")
    print(md)
    import os
    if x:=os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(x,"a",encoding="utf-8") as f:
            f.write(md)
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",default="downloaded_reference_predictions")
    p.add_argument("--out-dir",default="artifacts/is_SALL1_orthogonal")
    args=p.parse_args()
    run(Path(args.source),Path(args.out_dir))


if __name__=="__main__":
    try:
        main()
    except (ValueError,OSError,KeyError,TypeError) as exc:
        print("FAIL_CLOSED: "+str(exc),file=sys.stderr)
        raise
