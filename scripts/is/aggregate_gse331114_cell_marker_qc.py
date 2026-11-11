#!/usr/bin/env python3
"""Integrate 8 real GSE331114 cell-quality / marker-support QC artifacts.

No cell-type labels, no donor-genotype linkage and no DEG. A sample-wise
summary prevents cells from being treated as independent cases.
"""
import argparse
import csv
from datetime import datetime,timezone
import json
from pathlib import Path
import statistics
import sys

SOURCE_RUN_ID=43013961145
EXPECTED={(group,i):f"GSM{9739251 + i + (0 if group=='AIS' else 4)}"
          for group in ("AIS","Control") for i in range(4)}
PANELS=("T_cell_marker_support","B_cell_marker_support","NK_marker_support",
        "myeloid_marker_support","microglia_marker_support",
        "MULTIPLE_PANELS","UNASSIGNED")


def load(root):
    rows=[]
    for (group,index),gsm in EXPECTED.items():
        folder=root/f"is-gse331114-marker-{group}-{index}"
        f=folder/"marker_qc.json"
        if not f.exists():
            raise ValueError("Missing 8/8 source marker QC artifact: "+str(f))
        item=json.loads(f.read_text(encoding="utf-8"))
        if (item.get("source")!="GSE331114" or item.get("GSM")!=gsm
                or item.get("group")!=group or item.get("index")!=index
                or item.get("status")!="SAMPLE_CELL_QC_AND_MARKER_SUPPORT_PASS"
                or item.get("source_file_sha256_matches_prior") is not True
                or item.get("not_a_celltype_annotation") is not True
                or item.get("no_p_values_or_degs") is not True
                or item.get("celltype_pseudobulk")!="NOT_PERFORMED"
                or item.get("donor_independence")!="NOT_VERIFIED_WITH_INDIVIDUAL_IDENTIFIERS"):
            raise ValueError("Source QC identity or inference guardrail failed: "+str(f))
        qc=item["marker_panel_support_unvalidated_counts_after_QC"]
        n=int(item["cells_pass_qc"])
        if n<1 or sum(qc.values())!=n:
            raise ValueError("Marker QC category total does not sum to pass barcodes")
        if int(item["total_barcodes"])!=n+int(item["cells_fail_qc"]):
            raise ValueError("Barcode QC pass/fail counts inconsistent")
        row={"GSM":gsm,"group":group,"index":index,"all_barcodes":item["total_barcodes"],
             "qc_passing_barcodes":n,"pct_pass_qc":item["pct_pass_qc"],
             "median_mito_percent":item["median_mito_percent_all_barcodes"],
             "median_n_genes":item["median_n_genes_all_barcodes"],
             "median_n_umis":item["median_n_umis_all_barcodes"],
             "known_sample_characteristic_keys":"|".join(item.get("metadata_characteristic_field_names",[]))}
        for p in PANELS:
            count=int(qc.get(p,0))
            row[p+"_barcode_count"]=count
            row[p+"_fraction_of_QC"]=round(count/n,6)
        rows.append(row)
    if len({r["GSM"] for r in rows})!=8:
        raise ValueError("Nonunique GEO GSMs")
    return rows


def descriptive(rows):
    results=[]
    for p in PANELS:
        a=[r[p+"_fraction_of_QC"] for r in rows if r["group"]=="AIS"]
        c=[r[p+"_fraction_of_QC"] for r in rows if r["group"]=="Control"]
        if len(a)!=4 or len(c)!=4:
            raise ValueError("Missing sample-level groups")
        results.append({
            "marker_panel_not_validated_cell_type":p,
            "n_AIS_GSM":4,"n_Control_GSM":4,
            "AIS_median_support_fraction":round(statistics.median(a),6),
            "Control_median_support_fraction":round(statistics.median(c),6),
            "AIS_min":min(a),"AIS_max":max(a),"Control_min":min(c),"Control_max":max(c),
            "descriptive_delta_median_fraction":round(statistics.median(a)-statistics.median(c),6),
            "p_value":"NOT_TESTED","fdr":"NOT_TESTED",
            "celltype_validation":"NOT_PERFORMED",
        })
    return results


def write_tsv(path, rows):
    with open(path,"w",newline="",encoding="utf-8") as out:
        writer=csv.DictWriter(out,fieldnames=list(rows[0]),delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def process(source,out):
    rows=load(source)
    result=descriptive(rows)
    out.mkdir(parents=True,exist_ok=True)
    write_tsv(out/"sample_marker_quality.tsv",rows)
    write_tsv(out/"marker_group_descriptive.tsv",result)
    known_keys=sorted(set("|".join(r["known_sample_characteristic_keys"] for r in rows).split("|"))-set([""]))
    valid_annotations=False
    donor_unique_ids_verified=False
    ready=valid_annotations and donor_unique_ids_verified
    evidence={
        "as_of_utc":datetime.now(timezone.utc).isoformat(),
        "provenance_run":SOURCE_RUN_ID,
        "all_eight_actual_source_artifacts_validated":True,
        "n_GSM":8,
        "n_AIS":4,"n_Control":4,
        "clinical_sample_metadata_fields_present":known_keys,
        "independent_donor_id_verified":donor_unique_ids_verified,
        "validated_author_or_reference_celltype_annotations":valid_annotations,
        "cell_type_pseudobulk_ready":ready,
        "cell_type_DEG":"NOT_RUN",
        "inferred_cell_type_labels":"NONE",
        "metric":"exclusive expression of at least two marker genes per predefined panel among QC passing barcodes",
        "status":"MARKER_CONTEXT_DESCRIPTIVE_ONLY_CELLTYPE_PSEUDOBULK_HOLD",
        "missing_required_evidence":[
            "original validated per-cell immune celltype labels or independently benchmarked reference annotation",
            "unambiguous biological donor independence and sex/age/batch clinical metadata",
            "per-sample and per-celltype QC minimums verified before pseudobulk construction",
            "gene-wide raw count pseudobulk inputs with study-approved design",
        ],
    }
    (out/"evidence.json").write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    lines=["# Human IS CD45+ brain: 8-sample marker QC","",
           "**8/8 original source files checked, per-cell QC profiles summarized by GSM.**",
           "",
           "| Marker co-detection (not cell-type) | AIS median QC barcode fraction | Control median | Difference |",
           "|---|---:|---:|---:|"]
    for r in result:
        lines.append(f"| {r['marker_panel_not_validated_cell_type']} | "
                     f"{r['AIS_median_support_fraction']:.3f} | "
                     f"{r['Control_median_support_fraction']:.3f} | "
                     f"{r['descriptive_delta_median_fraction']:+.3f} |")
    lines+=["","- **NOT cell-type annotation:** marker co-detection, no reference mapping or independent label validation.",
            "- **NOT a pseudobulk DEG:** group-level fractions descriptive across 4 AIS and 4 control GSMs only.",
            "- Patient uniqueness, confounding covariates, and differential expression significance remain UNVERIFIED.",
            "- Microglia marker down-regulation during activation can lead to severe false negatives.",
            "- Complete single-cell data remain in temporary CI runner; outputs contain GSM aggregate metrics only."]
    md="\n".join(lines)+"\n"
    (out/"summary.md").write_text(md,encoding="utf-8")
    print(md)
    import os
    if path:=os.getenv("GITHUB_STEP_SUMMARY"):
        with open(path,"a",encoding="utf-8") as f:
            f.write(md)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",default="downloaded_gse331114_markers")
    p.add_argument("--out-dir",default="artifacts/is_gse331114_marker_synthesis")
    args=p.parse_args()
    process(Path(args.source),Path(args.out_dir))


if __name__=="__main__":
    try:
        main()
    except (ValueError,OSError,KeyError,TypeError) as e:
        print("FAIL_CLOSED: "+str(e),file=sys.stderr)
        raise
