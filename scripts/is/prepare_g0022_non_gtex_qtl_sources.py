#!/usr/bin/env python3
"""Shortlist independent public QTL studies, preserving EBI source manifest metadata.
Does not presume sampled cohort ancestry; prioritization based on tissue and N only.
"""
import argparse,csv,json
from pathlib import Path

INPUT=Path("/tmp/g0022_eqtl_catalog_paths.tsv")
OUT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alternate_qtl_sources_v1")
PICKS={
 "QTD000609":("OneK1K CD14 monocytes","single_cell_immune"),
 "QTD000620":("OneK1K NK cells","single_cell_immune"),
 "QTD000434":("ROSMAP brain DLPFC","brain"),
 "QTD000051":("BrainSeq brain","brain"),
 "QTD000021":("BLUEPRINT monocytes","myeloid"),
 "QTD000110":("GEUVADIS LCL","multi_population_LCL")
}
def build(source,out):
    with source.open() as f:
        items={r["dataset_id"]:r for r in csv.DictReader(f,delimiter="\t")}
    if not set(PICKS).issubset(items):raise ValueError("Study accessions missing")
    rows=[]
    for ds,(description,reason) in PICKS.items():
        r=items[ds]
        if r["quant_method"]!="ge" or not r["ftp_path"].endswith(".all.tsv.gz"):
            raise ValueError("Not gene-level nominal source "+ds)
        url=r["ftp_path"].replace("ftp://ftp.ebi.ac.uk/","https://ftp.ebi.ac.uk/")
        if not url.startswith("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"):
            raise ValueError("Unexpected source host")
        rows.append({
            "dataset_id":ds,"study_id":r["study_id"],"study_label":r["study_label"],
            "sample_group":r["sample_group"],"sample_size":r["sample_size"],
            "description":description,"priority_relevance":reason,
            "nominal_summary_uri":url,
            "ancestry_verified":"NO","publicly_geolocated_variant_source":"GRCh38",
            "study_wide_independence_from_GWAS":"UNVERIFIED",
            "measurement_status":"NOT_YET_ASSAYED"})
    out.mkdir(parents=True,exist_ok=True)
    p=out/"G0022_NON_GTEX_QTL_SOURCE_CANDIDATES.tsv"
    with p.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader();w.writerows(rows)
    summary={"candidate_studies":len(rows),
       "study_labels":sorted({r["study_label"] for r in rows}),
       "ancestry_confirmed_EAS_cohorts":0,
       "source":"eQTL_Catalogue_official_tabix_manifest",
       "warning":"These are QTL dataset candidates, not ancestry matched or replication-confirmed."}
    (out/"G0022_NON_GTEX_QTL_SOURCE_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    return rows
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,default=INPUT)
    p.add_argument("--out",type=Path,default=OUT)
    a=p.parse_args()
    build(a.manifest,a.out)
