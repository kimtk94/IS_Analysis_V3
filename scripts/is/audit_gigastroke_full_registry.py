#!/usr/bin/env python3
"""Audit full GIGASTROKE accession registry against actually usable local files.
Never assume unknown ancestry can be inferred from adjacent accession IDs.
"""
import argparse
import csv
import json
from pathlib import Path
from collections import Counter

META=Path("/srv/is-analysis/data/is/reference/gigastroke/metadata/GIGASTROKE_STUDY_MAP_V2.tsv")
DATA=Path("/srv/is-analysis/data/is")
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--registry",type=Path,default=META)
    p.add_argument("--data",type=Path,default=DATA)
    p.add_argument("--root",type=Path,default=ROOT)
    a=p.parse_args()
    with a.registry.open() as f: entries=list(csv.DictReader(f,delimiter="\t"))
    results=[]
    for r in entries:
        acc=r["accession"]
        trait=r["phenotype_class"]
        p=a.data/"processed/gigastroke/eas"/(acc+"_"+trait+"_GRCh37.canonical.tsv.gz")
        is_canonical=p.is_file()
        if is_canonical:
            pop="EAS_PATH_CONFIRMED_VERIFY_SOURCE_METADATA"
        elif r["ancestry_class"]=="UNKNOWN":
            pop="UNRESOLVED_SOURCE_METADATA"
        else:
            pop=r["ancestry_class"]+"_REGISTRY_LABEL_ONLY"
        results.append({"accession":acc,"phenotype":trait,
          "metadata_ancestry":r["ancestry_class"],"population_resolution":pop,
          "sample_size_reported":r.get("sample_size",""),
          "canonical_local":int(is_canonical),
          "canonical_path":str(p) if is_canonical else "",
          "publication_source":"GIGASTROKE_2022_GWAS_CATALOG",
          "genome_build_expected":"GRCh37",
          "fine_mapping_state":"NOT_READY_REQUIRE_LD_AND_QC",
          "download_stage":"CANONICAL_AVAILABLE" if is_canonical else "NOT_PRESENT_AS_EAS_CANONICAL",
          "phenotype_overlap_caution":"NOT_INDEPENDENT_BY_DEFAULT"})
    fields=list(results[0])
    with (a.root/"GIGASTROKE_FULL_STUDY_INVENTORY.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=fields)
        w.writeheader();w.writerows(results)
    summary={"registry_accessions":len(results),
      "local_eas_canonical":sum(x["canonical_local"] for x in results),
      "not_local_eas_canonical":sum(not x["canonical_local"] for x in results),
      "metadata_ancestry_classes":dict(Counter(x["metadata_ancestry"] for x in results)),
      "phenotypes":dict(Counter(x["phenotype"] for x in results)),
      "note":"No unknown-population row assigned an inferred ancestry; no new source downloaded."}
    (a.root/"GIGASTROKE_STUDY_INVENTORY_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=="__main__":main()
