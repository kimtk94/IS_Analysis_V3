#!/usr/bin/env python3
"""Source-of-truth GIGASTROKE ancestry/size audit from per-accession GWAS Catalog YAML.

Study families are correlated meta-analyses, NOT independent replication cohorts.
"""
import argparse,csv,json
from pathlib import Path
from collections import Counter
import yaml
DEFAULT_META=Path("/srv/is-analysis/data/is/reference/gigastroke/metadata")
DEFAULT_DATA=Path("/srv/is-analysis/data/is")
DEFAULT_OUT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
FAMILIES={"East Asian":"EAS","European":"EUR","South Asian":"SAS",
          "Hispanic or Latin American":"HIS","African American or Afro-Caribbean":"AFR"}
def build(meta,data,out):
    with (meta/"GIGASTROKE_STUDY_MAP_V2.tsv").open(newline="") as f:
        registry=list(csv.DictReader(f,delimiter="\t"))
    output=[]
    for old in registry:
        aid=old["accession"]
        source=meta/"yaml"/(aid+"_buildGRCh37.tsv.gz-meta.yaml")
        if not source.is_file():raise FileNotFoundError(source)
        info=yaml.safe_load(source.read_text())
        samples=info.get("samples") or []
        records=[]
        for item in samples:
            categories=item.get("sample_ancestry_category") or []
            if isinstance(categories,str):categories=[categories]
            for category in categories:
                if category not in FAMILIES:raise ValueError("Unknown category: "+str(category))
                records.append((FAMILIES[category],int(item["sample_size"])))
        if not records:raise ValueError("Missing per-study source ancestry: "+aid)
        ancestry="MULTI_ANCESTRY" if len(set(x[0] for x in records))>1 else records[0][0]
        actual=sum(x[1] for x in records)
        old_n=int(old.get("sample_size") or 0)
        canonical=data/"processed/gigastroke/eas"/f'{aid}_{old["phenotype_class"]}_GRCh37.canonical.tsv.gz'
        raw=data/"reference/gigastroke/eas"/f'{aid}_buildGRCh37.tsv.gz'
        source_md5=info.get("data_file_md5sum","")
        if len(str(source_md5))!=32:raise ValueError("missing file md5: "+aid)
        if info.get("genome_assembly")!="GRCh37":raise ValueError("Unexpected genome assembly: "+aid)
        url=("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
              f'GCST90104001-GCST90105000/{aid}/{aid}_buildGRCh37.tsv.gz')
        output.append({"accession":aid,"phenotype":old["phenotype_class"],
          "ancestry_verified":ancestry,"sample_ancestries":";".join(x[0] for x in records),
          "subpopulation_n":";".join(f'{x[0]}:{x[1]}' for x in records),
          "source_total_n":actual,"registry_legacy_n":old_n,
          "legacy_n_mismatch":int(old_n not in (0,actual)),
          "source_md5":source_md5,"genome_build":"GRCh37","source_url":url,
          "raw_local":int(raw.is_file()),"eas_canonical_local":int(canonical.is_file()),
          "acquisition_priority":("EXISTING_EAS" if canonical.is_file()
              else "P1_EUR_DISCOVERY" if ancestry=="EUR"
              else "P2_CROSS_ANCESTRY" if ancestry=="MULTI_ANCESTRY"
              else "P3_SECONDARY_ANCESTRY"),
          "independent_replication":"NO_OVERLAP_UNKNOWN",
          "source_yaml":str(source)})
    out.mkdir(parents=True,exist_ok=True)
    with (out/"GIGASTROKE_SOURCE_VERIFIED_V3.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=list(output[0]))
        w.writeheader();w.writerows(output)
    counter=Counter(x["ancestry_verified"] for x in output)
    summary={"studies":len(output),"ancestry_breakdown":dict(sorted(counter.items())),
             "mislabeled_or_unresolved_legacy_rows":sum(
                 x["ancestry_verified"]=="MULTI_ANCESTRY" or x["legacy_n_mismatch"] or
                 (x["ancestry_verified"] in ("EUR","EAS","SAS") and
                  x["accession"] in [r["accession"] for r in registry if r["ancestry_class"]=="UNKNOWN"])
                 for x in output),
             "n_sample_size_mismatches":sum(x["legacy_n_mismatch"] for x in output),
             "eas_canonical_loaded":sum(x["eas_canonical_local"] for x in output),
             "european_candidates":sum(x["ancestry_verified"]=="EUR" for x in output),
             "cross_ancestry_candidates":sum(x["ancestry_verified"]=="MULTI_ANCESTRY" for x in output),
             "warning":"GIGASTROKE cohorts overlap across trait and meta-analysis families. Never count as independent validation."}
    (out/"GIGASTROKE_SOURCE_VERIFIED_V3_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    return summary
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--meta",type=Path,default=DEFAULT_META)
    p.add_argument("--data",type=Path,default=DEFAULT_DATA)
    p.add_argument("--out",type=Path,default=DEFAULT_OUT)
    a=p.parse_args()
    build(a.meta,a.data,a.out)
