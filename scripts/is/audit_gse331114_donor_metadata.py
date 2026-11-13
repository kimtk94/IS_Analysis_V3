#!/usr/bin/env python3
"""Audit official human IS GEO sample metadata for donor/covariate completeness.

Only public GEO SOFT series metadata. No patient identifiers or values
retained in CI evidence. 8 different GSMs do NOT prove 8 unique donors.
Text fields are not a substitute for author-confirmed clinical study design.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re

from audit_stroke_versus_rosmap_sources import (
    geo_soft_url, download_metadata, parse_soft, derive_group,
)

FIELDS={
    "donor_id":re.compile(r"^(?:donor|donor[ _-]?id|participant[ _-]?id|patient[ _-]?id|subject[ _-]?id|individual[ _-]?id)$",re.I),
    "sex":re.compile(r"^(?:sex|gender|biological[ _-]?sex)$",re.I),
    "age":re.compile(r"^(?:age|age[ _-]?years|years[ _-]?old)$",re.I),
    "batch":re.compile(r"^(?:batch|sequencing[ _-]?batch|library[ _-]?batch|processing[ _-]?batch)$",re.I),
    "sampling_time":re.compile(r"^(?:time[ _-]?(?:after|from)[ _-]?(?:stroke|onset)|sampling[ _-]?time|stroke[ _-]?onset[ _-]?time|onset[ _-]?to[ _-]?collection)$",re.I),
    "lesion_site":re.compile(r"^(?:lesion[ _-]?site|brain[ _-]?region|tissue[ _-]?region|anatomic[ _-]?site)$",re.I),
}
MISSING={"none","null","n/a","na","unknown","not available","not applicable",""}
EXPECTED={"AIS":4,"Control":4}


def public_metadata_summary(info,samples):
    if info.get("accession")!="GSE331114" or len(samples)!=8:
        raise ValueError("Unexpected GEO accession/GSM count")
    records=[]
    counts=Counter()
    donor_values=[]
    keys=Counter()
    for gsm,sample in sorted(samples.items()):
        group=derive_group(sample)
        if group not in EXPECTED:
            raise ValueError("Unexpected GEO diagnosis label")
        counts[group]+=1
        values={}
        for raw in sample.get("characteristics",[]):
            if ":" not in raw:
                continue
            key,value=raw.split(":",1)
            k=key.strip().lower()
            keys[k]+=1
            for semantic,pattern in FIELDS.items():
                if pattern.fullmatch(k) and value.strip().lower() not in MISSING:
                    values.setdefault(semantic,[]).append(value.strip())
        donor=values.get("donor_id",[])
        donor_values.append(donor[0] if len(donor)==1 else None)
        records.append({
            "GSM":gsm,"group":group,
            "characteristic_key_names":[k.strip().lower() for k in (
                s.split(":",1)[0] for s in sample.get("characteristics",[]) if ":" in s)],
            "fields_with_explicit_values":sorted(values),
            "explicit_donor_id_present":len(donor)==1,
            "n_explicit_donor_id_values":len(donor),
            "n_sex_fields":len(values.get("sex",[])),
            "n_age_fields":len(values.get("age",[])),
            "n_batch_fields":len(values.get("batch",[])),
            "n_sampling_time_fields":len(values.get("sampling_time",[])),
            "n_lesion_site_fields":len(values.get("lesion_site",[])),
        })
    if dict(counts)!=EXPECTED:
        raise ValueError("GEO expected group counts are not 4 vs 4")
    n_explicit=sum(x is not None for x in donor_values)
    if n_explicit==8:
        n_unique=len(set(donor_values))
        donor_status=("EXPLICIT_IDS_UNIQUE_IN_GEO_NOT_BIOLOGICALLY_VERIFIED"
                      if n_unique==8 else "POSSIBLE_REPEATED_DONOR_IDS")
    else:
        n_unique=None
        donor_status=("NO_EXPLICIT_DONOR_IDS_IN_GEO" if n_explicit==0
                      else "PARTIAL_DONOR_IDS_IN_GEO")
    completeness={
        field:sum(field in r["fields_with_explicit_values"] for r in records)
        for field in FIELDS
    }
    return {
        "accession":"GSE331114",
        "source_organism":"Homo sapiens",
        "GEO_GSM_records":8,"GEO_groups":dict(counts),
        "metadata_schema_keys_with_occurrences":dict(keys),
        "available_clinical_field_GSM_counts":completeness,
        "n_explicit_donor_ids":n_explicit,
        "n_distinct_explicit_ids_without_values":n_unique,
        "independent_donor_evidence_status":donor_status,
        "biological_donor_independence_verified":False,
        "study_claim_biological_replicates_not_a_unique_donor_ID":True,
        "pseudobulk_DEG_permitted_from_this_metadata":False,
        "source_rows_without_identifier_values":records,
    }


def run(out):
    raw,compressed_bytes=download_metadata(geo_soft_url("GSE331114"))
    info,samples=parse_soft(raw)
    evidence=public_metadata_summary(info,samples)
    evidence["audit_utc"]=datetime.now(timezone.utc).isoformat()
    evidence["GEO_metadata_compressed_bytes"]=compressed_bytes
    out.mkdir(parents=True,exist_ok=True)
    (out/"sample_donor_covariate_coverage.json").write_text(json.dumps(evidence,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    lines=[
        "# GSE331114 GEO donor/covariate metadata completeness",
        "",
        f"**{evidence['independent_donor_evidence_status']}**; 8 GEO GSM records, AIS 4 / non-stroke 4.",
        "",
        "| Field | GSM records with explicit value |",
        "|---|---:|",
    ]
    for k,v in evidence["available_clinical_field_GSM_counts"].items():
        lines.append(f"| {k} | {v}/8 |")
    lines+=["",
        "Eight GSM accession numbers do not establish independent donors. "
        "Actual donor identifier values are deliberately omitted from artifacts.",
        "Pseudobulk DEG stays HOLD regardless of this metadata inspection; "
        "reference label validation and clinical/batch confounding remain unresolved.",
    ]
    summary="\n".join(lines)+"\n"
    (out/"summary.md").write_text(summary,encoding="utf-8")
    print(summary)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"],"a",encoding="utf-8") as f:
            f.write(summary)
    return evidence


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",default="artifacts/is_gse331114_donor_coverage")
    a=p.parse_args()
    run(Path(a.out_dir))


if __name__=="__main__":
    main()
