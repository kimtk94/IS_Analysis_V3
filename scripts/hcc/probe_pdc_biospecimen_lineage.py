#!/usr/bin/env python3
"""HCC PDC protein-phosphoprotein PUBLIC biospecimen lineage audit.

Uses documented PDC paginatedCasesSamplesAliquots and only bounded first-header
prefixes of quantitative TSV. No identifiers, protein expression, raw data, or
signed URLs written to disk/logs. This does NOT join RNA expression to protein.
"""
import argparse
import collections
import datetime as dt
import json
import pathlib
import re

from probe_pdc_quant_headers import PDC_FILES, get_file_link, bounded_prefix
from probe_quant_matrix_manifest import request_graphql

STUDIES = ("PDC000198", "PDC000199")
BATCH_SIZE = 250
MAX_CASES = 1500

def canon(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())

def get_cases(study):
    records = []
    offsets = []
    total = None
    for offset in range(0, MAX_CASES, BATCH_SIZE):
        query = (
            '{ paginatedCasesSamplesAliquots(pdc_study_id:"' + study
            + '" offset:' + str(offset) + ' limit:' + str(BATCH_SIZE)
            + '){total casesSamplesAliquots {case_id case_submitter_id'
            + ' samples {sample_id sample_submitter_id sample_type'
            + ' aliquots {aliquot_id aliquot_submitter_id'
            + ' aliquot_run_metadata {label experiment_number alias fraction replicate_number analyte}}}}}}'
        )
        data = request_graphql(query).get("paginatedCasesSamplesAliquots") or {}
        if isinstance(data, list):
            data = data[0] if data else {}
        if not isinstance(data, dict):
            raise ValueError("Unknown paginatedCasesSamplesAliquots response")
        batch = data.get("casesSamplesAliquots") or []
        if not isinstance(batch, list):
            raise ValueError("casesSamplesAliquots not a list")
        raw_total = data.get("total")
        if total is None:
            try:
                total = int(raw_total)
            except (ValueError, TypeError):
                raise ValueError("missing numeric total") from None
        records.extend(batch)
        offsets.append({"offset":offset,"records":len(batch)})
        if len(records) >= total or len(batch) < BATCH_SIZE:
            break
    if len(records) < total:
        raise ValueError("Incomplete pagination; refusing case-count claim")
    return records, {"repository_reported_cases":total,"returned_records":len(records),"pages":offsets}

def flatten(records):
    """Build temporary in-memory case/sample/aliquot aliases for lineage mapping."""
    links=[]
    for case in records:
        if not isinstance(case,dict):
            continue
        case_keys=[case.get("case_id"),case.get("case_submitter_id")]
        cid=case.get("case_id") or case.get("case_submitter_id")
        for sample in case.get("samples") or []:
            if not isinstance(sample,dict):
                continue
            sample_keys=[sample.get("sample_id"),sample.get("sample_submitter_id")]
            sid=sample.get("sample_id") or sample.get("sample_submitter_id")
            for aliquot in sample.get("aliquots") or []:
                if not isinstance(aliquot,dict):
                    continue
                ak=[aliquot.get("aliquot_id"),aliquot.get("aliquot_submitter_id")]
                for run in aliquot.get("aliquot_run_metadata") or []:
                    if isinstance(run,dict):
                        ak.extend([run.get("alias"),run.get("experiment_number")])
                links.append({"case":cid,"sample":sid,"aliquot":aliquot.get("aliquot_id") or aliquot.get("aliquot_submitter_id"),
                              "type":sample.get("sample_type"),
                              "case_keys":[x for x in case_keys if x],
                              "sample_keys":[x for x in sample_keys if x],
                              "aliquot_keys":[x for x in ak if x]})
    return links

def count_summary(records,links):
    all_cases={x.get("case_id") or x.get("case_submitter_id") for x in records}
    all_cases.discard(None)
    c=collections.Counter(x.get("type") or "MISSING_TYPE" for x in links)
    case_tumor={x["case"] for x in links if x.get("type")=="Primary Tumor"}
    case_normal={x["case"] for x in links if x.get("type")=="Solid Tissue Normal"}
    cases_with_aliquots={x["case"] for x in links}
    bytype={}
    for typ in sorted(c):
        matching=[x for x in links if (x.get("type") or "MISSING_TYPE")==typ]
        bytype[typ]={"unique_cases":len({x["case"] for x in matching}),
                     "unique_samples":len({x["sample"] for x in matching}),
                     "aliquots":len({x["aliquot"] for x in matching})}
    run_rows=0
    distinct_experiments=set()
    experiment_labels=set()
    normalized_labels=collections.Counter()
    no_experiment_aliquots=0
    for case in records:
        for sample in case.get("samples") or []:
            for aliquot in sample.get("aliquots") or []:
                runs=aliquot.get("aliquot_run_metadata") or []
                if not runs:
                    no_experiment_aliquots+=1
                for record in runs:
                    if not isinstance(record,dict):continue
                    run_rows+=1
                    experiment=record.get("experiment_number")
                    label=record.get("label")
                    if experiment is not None:distinct_experiments.add(str(experiment))
                    if experiment is not None and label:
                        experiment_labels.add((str(experiment),str(label)))
                    if label:
                        normalized_labels[str(label).lower()]+=1
    return {"distinct_cases":len(all_cases),"distinct_samples":len({x["sample"] for x in links}),
            "distinct_aliquots":len({x["aliquot"] for x in links}),
            "cases_with_aliquot":len(cases_with_aliquots),
            "cases_without_aliquot":len(all_cases-cases_with_aliquots),
            "cases_with_tumor_and_normal_sample":len(case_tumor & case_normal),
            "tumor_only_case_count":len(case_tumor-case_normal),
            "normal_only_case_count":len(case_normal-case_tumor),
            "run_metadata_records":run_rows,
            "aliquots_without_run_metadata":no_experiment_aliquots,
            "distinct_experiment_numbers":len(distinct_experiments),
            "distinct_experiment_label_pairs":len(experiment_labels),
            "tmt_labels_aggregated":dict(normalized_labels),
            "sample_types":bytype}

def header_linkage(study, links):
    fname=PDC_FILES[study][0]
    _,url=get_file_link(study,fname)
    if not url:
        raise ValueError("No bounded download URL")
    payload,_,_=bounded_prefix(url)
    head=payload.decode("utf-8-sig",errors="replace").splitlines()[0].split("\t")
    # Match full well-defined identifier aliases only; short labels would give false positives.
    matchsets={"case":set(),"sample":set(),"aliquot":set()}
    columns=collections.Counter()
    for j,h in enumerate(head):
        h_norm=canon(h)
        if not h_norm:
            continue
        matched_types=set()
        for link in links:
            for kind in ("case","sample","aliquot"):
                aliases=link[kind+"_keys"]
                if any(len(canon(alias))>=8 and canon(alias) in h_norm for alias in aliases):
                    matchsets[kind].add(link[kind])
                    matched_types.add(kind)
        for kind in matched_types:
            columns[kind]+=1
    return {"matrix_header_column_count":len(head),"matched_columns_by_alias_type":dict(columns),
            "unique_case_aliases_recognized":len(matchsets["case"]),
            "unique_sample_aliases_recognized":len(matchsets["sample"]),
            "unique_aliquot_aliases_recognized":len(matchsets["aliquot"]),
            "interpretation":"Full ID token containment in bounded matrix headers only; not a gene-by-case quant join",
            "matrix_gene_expression_values_verified":False,
            "rna_protein_join_verified":False}

def audit():
    runs={}
    memory={}
    for study in STUDIES:
        try:
            records,meta=get_cases(study)
            links=flatten(records)
            memory[study]=links
            runs[study]={"status":"PUBLIC_BIOSPECIMEN_LINEAGE_FETCHED",
                         "pagination":meta,"lineage":count_summary(records,links)}
            try:
                runs[study]["matrix_header_correspondence"]=header_linkage(study,links)
            except Exception as e:
                runs[study]["matrix_header_error_type"]=type(e).__name__
        except Exception as e:
            runs[study]={"status":"UNVERIFIED","error_type":type(e).__name__}
    if len(memory)==len(STUDIES):
        a=memory[STUDIES[0]]
        b=memory[STUDIES[1]]
        ca={x["case"] for x in a};cb={x["case"] for x in b}
        sa={x["sample"] for x in a};sb={x["sample"] for x in b}
        aa={x["aliquot"] for x in a};ab={x["aliquot"] for x in b}
        runs["cross_study_overlap"]={
            "distinct_cases_in_both":len(ca&cb),
            "distinct_samples_in_both":len(sa&sb),
            "distinct_aliquots_in_both":len(aa&ab),
            "note":"PDC protein↔phosphoprotein study overlap only; neither RNA nor quantitative joined values proven."
        }
    return runs

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/hcc/feasibility/pdc_biospecimen_lineage.json")
    args=p.parse_args()
    records=audit()
    summary={"checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
             "source":"https://pdc.cancer.gov/graphql",
             "audit_type":"PUBLIC_PDC_CASE_SAMPLE_ALIQUOT_LINEAGE",
             "records":records,
             "V02_RNA_PROTEIN_QUANT_JOIN":"NOT_VERIFIED",
             "warning":"PDC study overlap and header alias evidence cannot establish RNA/protein paired patient quantification."}
    target=pathlib.Path(args.output)
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
