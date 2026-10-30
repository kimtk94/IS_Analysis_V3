#!/usr/bin/env python3
"""Audit PDC multiplex mapping using public Experimental Design API.

The studyExperimentalDesign endpoint maps (plex, TMT channel) to aliquots.
Only aggregated counts and header structural fingerprints are saved. Does NOT
establish RNA–protein patient matching or infer independent patients from labels.
"""
import argparse
import collections
import datetime as dt
import json
import pathlib
import re
from probe_quant_matrix_manifest import request_graphql
from probe_pdc_biospecimen_lineage import STUDIES, get_cases, flatten
from probe_pdc_quant_headers import PDC_FILES, get_file_link, bounded_prefix

TMT = ("tmt_126","tmt_127n","tmt_127c","tmt_128n","tmt_128c",
       "tmt_129n","tmt_129c","tmt_130n","tmt_130c","tmt_131","tmt_131c")
Q = "{studyExperimentalDesign(pdc_study_id:%s){pdc_study_id study_run_metadata_id study_run_metadata_submitter_id experiment_number plex_dataset_name %s}}"
CHANNEL_FIELDS = " ".join(f"{x}{{aliquot_id aliquot_submitter_id aliquot_run_metadata_id}}" for x in TMT)

def load_design(study):
    result=request_graphql(Q % (json.dumps(study),CHANNEL_FIELDS))
    rows=result.get("studyExperimentalDesign") or []
    if isinstance(rows,dict):
        rows=[rows]
    if not isinstance(rows,list):
        raise ValueError("studyExperimentalDesign unexpected return type")
    if not rows:
        raise ValueError("No experimental design records")
    return rows

def nested_items(row,key):
    item=row.get(key) or []
    return item if isinstance(item,list) else [item] if isinstance(item,dict) else []

def flatten_design(rows):
    slots=[]
    for r in rows:
        experiment=str(r.get("experiment_number") or "")
        plex=str(r.get("plex_dataset_name") or r.get("study_run_metadata_submitter_id") or "")
        for channel in TMT:
            for a in nested_items(r,channel):
                if not isinstance(a,dict):continue
                aid=a.get("aliquot_id") or a.get("aliquot_submitter_id")
                if aid:
                    slots.append((experiment,plex,channel,str(aid)))
    return slots

def unique_design_counts(rows,links):
    slots=flatten_design(rows)
    to_case={}
    to_type={}
    for l in links:
        for alias in l["aliquot_keys"]:
            to_case[str(alias)]=l["case"]
            to_type[str(alias)]=l["type"]
    plex_slots=collections.defaultdict(set)
    slot_to_aliquots=collections.defaultdict(set)
    observed=set()
    unmatched=set()
    case_set=set()
    category=collections.Counter()
    for exp,plex,channel,aid in slots:
        slot_to_aliquots[(exp,plex,channel)].add(aid)
        plex_slots[(exp,plex)].add(channel)
        observed.add(channel)
        case=to_case.get(aid)
        if case is None:
            unmatched.add(aid)
        else:
            case_set.add(case)
            category[to_type.get(aid) or "MISSING_TYPE"]+=1
    return {
        "experimental_design_rows":len(rows),
        "number_experiment_plex_pairs":len(plex_slots),
        "number_channel_slots":len(slot_to_aliquots),
        "number_channel_assignments":len(slots),
        "unique_aliquots_in_design":len({aid for _,_,_,aid in slots}),
        "unique_aliquots_not_in_biospecimen":len(unmatched),
        "unique_cases_with_design_assignments":len(case_set),
        "biological_channel_sample_types_assignments":dict(category),
        "ambiguous_channel_slots_multi_aliquot":sum(len(x)>1 for x in slot_to_aliquots.values()),
        "duplicate_channel_assignments":len(slots)-len({(ex,px,ch,a) for ex,px,ch,a in slots}),
        "unique_label_types_count":len(observed),
        "label_counts":dict(sorted(collections.Counter(ch for _,_,ch,_ in slots).items())),
        "channel_counts_per_plex_histogram":dict(sorted(collections.Counter(len(x) for x in plex_slots.values()).items())),
        "warning":"Channel may be a common reference or analytical sample; neither protein abundance quality nor RNA joining proven.",
    }

def normalized_header_fingerprint(study,design_rows):
    _,url=get_file_link(study,PDC_FILES[study][0])
    if not url:
        raise ValueError("No public quantitative TSV signed link")
    raw,_,_=bounded_prefix(url)
    header=raw.decode("utf-8-sig",errors="replace").splitlines()[0].split("\t")
    columns=header[1:]
    terms=collections.Counter()
    label_types=collections.Counter()
    experiment_names=set(str(x.get("experiment_number") or "").lower() for x in design_rows)
    plexnames=set(str(x.get("plex_dataset_name") or "").lower() for x in design_rows if x.get("plex_dataset_name"))
    for col in columns:
        lc=col.lower()
        if "unshared" in lc: terms["unshared_log_ratio"]+=1
        elif "log ratio" in lc or "log_ratio" in lc: terms["log_ratio"]+=1
        if "median" in lc:terms["median"]+=1
        for tag in TMT:
            if re.search(r"(?<![a-z0-9])"+re.escape(tag[4:])+r"(?![a-z0-9])",lc):
                label_types[tag]+=1
        if any(re.search(r"(?<![0-9])"+re.escape(x)+r"(?![0-9])",lc) for x in experiment_names if x):
            terms["mentions_experiment_number"]+=1
        if any(x in lc for x in plexnames if len(x)>=5):
            terms["mentions_plex_dataset_name"]+=1
    return {
        "matrix_columns_total":len(header),
        "quantitative_header_columns":len(columns),
        "annotated_ratio_field_counts":dict(terms),
        "recognized_tmt_label_column_counts":dict(label_types),
        "header_fully_read":bool(raw.count(b"\n")),
        "only_aggregate_fingerprint_saved":True,
        "verified_column_to_case_assignments":0,
        "RNA_PROTEIN_JOIN_VERIFIED":False,
    }

def reconstruct_from_biospecimen_runs(cases):
    """Reconstruct candidate plex/channel assignment from PDC aliquot_run_metadata.

    Plex name is an explicit ambiguity: experiment_number alone may group multiple
    analytical mixtures. Compare experiment-only and experiment+alias groupings,
    never claim a final quantified aliquot mapping solely from this reconstruction.
    """
    raw=[]
    for case in cases:
        for sample in case.get("samples") or []:
            for aliquot in sample.get("aliquots") or []:
                aid=aliquot.get("aliquot_id") or aliquot.get("aliquot_submitter_id")
                for run in aliquot.get("aliquot_run_metadata") or []:
                    if not isinstance(run,dict):continue
                    label=str(run.get("label") or "").lower()
                    exp=str(run.get("experiment_number") or "")
                    alias=str(run.get("alias") or "")
                    if aid and exp and label in TMT:
                        raw.append((exp,alias,label,str(aid)))
    if not raw:
        raise ValueError("No aliquot_run_metadata experiment-label assignments")
    def groups(use_alias):
        d={}
        for exp,alias,label,aid in raw:
            key=(exp,alias if use_alias else "")
            if key not in d:
                d[key]={"experiment_number":exp,"plex_dataset_name":key[1]}
            d[key].setdefault(label,[]).append({"aliquot_id":aid})
        return list(d.values())
    exp_rows=groups(False)
    alias_rows=groups(True)
    # Aliases only resolve a multiplex split if most alias groups contain
    # multiple independent TMT channels; single-channel aliases are insufficient.
    multis=sum(sum(bool(x.get(c)) for c in TMT)>=8 for x in alias_rows)
    choose_alias=len(alias_rows)>len(exp_rows) and multis>=max(1,len(alias_rows)//2)
    chosen=alias_rows if choose_alias else exp_rows
    support={
        "run_metadata_rows_with_valid_tmt_labels":len(raw),
        "distinct_experiment_numbers":len(exp_rows),
        "distinct_experiment_alias_groups":len(alias_rows),
        "experiment_alias_groups_with_8plus_labels":multis,
        "grouping_mode":"EXPERIMENT_AND_ALIAS" if choose_alias else "EXPERIMENT_ONLY",
        "alias_grouping_unvalidated":True,
        "experiment_only_slot_count":len(flatten_design(exp_rows)),
        "alias_group_slot_count":len(flatten_design(alias_rows)),
        "interpretation":"Aliquot-run metadata supplies candidate TMT assignments, not verified multiplex/file column linkage."
    }
    return chosen,support

def audit():
    o={}
    for study in STUDIES:
        result={"status":"UNVERIFIED","V02_RNA_PROTEIN_JOIN":"BLOCKED"}
        try:
            bios,_=get_cases(study)
            links=flatten(bios)
            design,provenance=reconstruct_from_biospecimen_runs(bios)
            result["run_metadata_reconstruction"]=provenance
            result["design"]=unique_design_counts(design,links)
            result["status"]="PDC_BIOSPECIMEN_RUN_METADATA_RECONSTRUCTED"
            try:result["matrix_header_fingerprint"]=normalized_header_fingerprint(study,design)
            except Exception as e:result["header_fingerprint_error_type"]=type(e).__name__
        except Exception as e:
            result["error_type"]=type(e).__name__
        o[study]=result
    return o

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="results/hcc/feasibility/pdc_experimental_design.json")
    args=p.parse_args()
    report={"checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "API":"https://pdc.cancer.gov/graphql",
            "records":audit(),
            "note":"PDC proteomics experimental design only, not RNA/protein matched quantitative expressions"}
    target=pathlib.Path(args.out)
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
