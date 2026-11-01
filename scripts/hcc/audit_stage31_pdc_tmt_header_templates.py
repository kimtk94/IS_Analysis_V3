#!/usr/bin/env python3
"""Stage31 bounded PDC quant header↔source TMT channel identity template audit.

This is an exploratory metadata comparison, NOT certification of correct
sample identities. Only exact, unique, source-derived run×channel candidates
are allowed (no arbitrary row-order matching). IDs, header values, gene
values, and individualized hashes never leave ephemeral GitHub runner memory.
"""
import argparse,collections,datetime as dt,json,pathlib,re
from probe_stage30_gao_catalog import PDC,pdc_query
from probe_stage30_gao_quant_api_header_to_biospecimen import get_quant_header
from probe_stage30_gao_experimental_design import LABELS,FIELD,QUERY,BIOS

def variants(label):
    # Explicit TMT11 source reporter channel: e.g. tmt_127n.
    suffix=label.removeprefix("tmt_")
    return tuple(dict.fromkeys((label,suffix,suffix.upper(),"TMT-"+suffix.upper(),"TMT_"+suffix.upper())))

def candidates(group, channel):
    # The only permitted canonicalizations are delimiters; no fuzzy edit
    # distance, no inference from same matrix position or patient number.
    for tag in variants(channel):
        for sep in ("_", "-", ".", ":", "/", " ", ""):
            yield str(group)+sep+tag

def run_design(pid):
    data,digest=pdc_query(QUERY%(pid,FIELD))
    rows=data.get("studyExperimentalDesign")
    if not isinstance(rows,list) or not 10<=len(rows)<=100:
        raise ValueError("PDC original 33 run design missing")
    design=[]
    for row in rows:
        plex=str(row.get("plex_dataset_name") or "").strip()
        run=str(row.get("study_run_metadata_submitter_id") or "").strip()
        for channel in LABELS:
            contents=row.get(channel) or []
            if isinstance(contents,dict):contents=[contents]
            for al in contents:
                aid=al.get("aliquot_id")
                if not isinstance(aid,str) or not aid.strip():
                    raise ValueError("source design aliquot identity missing")
                design.append({"plex":plex,"run":run,"channel":channel,
                               "aliquot":aid.strip()})
    if len(design)!=363 or len({x["aliquot"] for x in design})!=331:
        raise ValueError("PDC original design membership changed")
    multiplicity=collections.Counter(x["aliquot"] for x in design)
    reference_ids={k for k,v in multiplicity.items() if v>1}
    if len(reference_ids)!=1 or sum(multiplicity.values())!=363:
        raise ValueError("PDC TMT reference expected one reused aliquot")
    return design,reference_ids,digest

def rule_candidates(row, rule):
    field,_,casefold=rule.partition("|")
    value=row.get(field,"")
    if not value:return ()
    vals=set(candidates(value,row["channel"]))
    if casefold=="casefold":return (x.casefold() for x in vals)
    return vals

def score(headers,design,refs):
    headers=list(headers)
    # A transformation may map many TMT reporter assignments to identical
    # headers; only a unique candidate->aliquot mapping can be actionable.
    rules=("plex|exact","run|exact","plex|casefold","run|casefold")
    details=[]
    for rule in rules:
        field,_,fold=rule.partition("|")
        h=[x.casefold() if fold=="casefold" else x for x in headers]
        header_set=set(h)
        alias_to_aliquot=collections.defaultdict(set)
        alias_to_assignments=collections.Counter()
        for row in design:
            for token in rule_candidates(row,rule):
                if token in header_set:
                    alias_to_aliquot[token].add(row["aliquot"])
                    alias_to_assignments[token]+=1
        unique_alias={a:next(iter(s)) for a,s in alias_to_aliquot.items() if len(s)==1}
        unique_nonref={a:v for a,v in unique_alias.items() if v not in refs}
        all_nonref_aliquots=set(x["aliquot"] for x in design if x["aliquot"] not in refs)
        details.append({
            "source_alias_rule":field+"_"+fold,
            "quant_header_column_n":len(headers),
            "quant_headers_with_template_match_n":len(alias_to_aliquot),
            "quant_headers_with_single_design_aliquot_n":len(unique_alias),
            "distinct_unique_nonreference_aliquots_n":len(set(unique_nonref.values())),
            "source_nonreference_aliquot_pool_n":len(all_nonref_aliquots),
            "matched_template_tokens_involving_reference_aliquot_n":
               sum(any(q in refs for q in vv) for vv in alias_to_aliquot.values()),
            "ambiguous_alias_multiple_aliquots_n":
               sum(len(vv)>1 for vv in alias_to_aliquot.values()),
            "one_to_one_all_330_nonreference_aliquots_candidate":
               len(alias_to_aliquot)==len(unique_nonref)==
               len(set(unique_nonref.values()))==len(all_nonref_aliquots)==330
               and not any(any(q in refs for q in vv) for vv in alias_to_aliquot.values()),
            "case_or_aliquot_ids_exported":False,
        })
    return details

def study(pid):
    o={"pdc_study_id":pid,"status":"ACCESS_HOLD",
       "quant_column_ids_exported":False,
       "case_or_aliquot_ids_exported":False,
       "source_file_or_gene_values_exported":False,
       "biological_match_certified":False}
    try:
        summary,header=get_quant_header(pid)
        design,refs,dig=run_design(pid)
        rules=score(header,design,refs)
        o.update({"status":"SOURCE_DERIVED_TMT_RULES_TESTED",
                  "source_experimental_design_SHA256":dig,
                  "quant_column_count_n":len(header),
                  "original_design_reporter_assignment_n":len(design),
                  "reused_reference_aliquot_n":len(refs),
                  "candidate_transformations":rules,
                  "one_to_one_source_identity_certified":False})
    except Exception as e:
        o["failure_category"]=type(e).__name__
    return o

def review(v):
    if v.get("schema")!="hcc-stage31-PDC-direct-quant-TMT-template-audit-v1":
        raise ValueError("source schema drift")
    rows=v.get("studies")
    if not isinstance(rows,list) or {q.get("pdc_study_id") for q in rows}!=set(PDC):
        raise ValueError("wrong PDC targets")
    for q in rows:
        for k in ("quant_column_ids_exported","case_or_aliquot_ids_exported",
                  "source_file_or_gene_values_exported","biological_match_certified"):
            if q.get(k) is not False:raise ValueError("unapproved biological identity claim")
        if q["status"]=="SOURCE_DERIVED_TMT_RULES_TESTED":
            if q.get("quant_column_count_n")!=330 or q.get("original_design_reporter_assignment_n")!=363:
                raise ValueError("source TMT design drift")
            for rule in q["candidate_transformations"]:
                if not 0<=rule["quant_headers_with_single_design_aliquot_n"]<=330:
                    raise ValueError("impossible template counts")
                if rule.get("case_or_aliquot_ids_exported") is not False:
                    raise ValueError("ID leak")
        elif q["status"]!="ACCESS_HOLD":raise ValueError("unexpected access status")
    return {"schema":"hcc-stage31-TMT-source-header-templates-postflight-v1",
      "studies_with_exploratory_templates_checked_n":
          sum(x["status"]=="SOURCE_DERIVED_TMT_RULES_TESTED" for x in rows),
      "no_permutation_or_fuzzy_ID_match_performed":True,
      "G30_B":"HOLD_UNVERIFIED_SOURCE_QUANT_COLUMNS_TO_BIOPSY_IDS",
      "G30_C":"HOLD_ORIGINAL_NODE_RNA_SAMPLE_ID",
      "G30_D":"STOP_NO_MATCHED_EXTERNAL_RNA_PROTEIN_COHORT",
      "external_patient_omics_effects_computed":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    a=p.parse_args()
    out={"schema":"hcc-stage31-PDC-direct-quant-TMT-template-audit-v1",
         "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
         "studies":[study(x) for x in PDC]}
    f=pathlib.Path(a.out);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(out,indent=2)+"\n")
    summary=review(out)
    g=pathlib.Path(a.review_out);g.parent.mkdir(parents=True,exist_ok=True)
    g.write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
    print("STAGE31_TMT_TEMPLATES",json.dumps(out["studies"]))
if __name__=="__main__":main()
