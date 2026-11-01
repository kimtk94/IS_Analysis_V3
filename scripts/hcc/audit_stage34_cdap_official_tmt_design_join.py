#!/usr/bin/env python3
"""Stage34: sample.txt AnalyticalSample ↔ official PDC TMT run/plex source join.

Only source-identical names, with separately reported casefold comparison,
may create provisional join candidates. Never join by row index. No patient ID,
file name, source row, sample alias or signed URL may be logged.
"""
import argparse
import collections
import datetime as dt
import hashlib
import json
import pathlib
import re

from probe_stage30_gao_catalog import PDC
from audit_stage33_original_cdap_sample_bridge import discover, read_bytes, parse_sample_table, RATIO_PATTERN
from audit_stage31_pdc_tmt_header_templates import run_design

def channel_key(label):
    x = str(label or "").strip().upper().replace("TMT_", "").replace("TMT-", "")
    return x if re.fullmatch(r"1\d\d[NC]?", x) else None

def design_index(design):
    by_plex = collections.defaultdict(dict)
    by_run = collections.defaultdict(dict)
    for x in design:
        k = channel_key(x["channel"])
        if not k:
            raise ValueError("unknown official TMT reporter")
        for group, dest in ((x["plex"], by_plex), (x["run"], by_run)):
            if not group:
                continue
            if k in dest[group] and dest[group][k] != x["aliquot"]:
                raise ValueError("official TMT reporter conflicts")
            dest[group][k] = x["aliquot"]
    return {"plex_exact": by_plex, "run_exact": by_run}

def analyze(rows, reporters, design, ref_ids):
    lookups = design_index(design)
    metric = {}
    for kind, index in lookups.items():
        for fold in (False, True):
            keyname = kind.replace("_exact","") + ("_casefold" if fold else "_exact")
            candidates = collections.defaultdict(set)
            for source_name in index:
                candidates[source_name.casefold() if fold else source_name].add(source_name)
            matched_groups, ambiguous_groups = set(), set()
            linked_aliquots = collections.Counter()
            reporter_source_labels = collections.defaultdict(set)
            official_source_reporter_positions_n = 0
            reference_assignments_n = 0
            missing_channel_n = 0
            ratios_parseable_n = 0
            source_labels_equal_aliquot_id_n = 0
            for row in rows:
                sample = (row.get("AnalyticalSample") or "").strip()
                possible = candidates.get(sample.casefold() if fold else sample,set())
                if len(possible) != 1:
                    if len(possible) > 1:
                        ambiguous_groups.add(sample)
                    continue
                name = next(iter(possible))
                matched_groups.add(name)
                source_channels = index[name]
                specs = [x.strip() for x in (row.get("Ratios") or "").split(",") if x.strip()]
                for spec in specs:
                    r = RATIO_PATTERN.fullmatch(spec)
                    if not r:
                        continue
                    ratios_parseable_n += 1
                    numerator = r.group(1).upper()
                    label_col = [c for c in reporters if channel_key(c) == numerator]
                    if len(label_col) != 1:
                        missing_channel_n += 1
                        continue
                    ali = source_channels.get(numerator)
                    if not ali:
                        missing_channel_n += 1
                        continue
                    official_source_reporter_positions_n += 1
                    linked_aliquots[ali] += 1
                    source_label = (row.get(label_col[0]) or "").strip()
                    if source_label:
                        reporter_source_labels[source_label].add(ali)
                        source_labels_equal_aliquot_id_n += (source_label == ali)
                    reference_assignments_n += ali in ref_ids
            unique_nonreference = {a for a, n in linked_aliquots.items() if a not in ref_ids and n==1}
            metric[keyname] = {
                "source_analytical_samples_uniquely_linked_to_official_run_or_plex_n":len(matched_groups),
                "source_analytical_samples_ambiguous_match_n":len(ambiguous_groups),
                "source_reporter_ratio_specs_parseable_on_joined_samples_n":ratios_parseable_n,
                "source_reporter_ratio_positions_with_official_aliquot_n":official_source_reporter_positions_n,
                "source_reporter_ratio_positions_with_missing_official_channel_n":missing_channel_n,
                "source_reporter_ratio_positions_assigning_reference_aliquot_n":reference_assignments_n,
                "distinct_official_nonreference_aliquots_at_unique_source_positions_n":len(unique_nonreference),
                "source_ID_reporter_labels_with_multiple_official_aliquots_n":
                    sum(len(z)>1 for z in reporter_source_labels.values()),
                "source_ID_labels_equal_official_aliquot_UUID_n":source_labels_equal_aliquot_id_n,
                "source_reporter_assignment_count_by_aliquot_not_exported": True,
                "provisional_330_of_330_run_channel_link_candidate":
                    len(matched_groups)==33 and official_source_reporter_positions_n==330
                    and len(unique_nonreference)==330 and reference_assignments_n==0
                    and missing_channel_n==0,
            }
    return metric

def one(pid):
    out = {
        "pdc_study_id":pid, "status":"ACCESS_HOLD",
        "source_names_ids_or_urls_exported":False, "quant_column_to_aliquot_certified":False,
        "RNA_protein_patient_join_certified":False,
    }
    try:
        files,metadata_sha = discover(pid)
        out["source_PDC_file_metadata_sha256"] = metadata_sha
        out["cdap_companion_candidate_files_n"]=len(files)
        if len(files)!=1:
            out["status"]="REVIEW_AMBIGUOUS_SOURCE_COMPANION"
            return out
        body=read_bytes(files[0])
        out["source_Canonical_CDAP_sample_content_sha256"] = hashlib.sha256(body).hexdigest()
        rows,reporters=parse_sample_table(body)
        design,refs,exp_sha=run_design(pid)
        out.update({
            "source_official_PDC_experimental_design_SHA256":exp_sha,
            "source_analytical_sample_rows_n":len(rows),
            "source_TMT_reporter_columns_n":len(reporters),
            "official_experiment_design_assignments_n":len(design),
            "official_reference_aliquots_n":len(refs),
            "comparison": analyze(rows,reporters,design,refs),
            "status":"CDAP_TO_OFFICIAL_TMT_DESIGN_IDENTITY_COMPARED",
        })
    except Exception as e:
        out["error_class"]=type(e).__name__
    return out

def verify(x):
    if x.get("schema")!="hcc-stage34-CDAP-to-PDC-TMT-source-identity-v1":
        raise ValueError("invalid source schema")
    d=x.get("source")
    if not isinstance(d,dict) or d.get("pdc_study_id") not in PDC:
        raise ValueError("wrong study")
    for k in ("source_names_ids_or_urls_exported","quant_column_to_aliquot_certified",
              "RNA_protein_patient_join_certified"):
        if d.get(k) is not False: raise ValueError("unsupported patient/source assertion")
    if d["status"]=="CDAP_TO_OFFICIAL_TMT_DESIGN_IDENTITY_COMPARED":
        if d["source_analytical_sample_rows_n"]!=33 or d["official_experiment_design_assignments_n"]!=363:
            raise ValueError("not the previously characterized Gao TMT structure")
        if set(d["comparison"])!={"plex_exact","plex_casefold","run_exact","run_casefold"}:
            raise ValueError("missing comparison rules")
        for v in d["comparison"].values():
            n=v["source_reporter_ratio_positions_with_official_aliquot_n"]
            if not 0<=n<=330 or not 0<=v["distinct_official_nonreference_aliquots_at_unique_source_positions_n"]<=330:
                raise ValueError("impossible match counts")
    elif d["status"] not in ("ACCESS_HOLD","REVIEW_AMBIGUOUS_SOURCE_COMPANION"):
        raise ValueError("invalid source access status")
    return {
        "schema":"hcc-stage34-CDAP-official-design-reviewed-v1",
        "direct_source_design_join_checked":d["status"]=="CDAP_TO_OFFICIAL_TMT_DESIGN_IDENTITY_COMPARED",
        "G30_B":"HOLD_NO_SOURCE_AUTHORIZED_330_QUANT_HEADER_TO_ALIQUOT_CERTIFICATION",
        "G30_C":"HOLD_NODE_RNA",
        "G30_D":"STOP_NO_MATCHED_EXTERNAL_RNA_PROTEIN",
        "biological_claims":False,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--study",choices=PDC,required=True)
    p.add_argument("--out",required=True)
    p.add_argument("--review-out",required=True)
    a=p.parse_args()
    v={"schema":"hcc-stage34-CDAP-to-PDC-TMT-source-identity-v1",
       "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
       "source":one(a.study)}
    status=verify(v)
    for filename,obj in ((a.out,v),(a.review_out,status)):
        path=pathlib.Path(filename)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(obj,indent=2)+"\n")
    print("HCC_STAGE34_SOURCE",json.dumps(v))
    print("HCC_STAGE34_REVIEW",json.dumps(status))

if __name__=="__main__":main()
