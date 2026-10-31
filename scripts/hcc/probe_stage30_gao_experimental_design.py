#!/usr/bin/env python3
"""Official PDC TMT plex/channel ↔ aliquot experimental design metadata.

No patient, plex or aliquot IDs exported or printed. Source assay design does
NOT prove normalized protein quantity or matched transcriptomic patient joins.
"""
import argparse,collections,datetime as dt,json,pathlib
from probe_stage30_gao_catalog import PDC,pdc_query
LABELS=("label_free","tmt_126","tmt_127n","tmt_127c","tmt_128n","tmt_128c",
        "tmt_129n","tmt_129c","tmt_130n","tmt_130c","tmt_131","tmt_131c",
        "tmt_132n","tmt_132c","tmt_133n","tmt_133c","tmt_134n","tmt_134c",
        "tmt_135n","itraq_113","itraq_114","itraq_115","itraq_116",
        "itraq_117","itraq_118","itraq_119","itraq_121")
FIELD=" ".join(x+"{aliquot_id aliquot_submitter_id aliquot_run_metadata_id}" for x in LABELS)
QUERY='{studyExperimentalDesign(pdc_study_id:"%s"){pdc_study_id plex_dataset_name study_run_metadata_submitter_id %s}}'
BIOS='{biospecimenPerStudy(pdc_study_id:"%s"){aliquot_id aliquot_submitter_id sample_id sample_submitter_id case_id case_submitter_id sample_type}}'

def source(pid):
    out={"pdc_study_id":pid,"status":"ACCESS_HOLD",
         "source_plex_channel_aliquot_identifiers_exported":False,
         "quant_file_column_to_channel_mapping_verified":False,
         "RNA_same_patient_cohort_join_verified":False}
    private=None
    try:
        b,_=pdc_query(BIOS%pid)
        bs=b.get("biospecimenPerStudy")
        if not isinstance(bs,list) or not 1<=len(bs)<=10000:raise ValueError("missing PDC biospecimen manifest")
        by_id={str(x.get("aliquot_id") or "").strip() for x in bs if x.get("aliquot_id")}
        by_submitter={str(x.get("aliquot_submitter_id") or "").strip() for x in bs if x.get("aliquot_submitter_id")}
        if len(by_id)!=331 or len(by_submitter)!=331:raise ValueError("PDC source aliquot IDs changed")
        data,digest=pdc_query(QUERY%(pid,FIELD))
        designs=data.get("studyExperimentalDesign")
        if not isinstance(designs,list) or not 0<len(designs)<5000:
            raise ValueError("PDC experimental design not returned")
        plex=set();run=set();channels=collections.Counter()
        seen_assigned_ids=set();seen_assigned_submitter=set()
        assignments=0;empty=0
        unique_channel_pairs=set()
        for d in designs:
            if not isinstance(d,dict):continue
            p=str(d.get("plex_dataset_name") or "").strip()
            rn=str(d.get("study_run_metadata_submitter_id") or "").strip()
            if p:plex.add(p)
            if rn:run.add(rn)
            for lab in LABELS:
                rows=d.get(lab) or []
                if isinstance(rows,dict):rows=[rows]
                if not isinstance(rows,list):raise ValueError("design channel schema")
                for x in rows:
                    if not isinstance(x,dict):continue
                    assignments+=1
                    channels[lab]+=1
                    ai=str(x.get("aliquot_id") or "").strip()
                    bi=str(x.get("aliquot_submitter_id") or "").strip()
                    if ai:seen_assigned_ids.add(ai)
                    if bi:seen_assigned_submitter.add(bi)
                    if not ai and not bi:empty+=1
                    unique_channel_pairs.add((p or rn,lab))
        private={"plexes":plex,"run_names":run,
                 "aliquot_ids":seen_assigned_ids,"aliquot_submitter_ids":seen_assigned_submitter,
                 "channel_assignments":unique_channel_pairs}
        out.update({"status":"OFFICIAL_PDC_EXPERIMENTAL_DESIGN_VERIFIED",
            "source_response_sha256":digest,"source_design_run_records_n":len(designs),
            "unique_plex_names_n":len(plex),
            "unique_run_metadata_names_n":len(run),
            "total_channel_aliquot_assignments_n":assignments,
            "unique_aliquot_IDs_in_design_n":len(seen_assigned_ids),
            "unique_aliquot_submitter_IDs_in_design_n":len(seen_assigned_submitter),
            "aliquot_IDs_also_in_official_biospecimen_n":len(seen_assigned_ids&by_id),
            "aliquot_submitter_IDs_also_in_biospecimen_n":len(seen_assigned_submitter&by_submitter),
            "distinct_plex_TMT_channel_design_positions_n":len(unique_channel_pairs),
            "missing_aliquot_identity_assignments_n":empty,
            "channel_assignment_counts_by_label":dict(sorted(channels.items())),
            "source_quant_matrix_channel_alias_verified":False})
    except Exception as exc:
        out["failure_type"]=type(exc).__name__
    return out,private

def compare(a,b):
    if a is None or b is None:return {"status":"CROSS_DESIGN_ACCESS_HOLD",
        "case_sample_aliquot_identifiers_exported":False,
        "independent_RNA_protein_join":False}
    return {"status":"PDC_198_199_DESIGN_AGGREGATE",
        "cross_study_common_plex_run_names_n":len(a["plexes"]&b["plexes"]),
        "cross_study_common_run_metadata_name_n":len(a["run_names"]&b["run_names"]),
        "cross_study_common_designed_aliquot_UUID_n":len(a["aliquot_ids"]&b["aliquot_ids"]),
        "cross_study_common_designed_aliquot_submitter_n":len(a["aliquot_submitter_ids"]&b["aliquot_submitter_ids"]),
        "cross_study_common_plex_channel_positions_n":len(a["channel_assignments"]&b["channel_assignments"]),
        "case_sample_aliquot_identifiers_exported":False,
        "independent_RNA_protein_join":False}

def verify(o):
    if (o.get("schema")!="hcc-stage30-gao-tmt-experimental-design-v1" or
        {x.get("pdc_study_id") for x in o.get("studies",[])}!=set(PDC)):
        raise ValueError("source design contract")
    for x in o["studies"]:
        if x.get("status") not in ("ACCESS_HOLD","OFFICIAL_PDC_EXPERIMENTAL_DESIGN_VERIFIED"):
            raise ValueError("invalid design status")
        for k in ("source_plex_channel_aliquot_identifiers_exported",
                  "quant_file_column_to_channel_mapping_verified",
                  "RNA_same_patient_cohort_join_verified",
                  "source_quant_matrix_channel_alias_verified"):
            if x.get(k) is True:raise ValueError("forbidden source join or ID leak")
        if x["status"]=="OFFICIAL_PDC_EXPERIMENTAL_DESIGN_VERIFIED":
            for key in ("total_channel_aliquot_assignments_n","aliquot_IDs_also_in_official_biospecimen_n",
                        "aliquot_submitter_IDs_also_in_biospecimen_n"):
                v=x.get(key)
                if type(v) is not int or not 0<=v<=20000:raise ValueError("source assignment count invalid")
            if x["aliquot_IDs_also_in_official_biospecimen_n"]>331 or x["aliquot_submitter_IDs_also_in_biospecimen_n"]>331:
                raise ValueError("source assignment exceeds biospecimen universe")
    d=o.get("cross_study_source_design")
    if not isinstance(d,dict) or d.get("case_sample_aliquot_identifiers_exported") is not False or d.get("independent_RNA_protein_join") is not False:
        raise ValueError("invalid source cross-design guard")
    return {"schema":"hcc-stage30-tmt-design-validation-v1",
        "source_designs_verified_n":sum(x["status"]=="OFFICIAL_PDC_EXPERIMENTAL_DESIGN_VERIFIED" for x in o["studies"]),
        "G30_B":"PARTIAL_SOURCE_PLEX_ALIQUOT_DESIGN_METADATA_ONLY",
        "G30_C":"HOLD_RNA_NODE_QUANT_MATRIX",
        "G30_D":"STOP_NO_PROTEOME_RNA_PATIENT_JOIN",
        "G30_E":"NOT_STARTED",
        "external_gene_effect_replication_performed":False}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    a=p.parse_args()
    studies=[source(pid) for pid in PDC]
    out={"schema":"hcc-stage30-gao-tmt-experimental-design-v1",
         "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
         "source_endpoint":"PDC GraphQL studyExperimentalDesign",
         "studies":[s[0] for s in studies],
         "cross_study_source_design":compare(studies[0][1],studies[1][1])}
    dst=pathlib.Path(a.out);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(out,indent=2)+"\n")
    reviewed=verify(out);dst=pathlib.Path(a.review_out);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(reviewed,indent=2)+"\n")
    print(json.dumps(reviewed,indent=2))
    print("STAGE30_DESIGN_AGGREGATES",json.dumps(out["studies"]))
    print("STAGE30_DESIGN_CROSS",json.dumps(out["cross_study_source_design"]))
if __name__=="__main__":main()
