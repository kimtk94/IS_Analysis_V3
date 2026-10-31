#!/usr/bin/env python3
"""Evidence *depth* report for HCC P0, distinct from GitHub execution status.

Reads only approved aggregate public-audit JSON from the six lane artifacts.
Does NOT serialize patient IDs, raw clinical fields, cohort rows or download
tokens. A successful script is not a verified quantitative RNA/protein join.
"""
import argparse
import datetime as dt
import json
import pathlib

MAX_BYTES=2_000_000

EVIDENCE_ITEMS=(
    ("ng_supplement", "ng-public", "ng_gao_metadata.json"),
    ("ng_quant_workbook_shape", "ng-public", "quant_matrix_manifest.json"),
    ("pdc_protein_sample_link", "pdc-mapping", "pdc_sample_mapping.json"),
    ("node_rna_access", "rna-access", "node_oep000321_access.json"),
    ("hccdb_rna_access", "rna-access", "hccdb25_access.json"),
    ("geo_metadata", "geo-clinical", "geo_metadata.json"),
    ("xing_clinical", "xing-clinical", "xing_clinical_complete_case.json"),
    ("xing_etiology", "xing-clinical", "xing_etiology_codebook_preflight.json"),
)

def public_json(root, lane, filename):
    p=pathlib.Path(root)/f"hcc-p0-{lane}"/filename
    if not p.is_file() or p.stat().st_size > MAX_BYTES:
        return None
    try:
        doc=json.loads(p.read_text(encoding="utf-8"))
        return doc if isinstance(doc,dict) else None
    except (OSError,UnicodeDecodeError,ValueError):
        return None

def summarize(root):
    docs={name:public_json(root,lane,filename) for name,lane,filename in EVIDENCE_ITEMS}
    ng=docs["ng_supplement"]
    ng_attempts=((ng.get("ng") or {}).get("attempts") or []) if ng else []
    ng_metadata=any(isinstance(a,dict) and
                    a.get("status") not in (None,"UNVERIFIED") and
                    isinstance(a.get("xlsx_schema"),dict) for a in ng_attempts)

    qm=docs["ng_quant_workbook_shape"]
    ng_workbooks=list((qm.get("ng_supplements") or {}).values()) if qm else []
    ng_shape_n=sum(isinstance(v,dict) and v.get("status") ==
                   "PUBLIC_XLSX_SHAPE_VERIFIED" for v in ng_workbooks)

    pm=docs["pdc_protein_sample_link"]
    studies=(pm.get("records") or {}) if pm else {}
    pdc_schema_n=sum(isinstance(v,dict) and v.get("status") ==
                     "PUBLIC_MATRIX_METADATA_CORRESPONDENCE_SCREENED"
                     for v in studies.values())
    pdc_any_header_match=sum(
        isinstance(v,dict) and
        (v.get("schema_comparison") or {}).get("any_exact_normalized_overlap") is True
        for v in studies.values())

    node=docs["node_rna_access"]
    node_record=(node.get("NODE") or {}) if node else {}
    node_landing=node_record.get("status")=="PUBLIC_LANDING_HTML_ONLY"
    node_quant_verified=node_record.get("downloaded_expression_matrix") is True

    hccdb=docs["hccdb_rna_access"]
    landing=list((hccdb.get("landing_checks") or {}).values()) if hccdb else []
    hccdb_landing=sum(isinstance(v,dict) and v.get("status")=="LANDING_HTML_ONLY"
                      for v in landing)
    hccdb_quant_verified=hccdb.get("matrix_download_status")=="VERIFIED" if hccdb else False

    geo=docs["geo_metadata"]
    geos=geo.get("records") if geo else None
    geo_fetched=sum(isinstance(v,dict) and
                    v.get("download_status")=="FETCHED_PUBLIC_SUPPLEMENT"
                    for v in geos) if isinstance(geos,list) else 0

    xing=docs["xing_clinical"]
    xa=(xing.get("audit") or {}) if xing else {}
    xing_152=(xa.get("status")=="MATCHED_COMPLETE_CASE_PRESCREEN_AUDITED"
              and xa.get("n_matched_structurally")==152
              and xa.get("source_data_sha256_verified") is True)
    xe=docs["xing_etiology"]
    et=(xe.get("audit") or {}) if xe else {}
    et_struct=(et.get("status")=="ETIOLOGY_CODEBOOK_PREFLIGHT_AUDITED"
               and et.get("matched_protein_to_clinical_rows")==152
               and et.get("source_checksums_verified") is True)
    et_codebook=(et.get("etiology_patient_group_assignment_confirmed") is True)

    evidence={
        "ng_supplement_structure": {
            "state":"SOURCE_STRUCTURE_VERIFIED" if ng_metadata else "UNVERIFIED",
            "observation":"Published public supplement workbook parser reached schema only",
            "patient_rna_protein_quantitative_join":False,
        },
        "ng_quant_workbook_index": {
            "state":"PUBLIC_WORKBOOK_SHAPE" if ng_shape_n else "UNVERIFIED",
            "workbooks_with_shape_verified":ng_shape_n,
            "quantitative_omics_join_verified":False,
        },
        "gao_pdc_protein_metadata": {
            "state":"MATRIX_SAMPLE_HEADERS_SCREENED" if pdc_schema_n else "UNVERIFIED",
            "studies_with_schema_checked":pdc_schema_n,
            "studies_with_exact_normalized_header_overlap":pdc_any_header_match,
            "same_patient_rna_and_protein_quantification":False,
        },
        "node_rna": {
            "state":"LANDING_ONLY" if node_landing else "UNVERIFIED",
            "downloaded_expression_matrix_claim_from_source":node_quant_verified,
            "patient_sample_join_established":False,
        },
        "hccdb_rna": {
            "state":"LANDING_ONLY" if hccdb_landing else "UNVERIFIED",
            "responsive_landing_count":hccdb_landing,
            "verified_processed_expression_matrix_claim_from_source":hccdb_quant_verified,
            "patient_sample_join_established":False,
        },
        "geo": {
            "state":"SOURCE_METADATA_FETCHED" if geo_fetched else "UNVERIFIED",
            "public_supplements_fetched":geo_fetched,
            "validated_recurrence_event_definition":False,
        },
        "xing_protein_clinical": {
            "state":"152_SOURCE_MATCHED_FORMAT_QC" if xing_152 else "UNVERIFIED",
            "structural_matched_pairs":152 if xing_152 else None,
            "clinical_endpoint_codebook_verified":False,
        },
        "xing_etiology": {
            "state":"152_SOURCE_CODE_SHAPES" if et_struct else "UNVERIFIED",
            "confirmed_etiology_patient_group_assignment":et_codebook,
        },
    }
    all_present=all(value is not None for value in docs.values())
    return {
        "schema_version":"hcc-p0-evidence-depth-v1",
        "created_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "expected_source_reports":len(EVIDENCE_ITEMS),
        "source_reports_read":sum(v is not None for v in docs.values()),
        "all_source_reports_read":all_present,
        "evidence":evidence,
        "research_claim_gates":{
            "V02_patient_matched_quantitative_rna_protein":"BLOCKED",
            "V04_source_validated_etiology_groups":"PARTIAL",
            "V12_event_coding_time_origin_censoring":"PARTIAL",
            "V13_external_validated_prediction":"BLOCKED",
        },
        "next_research_actions":[
            "Ng 2022: lawful patient-by-gene RNA quantitative matrix plus matching quantitative proteomics and matched distinct patient N",
            "Gao 2019: obtain verifiable OEP000321 processed RNA matrix and reconcile samples with PDC protein aliquots",
            "Xing 2023: confirm clinical RFS/OS 0/1 event coding, follow-up units, time origin and censoring from source codebook",
            "Xing 2023: confirm HBV/HCV marker positivity, etiology exclusivity and category mapping before group N or power claims",
        ],
        "interpretation":"These are source-file/metadata evidence levels. Even if every script executes, no same-patient quantitative omics join, true event N, or novelty is implied.",
        "patient_identifiers_or_raw_values_in_summary":False,
    }

def markdown(report):
    lines=["# HCC source-evidence depth (not execution status)","",
           f"- Source report files: {report['source_reports_read']}/{report['expected_source_reports']}",
           "- Research V02: BLOCKED; V04: PARTIAL; V12: PARTIAL; V13: BLOCKED.","",
           "| Investigation | Evidence depth |",
           "|---|---|"]
    for name,rec in report["evidence"].items():
        lines.append(f"| {name.replace('_',' ')} | {rec['state']} |")
    lines+=["","## Next unblock actions",""]
    for task in report["next_research_actions"]:
        lines.append(f"- {task}")
    lines+=["","**Do not equate CI execution PASS with verified RNA–protein patient pairing or an interpretable clinical endpoint.**",""]
    return "\n".join(lines)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--artifact-root",required=True)
    p.add_argument("--out-dir",default="results/hcc/parallel/batch")
    args=p.parse_args()
    report=summarize(args.artifact_root)
    out=pathlib.Path(args.out_dir)
    out.mkdir(parents=True,exist_ok=True)
    (out/"evidence_depth.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    (out/"evidence_depth.md").write_text(markdown(report),encoding="utf-8")
    print(markdown(report))
    # A missing evidence file must be visible as a CI failure, not silently a
    # claim of full source coverage; distinct from biological HOLD.
    if not report["all_source_reports_read"]:
        raise SystemExit(1)

if __name__=="__main__":
    main()
