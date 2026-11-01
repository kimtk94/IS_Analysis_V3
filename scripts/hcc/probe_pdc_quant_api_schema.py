#!/usr/bin/env python3
"""Bounded PDC quantDataMatrix probe: gene-by-aliquot schema only.

PDC's public quantDataMatrix interface advertises a gene/aliquot matrix.  Probe
the beginning of the JSON response (<= 768 KiB), inspect structural keys, and
compare candidate header keys to official PDC aliquot ids/submitter ids only in
memory. Patient/aliquot identifiers, expression values and full matrices are
never written to a public CI artifact. A first-row key match alone is only a
schema check: it does NOT validate all gene rows or any RNA expression join.
"""
import argparse
import datetime as dt
import json
import pathlib
import re
import urllib.error
import urllib.request
from probe_pdc_biospecimen_lineage import get_cases

ENDPOINT = "https://pdc.cancer.gov/graphql"
STUDIES = ("PDC000198", "PDC000199")
MAX_PREFIX = 768 * 1024


def request_prefix(study, max_prefix=MAX_PREFIX):
    query = '{ quantDataMatrix(pdc_study_id:"' + study + '", data_type:"log2_ratio") }'
    request = urllib.request.Request(
        ENDPOINT, data=json.dumps({"query":query}).encode("utf-8"),
        method="POST", headers={"Content-Type":"application/json",
                                "Accept":"application/json",
                                "Accept-Encoding":"identity",
                                "User-Agent":"HCC-Public-Quant-Matrix-Schema/1.0"})
    with urllib.request.urlopen(request, timeout=25) as response:
        code = response.status
        length = response.headers.get("Content-Length")
        raw = response.read(max_prefix + 1)
    try:
        length = int(length) if length else None
    except ValueError:
        length = None
    return raw[:max_prefix], {"http_status":code,
                              "reported_content_length":length,
                              "prefix_limit_reached":len(raw)>max_prefix,
                              "prefix_bytes_read":min(len(raw),max_prefix)}


def discover_keys(prefix):
    """Find JSON keys only from first row. No values or identifying keys emitted."""
    # Do not trust later gene rows (only the first array object for column keys).
    start = prefix.find(b'"quantDataMatrix"')
    if start < 0:
        return {"status":"MATRIX_FIELD_ABSENT","keys":[]}
    bracket = prefix.find(b"[",start)
    if bracket<0:
        return {"status":"JSON_LIST_START_ABSENT","keys":[]}
    brace=prefix.find(b"{",bracket)
    if brace<0:
        return {"status":"FIRST_RECORD_OBJECT_ABSENT","keys":[]}
    # An aliquot field has an identifier-like name but may include hyphens.
    # The first row can be larger than max_prefix. This only inventories keys in prefix.
    prefix_row=prefix[brace:]
    cutoff=prefix_row.find(b"}")
    status="FIRST_RECORD_PREFIX_ONLY" if cutoff<0 else "FIRST_RECORD_KEYS_INSPECTED"
    if cutoff>=0:
        prefix_row=prefix_row[:cutoff]
    names=[]
    for m in re.finditer(rb'"((?:[^"\\]|\\.){1,200})"\s*:', prefix_row):
        try:
            key=json.loads(b'"'+m.group(1)+b'"')
        except Exception:
            continue
        names.append(key)
    return {"status":status,"keys":list(dict.fromkeys(names))}


def audit_field_keys(keys, biospecimens):
    aliq_ids={}
    for case in biospecimens:
        for sample in case.get("samples") or []:
            sample_type=sample.get("sample_type") or "Not Reported"
            for aliq in sample.get("aliquots") or []:
                for ident in (aliq.get("aliquot_id"),aliq.get("aliquot_submitter_id")):
                    if ident:
                        aliq_ids[str(ident)]=sample_type
    matches=[key for key in keys if key in aliq_ids]
    types={}
    for m in matches:
        k=aliq_ids[m]
        types[k]=types.get(k,0)+1
    return {"first_record_number_of_keys_inspected":len(keys),
            "first_record_matched_unique_aliquot_fields":len(set(matches)),
            "matched_aliquot_sample_type_distribution":types,
            "biospecimens_distinct_aliquot_aliases_for_comparison":len(aliq_ids),
            "gene_by_aliquot_quantitative_values_verified":False,
            "RNA_PROTEIN_SAMPLE_JOIN_VERIFIED":False}


def audit():
    output={}
    for study in STUDIES:
        record={"status":"UNVERIFIED","pdc_study_id":study}
        try:
            prefix,http=request_prefix(study)
            record.update(http)
            scanned=discover_keys(prefix)
            record["status"]=scanned["status"]
            if scanned["keys"]:
                try:
                    bios,_=get_cases(study)
                    record["header_mapping"]=audit_field_keys(scanned["keys"],bios)
                except Exception as exc:
                    record["biospecimen_compare_error_type"]=type(exc).__name__
                    record["header_key_count_inspected"]=len(scanned["keys"])
            else:
                record["header_key_count_inspected"]=0
                # Do not persist GraphQL error text, which might contain identifiers.
                record["possibly_graphql_error"]=b'"errors"' in prefix[:4096]
        except urllib.error.HTTPError as exc:
            record["status"]="HTTP_ERROR"
            record["http_status"]=exc.code
        except Exception as exc:
            record["status"]="UNVERIFIED"
            record["error_type"]=type(exc).__name__
        output[study]=record
    return output


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="results/hcc/feasibility/pdc_quant_api_schema.json")
    args=p.parse_args()
    report={"checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "official_api":ENDPOINT,"candidate_data_type":"log2_ratio",
            "source":"PDC documented quantDataMatrix, https://pdc.cancer.gov/pdc/publicapi-documentation",
            "matrix_prefix_probe":audit(),
            "claims":{"full_proteome_matrix_join":False,"rna_to_protein_join":False},
            "warning":"Bounded first JSON row/schema is not full expression quantification or clinical inclusion validation."}
    path=pathlib.Path(args.out)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
