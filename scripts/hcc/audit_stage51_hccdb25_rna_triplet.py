#!/usr/bin/env python3
"""Stage51: reproducible HCCDB25 RNA triplet and fail-closed lineage audit.

The third-party published code references three local TSV files, NOT public
download links. A consistent triplet must not be misclassified as independently
source-authorized or as joined to official PDC aliquot / tumor case identities.
Never print SAMPLE_ID, PATIENT_ID, expression, clinical fields, individual
records or unverified URLs. Store only aggregate counts and file checksums.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import pathlib
import re
import urllib.parse
import urllib.request

EXPECTED={
    "expression":"HCCDB25_mRNA_level3.txt",
    "sample":"HCCDB25.sample.txt",
    "patient":"HCCDB25.patient.txt",
}
CODE_URL=("https://raw.githubusercontent.com/Cthfudan/hcc/"
          "f3d3cbcd8aa046c0f4d4b6a57b33932aeb502ffe/"
          "Fig6%26S11/survival_analysis_multivariate.R")
CODE_SHA="f3d3cbcd8aa046c0f4d4b6a57b33932aeb502ffe"
MAX_SOURCE=100000
MAX_BYTES=1024*1024*1024
MAX_COLUMNS=10000
MAX_ROWS=100000
MAX_CLINICAL=10000

class AuditError(ValueError):
    pass

def sha256_file(file):
    digest=hashlib.sha256()
    count=0
    with pathlib.Path(file).open("rb") as fh:
        while True:
            block=fh.read(1024*1024)
            if not block:break
            digest.update(block)
            count+=len(block)
            if count>MAX_BYTES:raise AuditError("SOURCE_FILE_BYTES_OVER_LIMIT")
    return digest.hexdigest(),count

def make_rows(file):
    fh=pathlib.Path(file).open("r",encoding="utf-8-sig",newline="")
    return fh,csv.reader(fh,delimiter="\t")

def parse_expression(file):
    handle,reader=make_rows(file)
    try:
        header=next(reader,None)
        if not header:raise AuditError("EXPRESSION_HEADER_MISSING")
        labels=[h.strip() for h in header]
        if labels.count("Symbol")!=1:raise AuditError("EXPRESSION_SYMBOL_SCHEMA")
        gene_col=labels.index("Symbol")
        sample_ids=labels[gene_col+1:]
        if not (1<=len(sample_ids)<=MAX_COLUMNS):raise AuditError("EXPRESSION_SAMPLE_COLUMNS_INVALID")
        if not all(sample_ids) or len(sample_ids)!=len(set(sample_ids)):
            raise AuditError("EXPRESSION_DUPLICATE_OR_MISSING_SAMPLE_IDS")
        if any(x in ("SAMPLE_ID","PATIENT_ID","TYPE") for x in sample_ids):
            raise AuditError("EXPRESSION_METADATA_HEADER_LEAK")
        n_genes=0
        genes=set()
        missing_cells=0
        numeric_cells=0
        for row in reader:
            if not any(v.strip() for v in row):continue
            if len(row)!=len(labels):raise AuditError("EXPRESSION_ROW_WIDTH_MISMATCH")
            gene=row[gene_col].strip()
            if not gene or gene in genes:raise AuditError("EXPRESSION_GENE_SYMBOL_MISSING_OR_DUPLICATED")
            genes.add(gene)
            n_genes+=1
            if n_genes>MAX_ROWS:raise AuditError("EXPRESSION_TOO_MANY_GENES")
            for value in row[gene_col+1:]:
                v=value.strip()
                if v=="" or v.upper() in {"NA","NAN","NULL"}:
                    missing_cells+=1
                    continue
                try:
                    v=float(v)
                except ValueError as exc:raise AuditError("EXPRESSION_VALUE_NONNUMERIC") from exc
                if not math.isfinite(v):raise AuditError("EXPRESSION_VALUE_NONFINITE")
                numeric_cells+=1
        if n_genes<2 or numeric_cells<2:raise AuditError("EXPRESSION_TOO_FEW_FEATURES")
        return {"sample_ids":sample_ids,"genes_n":n_genes,
                "observed_numeric_expression_cells_n":numeric_cells,
                "missing_expression_cells_n":missing_cells,
                "discarded_gene_metadata_columns_n":gene_col}
    finally:
        handle.close()

def parse_sample_transposed(file):
    """HCCDB25.sample.txt is field-by-column, per published R t() pipeline."""
    fh,reader=make_rows(file)
    try:
        data={}
        width=None
        for row in reader:
            if not any(v.strip() for v in row):continue
            if width is None:width=len(row)
            if len(row)!=width or width>MAX_COLUMNS+1 or width<2:
                raise AuditError("SAMPLE_TRANSPOSED_WIDTH_INVALID")
            key=row[0].strip()
            if not key or key in data:raise AuditError("SAMPLE_DUPLICATE_OR_EMPTY_METADATA_KEYS")
            data[key]=[v.strip() for v in row[1:]]
            if len(data)>100:raise AuditError("SAMPLE_TOO_MANY_METADATA_KEYS")
        required={"SAMPLE_ID","PATIENT_ID","TYPE"}
        if not required.issubset(data):raise AuditError("SAMPLE_REQUIRED_METADATA_KEYS_MISSING")
        samples=data["SAMPLE_ID"]
        patients=data["PATIENT_ID"]
        types=data["TYPE"]
        if any(not v for v in samples) or len(set(samples))!=len(samples):
            raise AuditError("SAMPLE_DUPLICATE_OR_EMPTY_SAMPLE_ID")
        if any(not v for v in patients):
            raise AuditError("SAMPLE_EMPTY_PATIENT_ID")
        if any(not v for v in types):
            raise AuditError("SAMPLE_EMPTY_TYPE")
        return {"sample_ids":samples,"patient_ids":patients,"types":types,
                "sample_metadata_fields_n":len(data)}
    finally:
        fh.close()

def parse_patient(file):
    fh,reader=make_rows(file)
    try:
        header=next(reader,None)
        if not header or header.count("PATIENT_ID")!=1:
            raise AuditError("PATIENT_TABLE_PATIENT_ID_REQUIRED")
        if len(header)>MAX_CLINICAL or len(set(header))!=len(header):
            raise AuditError("PATIENT_TABLE_BAD_HEADER")
        idx=header.index("PATIENT_ID")
        ids=[]
        for row in reader:
            if not any(v.strip() for v in row):continue
            if len(row)!=len(header):raise AuditError("PATIENT_ROW_WIDTH_MISMATCH")
            value=row[idx].strip()
            if not value:raise AuditError("PATIENT_EMPTY_PATIENT_ID")
            ids.append(value)
            if len(ids)>MAX_ROWS:raise AuditError("PATIENT_TOO_MANY_ROWS")
        if len(ids)!=len(set(ids)):raise AuditError("PATIENT_DUPLICATE_PATIENT_ID")
        return {"patient_ids":ids,"clinical_variables_n":len(header)-1}
    finally:
        fh.close()

def evaluate_triplet(directory):
    folder=pathlib.Path(directory)
    result={"status":"TRIPLET_UNVERIFIED","files_present_n":0,
            "source_authorized_hccdb25_download":False,
            "verified_PDC_official_aliquot_bridge":False,
            "raw_patient_ids_or_expression_values_exported":False}
    try:
        files={k:folder/name for k,name in EXPECTED.items()}
        result["files_present_n"]=sum(p.is_file() for p in files.values())
        if result["files_present_n"]!=3:
            result["status"]="TRIPLET_SOURCE_FILES_MISSING"
            return result
        result["source_file_metadata"]={}
        for name,file in files.items():
            digest,size=sha256_file(file)
            result["source_file_metadata"][name]={"file_sha256":digest,"file_size_bytes":size}
        expr=parse_expression(files["expression"])
        sample=parse_sample_transposed(files["sample"])
        patient=parse_patient(files["patient"])
        expression_ids=set(expr["sample_ids"])
        samples=set(sample["sample_ids"])
        patients=set(patient["patient_ids"])
        tumor_idx=[i for i,v in enumerate(sample["types"]) if v.upper()=="HCC"]
        tumor_samples={sample["sample_ids"][i] for i in tumor_idx}
        tumor_patients={sample["patient_ids"][i] for i in tumor_idx}
        linked_patients=set(sample["patient_ids"])
        result["aggregate"]={
            "genes_n":expr["genes_n"],
            "expression_sample_columns_n":len(expr["sample_ids"]),
            "sample_metadata_records_n":len(sample["sample_ids"]),
            "sample_metadata_fields_n":sample["sample_metadata_fields_n"],
            "unique_sample_patient_ids_n":len(linked_patients),
            "tumor_HCC_sample_metadata_n":len(tumor_idx),
            "tumor_HCC_unique_patient_ids_n":len(tumor_patients),
            "patient_clinical_rows_n":len(patient["patient_ids"]),
            "patient_clinical_variables_n":patient["clinical_variables_n"],
            "expression_samples_missing_metadata_n":len(expression_ids-samples),
            "sample_metadata_missing_expression_n":len(samples-expression_ids),
            "sample_patient_ids_missing_clinical_n":len(linked_patients-patients),
            "tumor_samples_missing_expression_n":len(tumor_samples-expression_ids),
            "tumor_patients_missing_clinical_n":len(tumor_patients-patients),
            "observed_numeric_expression_cells_n":expr["observed_numeric_expression_cells_n"],
            "missing_expression_cells_n":expr["missing_expression_cells_n"],
        }
        full=(expression_ids==samples and linked_patients.issubset(patients) and len(tumor_idx)>0)
        result["status"]=("TRIPLET_INTERNAL_CONSISTENCY_PASS_SOURCE_HOLD" if full
                          else "TRIPLET_INTERNAL_CROSSWALK_HOLD")
        result["internal_RNA_patient_mapping_consistent"]=full
    except (ValueError,OSError,csv.Error,UnicodeError) as exc:
        result["status"]="TRIPLET_VALIDATION_HOLD"
        result["error_class"]=type(exc).__name__
        result["validation_error_code"]=str(exc) if isinstance(exc,AuditError) else "SOURCE_PARSE_ERROR"
    return result

def source_code_reference(opener=urllib.request.urlopen):
    item={"source_repository":"Cthfudan/hcc","pinned_source_commit":CODE_SHA,
          "status":"SOURCE_CODE_REFERENCE_HOLD",
          "actual_RNA_file_download_verified":False,
          "public_download_endpoint_verified":False}
    try:
        request=urllib.request.Request(CODE_URL,headers={
            "User-Agent":"HCC-Stage51-public-code-reference/1.0",
            "Range":f"bytes=0-{MAX_SOURCE-1}"})
        with opener(request,timeout=20) as response:
            if response.geturl()!=CODE_URL:raise ValueError("unverified redirect")
            body=response.read(MAX_SOURCE+1)
            if len(body)>MAX_SOURCE:raise AuditError("CODE_REFERENCE_TOO_LARGE")
            if response.status not in (200,206):raise AuditError("CODE_REFERENCE_BAD_HTTP")
        text=body.decode("utf-8")
        counts={k:len(re.findall(re.escape(name),text)) for k,name in EXPECTED.items()}
        item["reference_source_code_sha256"]=hashlib.sha256(body).hexdigest()
        item["documented_local_path_file_reference_counts"]=counts
        item["documents_OEP000321_accession_literal"]="OEP000321" in text
        complete=(all(v>=1 for v in counts.values()) and item["documents_OEP000321_accession_literal"])
        item["status"]="RESEARCH_CODE_LOCAL_TRIPLET_REFERENCES_VERIFIED" if complete else "REFERENCE_INCOMPLETE_HOLD"
    except Exception as exc:
        item["error_class"]=type(exc).__name__
    return item

def review(verified):
    if verified.get("source_authorized_hccdb25_download") is not False or (
        verified.get("verified_PDC_official_aliquot_bridge") is not False):
        raise ValueError("unsupported original-data or cross-cohort assertion")
    return {"schema":"hcc-stage51-triplet-postflight-v1",
            "local_triplet_internal_consistency_checked":verified.get("status")=="TRIPLET_INTERNAL_CONSISTENCY_PASS_SOURCE_HOLD",
            "original_RNA_source_authorized":False,
            "official_RNA_to_PDC_same_biopsy_join":False,
            "G30_B":"HOLD_PDC_OFFICIAL_ALIQUOT_JOIN",
            "G30_C":"HOLD_RNA_SOURCE_PROVENANCE_OR_ORIGINAL_TRIPLET",
            "G30_D":"STOP_UNVERIFIED_PATIENT_MATCHED_RNA_PROTEIN",
            "biological_claims_allowed":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--data-dir",default=None)
    p.add_argument("--out",required=True)
    p.add_argument("--reference-only",action="store_true")
    args=p.parse_args()
    if args.reference_only:
        source=source_code_reference()
        report={"schema":"hcc-stage51-research-code-reference-v1",
                "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
                "reference":source}
        print("HCC_STAGE51_REFERENCE",json.dumps(report))
    else:
        source=(evaluate_triplet(args.data_dir) if args.data_dir else
                {"status":"TRIPLET_SOURCE_FILES_NOT_PROVIDED",
                 "source_authorized_hccdb25_download":False,
                 "verified_PDC_official_aliquot_bridge":False})
        report={"schema":"hcc-stage51-rna-triplet-identity-audit-v1",
                "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
                "triplet":source,"review":review(source)}
        print("HCC_STAGE51_RNA",json.dumps(report))
    target=pathlib.Path(args.out);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")

if __name__=="__main__":main()
