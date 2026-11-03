#!/usr/bin/env python3
"""Compare 20,689 archived QTD000161 summary statistics with live official QTD.

The archive was downloaded from connected Google Drive and independently
verified against 2026-10-11 archive SHA-256 before its data were frozen into
two tracked keyed TSV files. This is an ARCHIVED RECONSTRUCTED TSV identity
test, not full gzip-byte or EBI release-6 provenance certification.

No server SSH, paid downloads, or coloc. Fails closed on missing variants,
key collisions, allele differences, beta/SE/p mismatches or source read errors.
"""
import argparse
from collections import Counter
import csv
import importlib.util
import json
import math
import os
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SPEC=importlib.util.spec_from_file_location("qtd_fullcis_reader",HERE/"audit_gtex_v8_qtd_four_gene_fullcis.py")
QTL=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(QTL)
BASE=QTL.BASE

SOURCE_SHA256_ARCHIVED_20261011="a80d5a6903aea551d5ee6f8bebf12772aea411a7cd29437bd545d5d2452768a3"
ARCHIVE={
 "4":"research/is/evidence/is_priority4_archived_qtd000161_chr4_source20261011.tsv",
 "10":"research/is/evidence/is_priority4_archived_qtd000161_chr10_source20261011.tsv",
}
EXPECTED={"FGF5":5511,"C4orf22":5490,"CALHM2":4868,"NEURL":4820}
COLUMNS=["gene_id","chromosome","position_GRCh38","REF_GRCh38","ALT_GRCh38",
         "QTD000161_nominal_QTL_p","QTD000161_beta","QTD000161_se"]
MAX_MISMATCH_EXAMPLES=24

def load_archived(root=ROOT):
    by_key={}
    count=Counter()
    allowed={v:k for k,v in BASE.TARGET_GENES.items()}
    for chrom,relative in ARCHIVE.items():
        path=Path(root)/relative
        with path.open("rt",encoding="utf-8",newline="") as f:
            reader=csv.DictReader(f,delimiter="\t")
            if list(reader.fieldnames or [])!=COLUMNS:
                raise ValueError("Archived TSV frozen columns changed")
            for item in reader:
                gene=item["gene_id"]
                if gene not in allowed or item["chromosome"]!=chrom:
                    raise ValueError("Archive gene or chromosome is invalid")
                pos=int(item["position_GRCh38"])
                ref=item["REF_GRCh38"];alt=item["ALT_GRCh38"]
                if (pos<1 or not ref or not alt or ref==alt or
                    not set(ref+alt).issubset(set("ACGT"))):
                    raise ValueError("Archive variant allele/coordinate invalid")
                key=(gene,chrom,pos,ref,alt)
                if key in by_key:raise ValueError("Duplicate historical variant–gene")
                p,b,se=(float(item[k]) for k in (
                     "QTD000161_nominal_QTL_p","QTD000161_beta","QTD000161_se"))
                if not all(map(math.isfinite,(p,b,se))) or not 0<=p<=1 or se<=0:
                    raise ValueError("Invalid historical p/beta/se")
                by_key[key]=(b,se,p)
                count[allowed[gene]]+=1
    if dict(count)!=EXPECTED or len(by_key)!=20689:
        raise ValueError("Archive gene-specific fullcis denominator no longer matches")
    return by_key

def compare_archived_with_live(archived,live, *, capture=MAX_MISMATCH_EXAMPLES):
    original=set(archived)
    current=set(live)
    matched=original&current
    summary=Counter()
    by_gene={k:Counter() for k in EXPECTED}
    examples=[]
    for name in EXPECTED:
        gene=BASE.TARGET_GENES[name]
        by_gene[name]["historic_records"]=sum(k[0]==gene for k in original)
        by_gene[name]["live_records"]=sum(k[0]==gene for k in current)
        by_gene[name]["matched_variant_keys"]=sum(k[0]==gene for k in matched)
    for key in sorted(original|current):
        name=next((k for k,v in BASE.TARGET_GENES.items() if v==key[0]),None)
        if name is None:raise ValueError("Unexpected noncandidate gene in the same-source comparison")
        if key not in archived:
            by_gene[name]["live_only"]+=1
            summary["live_only"]+=1
            if len(examples)<capture:
                examples.append({"key":key,"reason":["live_only"]})
        elif key not in live:
            by_gene[name]["archive_only"]+=1
            summary["archive_only"]+=1
            if len(examples)<capture:
                examples.append({"key":key,"reason":["archive_only"]})
        else:
            b1,se1,p1=archived[key]
            b2,se2,p2=live[key]
            beta_ok=math.isclose(b1,b2,rel_tol=0.00005,abs_tol=0.000001)
            se_ok=math.isclose(se1,se2,rel_tol=0.00005,abs_tol=0.000001)
            p_ok=(p1==p2 if p1<=0 or p2<=0 else
                 abs(math.log10(p1)-math.log10(p2))<=0.0001)
            if beta_ok and se_ok and p_ok:
                by_gene[name]["numerical_pass"]+=1
                summary["numerical_pass"]+=1
            else:
                by_gene[name]["numerical_fail"]+=1
                summary["numerical_fail"]+=1
                if not beta_ok:summary["beta_discordance"]+=1
                if not se_ok:summary["se_discordance"]+=1
                if not p_ok:summary["p_discordance"]+=1
                if len(examples)<capture:
                    examples.append({"key":key,
                       "reason":[n for n,valid in (
                          ("beta_ALT",beta_ok),("se",se_ok),("p",p_ok)) if not valid],
                       "historical":[b1,se1,p1],"live":[b2,se2,p2]})
    all_good=(summary["numerical_pass"]==len(archived)==len(live)==20689
              and summary["numerical_fail"]==0
              and summary["archive_only"]==0 and summary["live_only"]==0)
    return {
        "archived_gene_variant_pairs":len(archived),
        "live_gene_variant_pairs":len(live),
        "variant_keys_in_both":len(matched),
        "numerically_concordant_beta_se_p":summary["numerical_pass"],
        "by_gene":{name:dict(values) for name,values in by_gene.items()},
        "difference_classes":dict(summary),
        "bounded_difference_examples":examples,
        "archived_to_current_QTD_stats_gate":(
            "20689_HISTORICAL_RECONSTRUCTED_QTL_VARIANT_STATISTICS_MATCH_LIVE_QTD"
            if all_good else "HOLD_HISTORICAL_QTD_STATISTICS_CHANGED_OR_COVERAGE_INCOMPLETE"),
        "original_2023_EBI_release6_gzip_sha256_authenticated":False,
        "original_GTEx_v8_FastQTL_equivalence_verified":False,
        "coloc_performed":False,
        "scientific_status":"HOLD",
    }

def live_audit(backend=None,header_provider=None,wait=QTL.WAIT_SECONDS,archive_root=ROOT):
    outcome={"schema":"IS_QTD000161_20261011_archive_vs_20261103_official_source_v1",
        "frozen_archived_GDrive_source_SHA256":SOURCE_SHA256_ARCHIVED_20261011,
        "frozen_archived_source_GDrive_ID":"12g0HrusD3ONV2v9ormnCDMfLNGdXqibU",
        "official_live_source":BASE.REPROCESSED,
        "historic_file_byte_identity_to_manifest_confirmed_locally":True,
        "historic_raw_EBI_QTD_gzip_hash_known":False,
        "original_GTEx_model_authenticated":False,
        "coloc_performed":False,"scientific_status":"HOLD"}
    try:
        archived=load_archived(archive_root)
    except (OSError,ValueError,csv.Error) as exc:
        outcome["archive_to_live_gate"]="HOLD_ARCHIVE_SOURCE_INVALID"
        outcome["archive_failure_class"]=type(exc).__name__
        return outcome
    if backend is None:
        os.environ["CURL_CA_BUNDLE"]="/etc/ssl/certs/ca-certificates.crt"
        os.environ["SSL_CERT_FILE"]="/etc/ssl/certs/ca-certificates.crt"
        import pysam
        backend=pysam
    handle=None
    try:
        handle=backend.TabixFile(BASE.REPROCESSED)
        header=(header_provider(BASE.REPROCESSED)["header"] if header_provider
                else BASE.original_bgzf_first_line(BASE.REPROCESSED)["header"])
        rows={}
        source_scan={}
        from time import sleep
        for region in QTL.FULL_WINDOWS:
            if wait:sleep(wait)
            current=QTL.read_source_rows(handle,header,region)
            chromosome=region[0]
            source_scan[chromosome]={k:v for k,v in current.items()
                          if k not in ("rows","genotype_qc_by_key","invalid_source_rows_bounded")}
            if current["status"]!="SCREEN_VALID":
                outcome["archive_to_live_gate"]="HOLD_QTD_REMOTE_SCREEN_FAILED"
                outcome["source_scans"]=source_scan
                return outcome
            for key,v in current["rows"].items():
                if key in rows:raise ValueError("Variant is duplicated across source intervals")
                rows[key]=v
        outcome["source_scans"]=source_scan
        outcome["live_statistic_comparison"]=compare_archived_with_live(archived,rows)
        outcome["archive_to_live_gate"]=outcome["live_statistic_comparison"]["archived_to_current_QTD_stats_gate"]
    except (OSError,ValueError,KeyError,TypeError) as exc:
        outcome["archive_to_live_gate"]="HOLD_LIVE_SOURCE_OR_STATS_QUERY_FAILURE"
        outcome["source_failure_class"]=type(exc).__name__
        outcome["source_failure_detail_bounded"]=str(exc)[-160:]
    finally:
        if handle is not None:handle.close()
    return outcome

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    report=live_audit()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("IS_QTD_ARCHIVED_VS_LIVE="+json.dumps(report,sort_keys=True,separators=(",",":")))
if __name__=="__main__":
    main()
