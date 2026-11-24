#!/usr/bin/env python3
"""Audit historic GTEx_V8-labeled QTL against EBI-imported *original* GTEx V8.

Inputs: read-only 6,242 archived gene-SNP rows frozen from a Google Drive file
dated 2026-10-11, whose original SHA256 was independently checked. The source
provenance decision is limited to those 6,242 exact gene+GRCh38+REF/ALT keys:
matching beta & p implies matching *statistics*, not identical complete source
bytes or complete original v8 cis universe.

Downloads no requester-pays objects. Two bounded indexed HTTPS queries using
the existing official original-imported GTEx V8 gene-expression source and
predeclared full cis comparison windows.
"""
import argparse
from collections import Counter
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import statistics

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
spec=importlib.util.spec_from_file_location("gtex_legacy_fullcis",HERE/"audit_gtex_v8_qtd_four_gene_fullcis.py")
F=importlib.util.module_from_spec(spec)
spec.loader.exec_module(F)
BASE=F.BASE

FROZEN={
 "4":("research/is/evidence/is_legacy_gtex_v8_labeled_allele_lineage_chr4_20261011.tsv",
      "ddace217c5b92d0503ab9261fef65c26ac3df9dd33a83276c79c75320ea527a2"),
 "10":("research/is/evidence/is_legacy_gtex_v8_labeled_allele_lineage_chr10_20261011.tsv",
       "0ca7fe24369a9506a4e8dae76dd26d6da4565ee67d2012f03b082eb7608a52a5"),
}
ORIGINAL_DRIVE_ID="1d-Pnz3NPQ4_pwgKy-ZbtCph7P0OIatHd"
ORIGINAL_DRIVE_SHA256="e2eda6ee41a8bd4d2437df24f8a64485f0ab8762bdf0649036c75327f59edffd"
ARCHIVE_EXPECTED={"FGF5":1568,"C4orf22":1568,"CALHM2":1553,"NEURL":1553}
COLUMNS=["gene_id","chromosome","position_GRCh38","ref","alt",
        "archive_source_p","archive_source_beta","archive_source_maf",
        "archive_source_ac","archive_source_an","archive_QTD_p","archive_QTD_beta"]
MAX_EXAMPLES=25

def load_legacy(root=ROOT):
    answer={}
    counts=Counter()
    id_to_name={v:k for k,v in BASE.TARGET_GENES.items()}
    for chrom,(rel,expected_sha) in FROZEN.items():
        source=Path(root)/rel
        if hashlib.sha256(source.read_bytes()).hexdigest()!=expected_sha:
            raise ValueError("Historic GTEx source lineage frozen panel SHA256 failed")
        with source.open("rt",encoding="utf-8",newline="") as f:
            reader=csv.DictReader(f,delimiter="\t")
            if reader.fieldnames!=COLUMNS:
                raise ValueError("Historic QTL lineage columns changed")
            for row in reader:
                gene=row["gene_id"]
                if gene not in id_to_name or row["chromosome"]!=chrom:
                    raise ValueError("Historic QTL gene/chr mismatches")
                pos=int(row["position_GRCh38"])
                ref=row["ref"].upper();alt=row["alt"].upper()
                if not pos>0 or not ref or not alt or ref==alt or not set(ref+alt)<=set("ACGT"):
                    raise ValueError("Invalid legacy QTL key")
                p,b,maf,ac,an=map(float,[row[x] for x in
                    ("archive_source_p","archive_source_beta","archive_source_maf",
                     "archive_source_ac","archive_source_an")])
                qp,qbeta=map(float,[row[x] for x in ("archive_QTD_p","archive_QTD_beta")])
                if not all(map(math.isfinite,[p,b,maf,ac,an,qp,qbeta])) or not(
                    0<=p<=1 and 0<=qp<=1 and 0<=maf<=.5 and 0<=ac<=an and an>0):
                    raise ValueError("Invalid archived QTL statistics or QC attributes")
                key=(gene,chrom,pos,ref,alt)
                if key in answer:raise ValueError("Duplicate legacy gene-variant row")
                answer[key]={"p":p,"beta":b,"maf":maf,"ac":ac,"an":an,
                             "archive_qtd_p":qp,"archive_qtd_beta":qbeta}
                counts[id_to_name[gene]]+=1
    if dict(counts)!=ARCHIVE_EXPECTED or len(answer)!=6242:
        raise ValueError("Historic complete 6242 gene-variant audit missing")
    return answer

def p_close(a,b):
    return a==b if a<=0 or b<=0 else abs(math.log10(a)-math.log10(b))<=0.0001

def stat_close(a,b):
    return math.isclose(a,b,rel_tol=.00005,abs_tol=.000001)

def compare(archived,live,qc):
    src_keys=set(archived)
    live_keys=set(live)
    counters={g:Counter() for g in ARCHIVE_EXPECTED}
    mismatch_counts=Counter()
    mismatches=[]
    d_beta=[]
    d_logp=[]
    for key in sorted(src_keys):
        gene=next((g for g,v in BASE.TARGET_GENES.items() if v==key[0]),None)
        if gene is None:raise ValueError("Unrecognized gene in archived records")
        metrics=counters[gene]
        metrics["archive_rows"]+=1
        if key not in live:
            metrics["missing_official_original"]+=1
            mismatch_counts["missing_official_original"]+=1
            if len(mismatches)<MAX_EXAMPLES:
                mismatches.append({"gene":gene,"key":key,"reason":["NOT_IN_OFFICIAL_ORIGINAL"]})
            continue
        metrics["exact_key_matched"]+=1
        ob,se,op=live[key]
        a=archived[key]
        b_ok=stat_close(a["beta"],ob)
        p_ok=p_close(a["p"],op)
        if b_ok:metrics["beta_pass"]+=1
        else: mismatch_counts["beta_fail"]+=1
        if p_ok:metrics["p_pass"]+=1
        else: mismatch_counts["p_fail"]+=1
        if b_ok and p_ok:metrics["joint_numeric_pass"]+=1
        else:
            metrics["joint_numeric_fail"]+=1
            if len(mismatches)<MAX_EXAMPLES:
                mismatches.append({"gene":gene,"key":key,
                   "reason":[x for x,ok in (("beta",b_ok),("p",p_ok)) if not ok],
                   "archive_source_beta":a["beta"],"official_original_beta":ob,
                   "archive_source_p":a["p"],"official_original_p":op})
        d_beta.append(abs(a["beta"]-ob))
        if 0<a["p"]<=1 and 0<op<=1:
            d_logp.append(abs(math.log10(a["p"])-math.log10(op)))
        f=qc.get(key,{})
        maf=f.get("maf")
        ac=f.get("ac")
        an=f.get("an")
        if maf is None:
            metrics["official_maf_missing"]+=1
        elif math.isclose(maf,a["maf"],rel_tol=0,abs_tol=2.1e-6):
            metrics["maf_pass"]+=1
        else:
            metrics["maf_fail"]+=1
        if an is None:
            metrics["official_an_missing"]+=1
        elif math.isclose(an,a["an"],rel_tol=0,abs_tol=.001):
            metrics["an_pass"]+=1
        else:
            metrics["an_fail"]+=1
        if ac is None or an is None:
            metrics["official_ac_missing"]+=1
        elif math.isclose(ac,a["ac"],abs_tol=.001,rel_tol=0):
            metrics["ac_same"]+=1
        elif math.isclose(an-ac,a["ac"],abs_tol=.001,rel_tol=0):
            metrics["ac_complement"]+=1
        else:
            metrics["ac_other"]+=1
    full={gene:dict(v) for gene,v in counters.items()}
    count_joint=sum(c["joint_numeric_pass"] for c in counters.values())
    identity=count_joint==len(archived)==6242 and not mismatch_counts
    return {"archived_gene_variant_rows":len(archived),
        "official_source_target_variant_rows":len(live),
        "exact_gene_variant_rows_observed":sum(v["exact_key_matched"] for v in counters.values()),
        "archived_source_p_and_beta_match_official_original":count_joint,
        "archived_source_match_rate":count_joint/len(archived) if archived else 0,
        "median_absolute_beta_diff":statistics.median(d_beta) if d_beta else None,
        "median_absolute_log10p_diff":statistics.median(d_logp) if d_logp else None,
        "by_gene":full,"difference_classes":dict(mismatch_counts),
        "bounded_discordant_examples":mismatches,
        "legacy_QTL_statistical_source_identity_gate":(
            "6242_ARCHIVED_GTEX_V8_QTL_STATS_MATCH_EBI_IMPORTED_ORIGINAL"
            if identity else "HOLD_ARCHIVED_GTEX_V8_SOURCE_NOT_PROVEN_MATCH_TO_IMPORTED_ORIGINAL"),
        "does_not_authenticate_original_historic_source_file_bytes":True,
        "does_not_authenticate_full_original_GTEx_cis_universe":True,
        "no_allele_beta_flipping":True,
        "MR_or_coloc_performed":False,
        "scientific_status":"HOLD",
    }

def run(backend=None,header_provider=None,wait=2.5,archive_root=ROOT):
    result={"schema":"IS_Archived_GTExV8_6242_vs_Official_Imported_GTExv8_v1",
        "archive_GDrive_original_file_id":ORIGINAL_DRIVE_ID,
        "archive_GDrive_original_SHA256":ORIGINAL_DRIVE_SHA256,
        "pinned_ARCHIVED_CHR_panels_SHA256":{k:v[1] for k,v in FROZEN.items()},
        "live_official_ebi_imported_original_GTEx_v8":BASE.ORIGINAL,
        "gene_samples_are_not_independent_replications":True,
        "all_2225_gene_screen_unaltered":True,
        "formal_MR_or_coloc_performed":False,"scientific_status":"HOLD"}
    try:
        historic=load_legacy(archive_root)
    except (OSError,ValueError,csv.Error) as exc:
        result["legacy_source_gate"]="HOLD_ARCHIVED_SOURCE_INVALID"
        result["archive_error_class"]=type(exc).__name__
        return result
    if backend is None:
        os.environ["CURL_CA_BUNDLE"]="/etc/ssl/certs/ca-certificates.crt"
        os.environ["SSL_CERT_FILE"]="/etc/ssl/certs/ca-certificates.crt"
        import pysam
        backend=pysam
    handle=None
    try:
        handle=backend.TabixFile(BASE.ORIGINAL)
        header=(header_provider(BASE.ORIGINAL)["header"] if header_provider
                else BASE.original_bgzf_first_line(BASE.ORIGINAL)["header"])
        live={}
        qc={}
        scans={}
        from time import sleep
        for region in F.FULL_WINDOWS:
            if wait:sleep(wait)
            data=F.read_source_rows(handle,header,region)
            chrom=region[0]
            scans[chrom]={k:v for k,v in data.items()
                           if k not in ("rows","genotype_qc_by_key","invalid_source_rows_bounded")}
            if data["status"]!="SCREEN_VALID":
                result["legacy_source_gate"]="HOLD_LIVE_ORIGINAL_SOURCE_WINDOW"
                result["window_audit"]=scans
                return result
            for key,value in data["rows"].items():
                if key in live:raise ValueError("Source window duplicate genomic key")
                live[key]=value
            qc.update(data.get("genotype_qc_by_key",{}))
        result["window_audit"]=scans
        result["legacy_source_comparison"]=compare(historic,live,qc)
        result["legacy_source_gate"]=result["legacy_source_comparison"]["legacy_QTL_statistical_source_identity_gate"]
    except (OSError,ValueError,KeyError,TypeError) as exc:
        result["legacy_source_gate"]="HOLD_SOURCE_NETWORK_OR_SCHEMA_FAILURE"
        result["source_failure_class"]=type(exc).__name__
        result["bounded_error"]=str(exc)[-200:]
    finally:
        if handle is not None:handle.close()
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    x=run()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(x,sort_keys=True,indent=2)+"\n")
    print("IS_LEGACY_GTEX_ORIGINAL="+json.dumps(x,sort_keys=True,separators=(",",":")))
if __name__=="__main__":
    main()
