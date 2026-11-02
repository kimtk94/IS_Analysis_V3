#!/usr/bin/env python3
"""Audit GTEx V8 actual collapsed gene model TSS for two original-QTL gaps.

Uses ONLY the GTEx consortium public reference gene model, not the
3.4GB requester-pays original allpairs associations. No tokens or billing.
If the public object is unavailable or needs requester-pays, report HOLD.
Gene model can contain transcript/exon-level features; use collapsed exon
union to identify strand-specific TSS only if well-formed and unique.
"""
import argparse
import hashlib
import json
import re
import tempfile
import urllib.request
from pathlib import Path

URL="https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf"
MAX_SOURCE_BYTES=175_000_000
GENES={
    "FGF5":("ENSG00000138675","4",81266618,"T","G"),
    "NEURL":("ENSG00000107954","10",102493750,"T","C"),
}
ATTR_GENE=re.compile(r'(?:^|;\s*)gene_id\s+"([^"]+)"')
CIS_RADIUS=1_000_000
HISTORICAL_TSS={"FGF5":80266599,"NEURL":103493979}

def inspect_gtf(path):
    path=Path(path)
    count=0
    hash=hashlib.sha256()
    gene_track={}
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):
            count+=len(chunk)
            if count>MAX_SOURCE_BYTES:
                raise ValueError("Public GTEx collapsed gene model exceeds hard source bound")
            hash.update(chunk)
    with path.open("rt",encoding="utf-8",newline="") as f:
        for line in f:
            if line.startswith("#"):continue
            cols=line.rstrip("\n").split("\t")
            if len(cols)!=9:continue
            found=ATTR_GENE.search(cols[8])
            if not found:continue
            stable=found.group(1).split(".")[0]
            if stable not in {x[0] for x in GENES.values()}:continue
            chrom=cols[0].removeprefix("chr")
            strand=cols[6]
            start,end=int(cols[3]),int(cols[4])
            if strand not in ("+","-") or start<1 or end<start:
                raise ValueError("Invalid GTEx v8 collapsed gene model geometry")
            current=gene_track.setdefault(stable,{
                "chromosome":chrom,"strand":strand,"min_start":start,"max_end":end,
                "feature_counts":{},"gene_model_versions":[]})
            if current["chromosome"]!=chrom or current["strand"]!=strand:
                raise ValueError("Conflicting source gene features")
            current["min_start"]=min(start,current["min_start"])
            current["max_end"]=max(end,current["max_end"])
            current["feature_counts"][cols[2]]=current["feature_counts"].get(cols[2],0)+1
            if len(current["gene_model_versions"])<3 and found.group(1) not in current["gene_model_versions"]:
                current["gene_model_versions"].append(found.group(1))
    if set(gene_track)!={x[0] for x in GENES.values()}:
        raise ValueError("Missing target genes in original GTEx collapsed model")
    details={}
    for name,(gene,chrom,variant_pos,ref,alt) in GENES.items():
        entry=gene_track[gene]
        if entry["chromosome"]!=chrom:
            raise ValueError("GTEx collapsed gene contig mismatch")
        tss=entry["min_start"] if entry["strand"]=="+" else entry["max_end"]
        lower=tss-CIS_RADIUS
        upper=tss+CIS_RADIUS
        bp_before=lower-variant_pos if variant_pos<lower else 0
        bp_after=variant_pos-upper if variant_pos>upper else 0
        details[name]={
            **entry,"stable_gene_id":gene,
            "GTEx_collapsed_TSS_GRCh38":tss,
            "GENCODE_v26_regular_TSS_GRCh38":HISTORICAL_TSS[name],
            "collapsed_vs_regular_TSS_delta_bp":tss-HISTORICAL_TSS[name],
            "GTEx_v8_cis_lower_TSS_minus1Mb":lower,
            "GTEx_v8_cis_upper_TSS_plus1Mb":upper,
            "QTD_only_SNP_GRCh38":f"{chrom}:{variant_pos}:{ref}:{alt}",
            "QTD_only_position_inside_collapsed_cis":lower<=variant_pos<=upper,
            "QTD_only_position_bp_before_lower":bp_before,
            "QTD_only_position_bp_after_upper":bp_after,
            "GTEx_original_source_nominal_association_found":False,
            "specific_variant_QC_filters_independently_proven":False,
        }
    outside=all(not d["QTD_only_position_inside_collapsed_cis"] for d in details.values())
    return {
        "schema":"GTExV8_actual_collapsed_gene_model_QTD_two_unmatched_QTL_variants_v1",
        "source_url":URL,"source_file_size_bytes":count,
        "whole_collapsed_model_SHA256":hash.hexdigest(),
        "model":"Original GTEx v8 public v26 collapsed genes GTF",
        "method":"GTEx v8 original FastQTL TSS ±1Mb",
        "two_target_SNPs_outside_actual_GTEx_v8_collapsed_cis_window":outside,
        "source_model_QTL_boundary_gate":(
            "ACTUAL_GTEXV8_COLLAPSED_GENE_MODEL_SUPPORTS_TWO_CIS_EXCLUSIONS"
            if outside else "HOLD_COLLAPSED_TSS_DISAGREES_WITH_QTL_EXCLUSION_HYPOTHESIS"),
        "genes":details,
        "not_evidence_of_independent_QTL_replication":True,
        "exact_full_GTEx_original_file_MD5_recomputed":False,
        "QTL_GWAS_LD_N_approved":False,
        "coloc_performed":False,"scientific_status":"HOLD",
    }

def inspect_download():
    req=urllib.request.Request(URL,headers={"User-Agent":"IS-GTEx-v8-Public-Collapsed-Model-Cis/1.0"})
    with tempfile.TemporaryDirectory(prefix="gtexv8gene_") as folder:
        dest=Path(folder)/"collapsed.gtf"
        total=0
        with urllib.request.urlopen(req,timeout=50) as response,dest.open("wb") as out:
            if response.url!=URL:
                raise ValueError("Unexpected GTEx official public model redirect")
            while True:
                chunk=response.read(1024*1024)
                if not chunk:break
                total+=len(chunk)
                if total>MAX_SOURCE_BYTES:
                    raise ValueError("Public GTEx model exceeds byte budget")
                out.write(chunk)
        return inspect_gtf(dest)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--local-gtf",type=Path,default=None)
    args=p.parse_args()
    try:
        result=inspect_gtf(args.local_gtf) if args.local_gtf else inspect_download()
    except (OSError,ValueError,UnicodeError) as exc:
        result={"schema":"GTExV8_actual_collapsed_gene_model_QTD_two_unmatched_QTL_variants_v1",
                "source_model_QTL_boundary_gate":"HOLD_PUBLIC_MODEL_SOURCE_UNAVAILABLE",
                "failure_class":type(exc).__name__,
                "error_bounded":str(exc)[-160:],
                "scientific_status":"HOLD","coloc_performed":False}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("IS_GTEX_COLLAPSED_CIS="+json.dumps(result,sort_keys=True,separators=(",",":")))

if __name__=="__main__":
    main()
