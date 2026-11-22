#!/usr/bin/env python3
"""One indexed chr10 survey per official eQTL Catalogue source; no raw SNP uploads.

A 2.05Mb P3A positional window is NOT necessarily each gene's full cis domain.
REF/ALT and signed source beta/SE are read, but build-matched GWAS/N are HOLD.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib
import io
import json
import math
import subprocess
import tempfile
import time
from collections import Counter
from pathlib import Path

QTL=importlib.import_module("scripts.is.audit_official_p3a_qtl_catalogue_provenance")
UCSC=importlib.import_module("scripts.is.audit_p3a_chr10_ucsc_liftover_fasta")
MAX_LINES=2_000_000
MAX_BYTES=260_000_000
COLUMNS=("chromosome","position","ref","alt","variant","gene_id",
         "molecular_trait_id","beta","se","maf","pvalue")
EXPECTED_N={"QTD000131":387,"QTD000171":205,"QTD000434":560}
OFFICIAL_COLUMN_ORDER=("variant","r2","pvalue","molecular_trait_object_id",
    "molecular_trait_id","maf","gene_id","median_tpm","beta","se",
    "an","ac","ma_samples","chromosome","position","ref","alt","type","rsid")


def require(ok,msg):
    if not ok: raise ValueError(msg)


def verify_inputs(m,official,lift):
    gr=m.get("groups",[])
    require(m.get("schema")=="is-p3a-36-real-group-keys-from-legacy-393-v1"
            and m.get("scientific_status")=="HOLD"
            and len(gr)==36
            and len({(x["region_id"],x["QTL_dataset"],x["gene_id"]) for x in gr})==36
            and len({x["gene_id"] for x in gr})==20
            and Counter(x["QTL_dataset"] for x in gr)==
                {"QTD000131":10,"QTD000171":11,"QTD000434":15},
            "authentic archived candidate group manifest changed")
    require(official.get("schema")=="is-p3a-official-qtl-study-metadata-v1"
            and official.get("technical_status")==
                "PASS_3_SOURCE_IDENTITIES_IN_IMMUTABLE_CATALOGUE"
            and official.get("scientific_status")=="HOLD"
            and official.get("n_original_assay_groups")==36
            and official.get("qtl_reference_build")=="GRCh38"
            and official.get("metadata_sample_size_automatically_approved_per_gene") is False
            and set(official.get("QTL_datasets",{}))==set(EXPECTED_N),
            "official eQTL Catalogue source metadata not independently verified")
    for dataset,sample_n in EXPECTED_N.items():
        x=official["QTL_datasets"][dataset]
        require(x["nominal_catalogue_sample_size"]==sample_n
                and x["QTL_effect_allele"]=="ALT_GRCh38"
                and x["source_HEAD_probe"]["http_status"]==200
                and x["tabix_TBI_HEAD_probe"]["http_status"]==200
                and x["source_file_variant_rows_verified"] is False
                and x["official_individual_N_approved"] is False
                and x["source_https_probe"]==
                 "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"+
                 x["study_id"]+"/"+dataset+"/"+dataset+".all.tsv.gz",
                "source DOI/study endpoint or metadata wrongly adopted: "+dataset)
    refs=lift.get("source_provenance",{})
    a=lift.get("liftover",{})
    require(lift.get("schema")=="is-p3a-chr10-official-ucsc-grch37-to-grch38-anchor-ref-audit-v1"
            and lift.get("technical_status")=="PASS_SOURCE_MD5_AND_ANCHOR_REF_VALIDATION"
            and lift.get("scientific_status")=="HOLD"
            and lift.get("n_real_group_keys")==36
            and lift.get("formal_coloc_abf_eligible")==0
            and a.get("anchor_status")=="PASS_SIX_ANCHORS_AND_TWO_FASTAS"
            and a.get("all_unique_hg38_chr10_plus_monotonic_REF_concordant") is True
            and a.get("full_cis_region_equivalence_verified") is False
            and refs.get("hg38",{}).get("published_UCSC_MD5")==
                UCSC.FILES["hg38"]["md5"],
            "UCSC real FASTA or official chr10 mapped six anchors not approved")
    region=a.get("candidate_query_region_GRCh38_1based",{})
    require(region.get("chrom")=="10"
            and region.get("start")==102403919
            and region.get("end")==104460593
            and region.get("is_lower_bound_of_full_cis_with_guaranteed_coverage") is False,
            "wrong or falsely complete chr10 QTL query interval")
    return region


def new_stats():
    return {"region_rows_scanned":0,"candidate_gene_rows_seen":0,
            "candidate_unique_SNP_gene_pairs":0,
            "valid_unique_nonpal_SNP_gene_pairs":0,
            "duplicated_variant_rsid_rows":0,"conflicting_duplicates":0,
            "bad_REF_against_UCSC_GRCh38":0,"invalid_signed_beta_SE_MAF":0,
            "invalid_allele_or_variant_key":0,"non_SNP_or_palindromic":0,
            "p_zero_underflow":0,"noncandidate_rows":0}


def read_rows(lines,seq,genes,expected_region):
    """All original SNP effect values stay in process memory; emit aggregate QC only."""
    per_gene={g:{"candidate_rows":0,"valid_unique_SNPs":0,
                 "valid_nonpal_SNPs":0,"bad_ref":0,"bad_effect_SE_MAF":0,
                 "duplicate_rsid_records":0,"conflicting_duplicates":0} for g in genes}
    seen={g:{} for g in genes}
    stats=new_stats()
    total=0
    header=None
    digest=hashlib.sha256()
    nbytes=0
    for line in lines:
        nbytes+=len(line.encode("utf8"))
        require(nbytes<=MAX_BYTES,"indexed source returned unexpectedly large output")
        digest.update(line.encode("utf8"))
        row=line.rstrip("\r\n")
        if not row:continue
        cols=row.lstrip("#").split("\t")
        if header is None:
            if set(COLUMNS).issubset(set(cols)):
                header=cols
                require(len(set(cols))==len(cols),"duplicate official tabix column")
            elif row.startswith("#"):
                continue
            else:
                raise ValueError("missing official tabix -h column header")
            continue
        if row.startswith("#"):continue
        total+=1
        require(total<=MAX_LINES,"bounded region query result exceeded maximum")
        require(len(cols)==len(header),"QTL sumstats TSV truncated/incomplete")
        x=dict(zip(header,cols))
        stats["region_rows_scanned"]+=1
        chrom=x["chromosome"].removeprefix("chr")
        pos=int(x["position"])
        require(chrom=="10" and expected_region["start"]<=pos<=expected_region["end"],
                "remote tabix SNP outside mapped GRCh38 chr10 region")
        gene=x["gene_id"]
        if gene not in genes:
            stats["noncandidate_rows"]+=1
            continue
        stats["candidate_gene_rows_seen"]+=1
        pg=per_gene[gene]
        pg["candidate_rows"]+=1
        ref=x["ref"].upper()
        alt=x["alt"].upper()
        key=f"chr10_{pos}_{ref}_{alt}"
        if (x["variant"]!=key or len(ref)!=1 or len(alt)!=1
            or ref not in "ACGT" or alt not in "ACGT" or ref==alt):
            stats["invalid_allele_or_variant_key"]+=1
            continue
        if chr(seq[pos-1])!=ref:
            stats["bad_REF_against_UCSC_GRCh38"]+=1
            pg["bad_ref"]+=1
            continue
        try:
            beta=float(x["beta"]);se=float(x["se"]);maf=float(x["maf"])
            p=float(x["pvalue"])
        except (TypeError,ValueError):
            stats["invalid_signed_beta_SE_MAF"]+=1
            pg["bad_effect_SE_MAF"]+=1
            continue
        if not (all(math.isfinite(v) for v in (beta,se,maf,p))
                and se>0 and 0<maf<=.5 and 0<=p<=1):
            stats["invalid_signed_beta_SE_MAF"]+=1
            pg["bad_effect_SE_MAF"]+=1
            continue
        is_pal=frozenset((ref,alt)) in (frozenset(("A","T")),frozenset(("C","G")))
        if is_pal or x.get("type","SNP")!="SNP":
            stats["non_SNP_or_palindromic"]+=1
        if p==0:stats["p_zero_underflow"]+=1
        unique=(pos,ref,alt,x["molecular_trait_id"])
        signed=(x["beta"],x["se"],x["maf"],x["pvalue"])
        prev=seen[gene].get(unique)
        if prev is not None:
            stats["duplicated_variant_rsid_rows"]+=1
            pg["duplicate_rsid_records"]+=1
            if prev!=signed:
                stats["conflicting_duplicates"]+=1
                pg["conflicting_duplicates"]+=1
            continue
        seen[gene][unique]=signed
        stats["candidate_unique_SNP_gene_pairs"]+=1
        pg["valid_unique_SNPs"]+=1
        if not is_pal and x.get("type","SNP")=="SNP":
            stats["valid_unique_nonpal_SNP_gene_pairs"]+=1
            pg["valid_nonpal_SNPs"]+=1
    require(header is not None,"indexed QTL file did not return header")
    return {"metadata_original_header":header,
            "technical_status":"PASS_INDEXED_REGION_SOURCE_READ",
            "interval_rows":stats,
            "gene_rows":per_gene,"region_stream_SHA256":digest.hexdigest(),
            "processed_lines":total,
            "stream_bytes":nbytes}


def stream_tabix(url,region,seq,genes,open_tabix=None):
    """Use HTSlib bundled in pysam (network-enabled wheel), one indexed range.

    Ubuntu's apt tabix package may be built without HTTPS libcurl support.
    Official release-7 Column.md supplies the 19-field order if the tabix
    index has no commented header; column content and REF/variant are checked.
    """
    if open_tabix is None:
        import pysam
        open_tabix=pysam.TabixFile
    with open_tabix(url) as tabix:
        contigs=set(tabix.contigs)
        chrom="10" if "10" in contigs else "chr10" if "chr10" in contigs else None
        require(chrom is not None,"official Tabix index missing chromosome 10")
        raw_header=[v.decode("utf-8") if isinstance(v,bytes) else str(v)
                    for v in tabix.header]
        declared=any(set(COLUMNS).issubset(set(x.lstrip("#").split("\t")))
                     for x in raw_header)
        # No extra incomplete remote HTTP range request for header: official
        # 19 columns documented in the frozen eQTL Catalogue release source.
        incoming=(raw_header if declared else ["\t".join(OFFICIAL_COLUMN_ORDER)])
        def rows():
            for line in incoming:
                yield line+"\n"
            for line in tabix.fetch(chrom,region["start"]-1,region["end"]):
                yield (line.decode("utf-8") if isinstance(line,bytes) else line)+"\n"
        data=read_rows(rows(),seq,genes,region)
        require(data["interval_rows"]["region_rows_scanned"]>0,
                "original indexed eQTL returned no chr10 association rows")
        data["remote_contig_used"]=chrom
        data["header_origin"]=("ORIGINAL_TABIX_COMMENTED_HEADER"
                                if declared else
                                "PINNED_EQTL_CATALOGUE_COLUMN_DOCUMENTATION")
        return data


def main():
    parser=argparse.ArgumentParser()
    for name in ("manifest","official","liftover","out"):
        parser.add_argument("--"+name,type=Path,required=True)
    args=parser.parse_args()
    result={"schema":"is-p3a-36-original-eqtl-catalogue-indexed-chr10-audit-v1",
            "technical_status":"HOLD_INPUT_OR_REMOTE_SOURCE",
            "scientific_status":"HOLD","coloc_abf_eligible":0,
            "coloc_abf_executed":0,"coloc_susie_eligible":0,
            "coloc_susie_executed":0,"QTL_datasets":{},"groups":[]}
    try:
        manifest=json.loads(args.manifest.read_text())
        source=json.loads(args.official.read_text())
        lift=json.loads(args.liftover.read_text())
        region=verify_inputs(manifest,source,lift)
        by_dataset={}
        with tempfile.TemporaryDirectory(prefix="IS_EBI_QTL_temp") as temp:
            dest=Path(temp)/"hg38.chr10.fa.gz"
            ref=UCSC.fetch_checked("hg38",dest)
            require(ref["observed_SHA256"]==
                    lift["source_provenance"]["hg38"]["observed_SHA256"],
                    "reference FASTA differs from approved 6-anchor run")
            chrom_seq=UCSC.chromosome_fasta_sequence(dest)
            require(hashlib.sha256(chrom_seq).hexdigest()==
                    lift["chr10_uncompressed_fasta_base_SHA256"]["GRCh38_chr10"],
                    "reference bases not identical to official liftover run")
            for i,(dataset,n) in enumerate(EXPECTED_N.items()):
                if i:time.sleep(8)
                genes={g["gene_id"] for g in manifest["groups"]
                       if g["QTL_dataset"]==dataset}
                url=source["QTL_datasets"][dataset]["source_https_probe"]
                try:
                    data=stream_tabix(url,region,chrom_seq,genes)
                    data["source_url"]=url
                    data["dataset_sample_N_from_catalogue_not_gene_effective_N"]=n
                    by_dataset[dataset]=data
                except (ValueError,OSError,TimeoutError,UnicodeError,
                        subprocess.SubprocessError) as exc:
                    by_dataset[dataset]={
                        "technical_status":"HOLD_TABIX_OR_ORIGINAL_DATA",
                        "source_url":url,"error_type":type(exc).__name__,
                        "error":str(exc)[:280],"gene_rows":{}}
                    for later in list(EXPECTED_N)[i+1:]:
                        by_dataset[later]={
                            "technical_status":"HOLD_NOT_QUERIED_AFTER_FIRST_SOURCE_ERROR",
                            "source_url":source["QTL_datasets"][later]["source_https_probe"],
                            "gene_rows":{}}
                    break
        result["QTL_datasets"]=by_dataset
        result["converted_GRCh38_queried_context"]=region
        result["original_full_11GB_source_downloaded"]=False
        result["original_liftover_artifact_fullcis_not_verified"]=True
        result["gene_full_cis_window_coverage_verified"]=False
        result["all_variant_GRCh38_to_GRCh37_checked"]=False
        groups=[]
        for g in manifest["groups"]:
            qtl=g["QTL_dataset"]
            obj=by_dataset[qtl]
            record=obj.get("gene_rows",{}).get(g["gene_id"],{})
            groups.append({
               "region_id":g["region_id"],"QTL_dataset":qtl,
               "gene_id":g["gene_id"],"gene_symbol":g["gene_symbol"],
               "GRCh38_context_gene_associations_seen":record.get("candidate_rows"),
               "GRCh38_context_REF_beta_SE_MAF_valid_unique":
                   record.get("valid_unique_SNPs"),
               "GRCh38_context_valid_nonpal_unique":
                   record.get("valid_nonpal_SNPs"),
               "gene_full_cis_per_variant_proven":False,
               "GWAS_GRCh37_matched_source_N_and_allele":False,
               "coloc_abf_eligible":False,"coloc_susie_eligible":False,
               "status":"HOLD_INCOMPLETE_GENE_CIS_AND_GWAS_SOURCE_N",
            })
        result["groups"]=groups
        require(len(groups)==36,"all real assay groups must appear")
        n_pass=sum(x["technical_status"]=="PASS_INDEXED_REGION_SOURCE_READ"
                    for x in by_dataset.values())
        result["source_datasets_with_actual_indexed_rows_read"]=n_pass
        result["technical_status"]=("PASS_3_OFFICIAL_INDEXED_QTL_CONTEXT_REPORTS"
                                   if n_pass==3 else "HOLD_ONE_OR_MORE_TABIX_SOURCES")
    except (ValueError,KeyError,TypeError,OSError,ImportError) as exc:
        result["technical_status"]="HOLD_INPUT_OR_REFERENCE_INTEGRITY"
        result["error"]=str(exc)[:280]
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({
       "technical_status":result["technical_status"],
       "scientific_status":"HOLD",
       "source_datasets_with_actual_indexed_rows_read":
            result.get("source_datasets_with_actual_indexed_rows_read",0),
       "total_groups":len(result["groups"]),
       "source_status":{k:v["technical_status"]
                        for k,v in result["QTL_datasets"].items()},
       "error":result.get("error"),
    },ensure_ascii=False))
    return 0 if result["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
