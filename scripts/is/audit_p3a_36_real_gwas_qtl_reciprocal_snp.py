#!/usr/bin/env python3
"""P3A: bidirectional FASTA-backed single-SNP QTL--GWAS source harmonization.

Scan all 36 genuinely recovered QTL gene-tissue groups, not just top eGenes.
eQTL Catalogue nominal source: GRCh38 ALT effect; original GIGASTROKE
published EAS/EUR GWAS: GRCh37 effect_allele (not necessarily ALT).
Input sources and reference chains are pinned to original publisher digests.
Only aggregate QC counts are saved; no SNP effect rows leave the runner.

No approved GWAS case/control N, allele-effective QTL N, study-specific
complete cis test enumeration, INDEL harmonization, or coloc calculations.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib
import json
import math
import tempfile
import time
import urllib.request
from collections import Counter,defaultdict
from pathlib import Path

QTL=importlib.import_module("scripts.is.audit_p3a_indexed_qtl_real_chr10")
CIS=importlib.import_module("scripts.is.audit_p3a_official_cis_union_36")
UCSC=importlib.import_module("scripts.is.audit_p3a_chr10_ucsc_liftover_fasta")
GWAS=importlib.import_module("scripts.is.audit_p3a_chr10_original_publisher_gwas")
PINNED=importlib.import_module("scripts.is.p3a_publisher_md5_gate")

REVERSE_CHAIN_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
REVERSE_CHAIN_MD5="ff3031d93792f4cbb86af44055efd903"
MAX_REVERSE_CHAIN_BYTES=3_000_000
EXPECTED_DATASETS=CIS.EXPECTED_GROUPS
REFERENCE_UNION=CIS.EXPECTED_WINDOW
GWAS_ANCESTRIES=("EAS","EUR")
SOURCE_HEADER=QTL.OFFICIAL_COLUMN_ORDER
COMP=str.maketrans("ACGT","TGCA")
PASS_STATUS="PASS_QTL_GWAS_SOURCE_SNP_HARMONIZATION_AGGREGATES_ONLY"


def require(x,reason):
    if not x:raise ValueError(reason)


def fetch_reverse_chain(dest,fetch=urllib.request.urlopen):
    req=urllib.request.Request(REVERSE_CHAIN_URL,
         headers={"User-Agent":"IS-P3A-Bidirectional-liftover-genomic-source-QC/1.0"})
    md5=hashlib.md5(usedforsecurity=False)
    sha=hashlib.sha256()
    n=0
    with fetch(req,timeout=90) as response:
        require(getattr(response,"status",200)==200,
                "official UCSC reverse chain HTTP status not 200")
        with dest.open("wb") as o:
            while True:
                b=response.read(1024*1024)
                if not b:break
                n+=len(b)
                require(n<=MAX_REVERSE_CHAIN_BYTES,"reverse chain size bound exceeded")
                md5.update(b);sha.update(b);o.write(b)
    require(n>300_000 and md5.hexdigest()==REVERSE_CHAIN_MD5,
            "reverse GRCh38→GRCh37 chain official MD5 mismatch: "+
            md5.hexdigest())
    return {"url":REVERSE_CHAIN_URL,
            "official_publisher_md5":REVERSE_CHAIN_MD5,
            "observed_md5":md5.hexdigest(),
            "observed_sha256":sha.hexdigest(),
            "compressed_bytes":n,"verified":True}


def unique_position(converter,pos):
    mappings=converter.convert_coordinate("chr10",pos-1) or []
    if len(mappings)!=1:
        return None,"HOLD_NO_UNIQUE_CHAIN_MAPPING"
    chrom,coordinate,strand,_=mappings[0]
    if chrom!="chr10" or strand!="+" or not isinstance(coordinate,(int,float)) \
            or not math.isfinite(coordinate) or coordinate!=int(coordinate):
        return None,"HOLD_NONCHR10_OR_NEGATIVE_STRAND"
    return int(coordinate)+1,None


def reciprocal_allele_check(chr38_pos,ref,old_seq,new_seq,to37,to38):
    require(ref in "ACGT" and len(ref)==1 and chr38_pos>0,
            "invalid QTL REF for a single-SNP liftover")
    if chr38_pos>len(new_seq) or chr(new_seq[chr38_pos-1])!=ref:
        return None,"HOLD_GRCh38_REFERENCE_MISMATCH"
    pos37,reason=unique_position(to37,chr38_pos)
    if reason:return None,reason
    if pos37<=0 or pos37>len(old_seq) or chr(old_seq[pos37-1])!=ref:
        return None,"HOLD_GRCh37_REFERENCE_MISMATCH"
    verify38,reason=unique_position(to38,pos37)
    if reason or verify38!=chr38_pos:
        return None,"HOLD_RECIPROCAL_LIFTOVER_DISCORDANCE"
    return pos37,None


def harmonize(qtl_ref,qtl_alt,gw_ea,gw_oa,gw_beta,gw_eaf):
    """Return signed GWAS effect on molecular QTL *ALT* dosage.

    The SNP orientation is source-anchored to GRCh37 and GRCh38 FASTAs;
    complemented-only pairs and palindrome are never independently approved.
    """
    ref=qtl_ref.upper();alt=qtl_alt.upper()
    ea=gw_ea.upper();oa=gw_oa.upper()
    if len({ref,alt})!=2 or any(z not in "ACGT" for z in (ref,alt,ea,oa)):
        return "HOLD_NONCANONICAL_ALLELES",None,None
    if {ref,alt} in ({"A","T"},{"C","G"}):
        return "HOLD_PALINDROMIC_STRAND",None,None
    if {ea,oa}!={ref,alt}:
        if {ea.translate(COMP),oa.translate(COMP)}=={ref,alt}:
            return "HOLD_COMPLEMENT_ONLY_NO_ORIENTATION_APPROVAL",None,None
        return "HOLD_ALLELE_PAIR_MISMATCH",None,None
    if not (math.isfinite(gw_beta) and math.isfinite(gw_eaf) and 0<=gw_eaf<=1):
        return "HOLD_INVALID_GWAS_BETA_OR_EAF",None,None
    if ea==alt:
        return "MATCH_ALT_EFFECT",gw_beta,gw_eaf
    return "MATCH_REF_EFFECT_SIGN_FLIPPED",-gw_beta,1-gw_eaf


def qtl_record(fields):
    try:
        chromosome=fields["chromosome"].removeprefix("chr")
        pos=int(fields["position"])
    except (ValueError,TypeError):
        return None
    if chromosome!="10" or not (REFERENCE_UNION["start"]<=pos<=REFERENCE_UNION["end"]):
        raise ValueError("nominal QTL original variant outside expected GRCh38 source")
    ref=fields["ref"].upper()
    alt=fields["alt"].upper()
    if fields["variant"]!=f"chr10_{pos}_{ref}_{alt}" or len(ref)!=1 or len(alt)!=1 \
            or ref not in "ACGT" or alt not in "ACGT" or ref==alt or fields["type"]!="SNP":
        return None
    try:
        beta=float(fields["beta"]);se=float(fields["se"])
        maf=float(fields["maf"]);p=float(fields["pvalue"])
    except (TypeError,ValueError):
        return None
    if not (all(math.isfinite(x) for x in (beta,se,maf,p))
             and se>0 and 0<maf<=.5 and 0<=p<=1):
        return None
    if {ref,alt} in ({"A","T"},{"C","G"}):
        return None
    try:
        an=float(fields.get("an",""));ac=float(fields.get("ac",""))
        af=ac/an if math.isfinite(an) and math.isfinite(ac) and an>0 and 0<=ac<=an else None
    except (ValueError,ZeroDivisionError):
        af=None
    # ac/an is within-QTL ALT frequency, not reference MAF or approved N.
    return {"pos38":pos,"ref":ref,"alt":alt,"beta":beta,"se":se,
            "maf":maf,"alt_af":af,"molecular_trait_id":fields["molecular_trait_id"]}


def gather_qtl(url,groupgenes,seq37,seq38,to37,to38,open_tabix=None):
    if open_tabix is None:
        import pysam
        open_tabix=pysam.TabixFile
    by_gene={gene:{"n_original_qtl_valid_nonpal":0,
                   "unique_qtl_lifted":0,
                   "all_true_SNP_QTL_allele_pair_ready":False,
                   "n_liftover_blockers":{},
                   "n_QTL_ALT_AF_from_observed_ac_an":0}
             for gene in groupgenes}
    source_counts={"candidate_SNP_rows_visited":0,
                   "candidate_unique_nonpal_SNP_rows":0,
                   "unique_liftover_QTL_variants":0,
                   "duplicate_original_rsid_records":0,
                   "QTL_reference_liftover_reasons":Counter(),
                   "QTL_source_sha256_stream":hashlib.sha256()}
    snps_by_gene={}
    cache={}
    seen={gene:set() for gene in groupgenes}
    with open_tabix(url) as reader:
        names=set(reader.contigs)
        chrom="10" if "10" in names else "chr10" if "chr10" in names else None
        require(chrom is not None,"QTL Tabix index missing chromosome 10")
        for line in reader.fetch(chrom,REFERENCE_UNION["start"]-1,REFERENCE_UNION["end"]):
            if isinstance(line,bytes):line=line.decode("utf8")
            source_counts["QTL_source_sha256_stream"].update((line.rstrip("\r\n")+"\n").encode("utf8"))
            cells=line.rstrip("\r\n").split("\t")
            require(len(cells)==len(SOURCE_HEADER),"original QTL source physical TSV columns changed")
            source=dict(zip(SOURCE_HEADER,cells))
            gene=source["gene_id"]
            if gene not in groupgenes:continue
            item=qtl_record(source)
            if item is None:continue
            source_counts["candidate_SNP_rows_visited"]+=1
            pos=item["pos38"]
            key=(pos,item["ref"],item["alt"],item["molecular_trait_id"])
            if key in seen[gene]:
                source_counts["duplicate_original_rsid_records"]+=1
                continue
            seen[gene].add(key)
            source_counts["candidate_unique_nonpal_SNP_rows"]+=1
            by_gene[gene]["n_original_qtl_valid_nonpal"]+=1
            posref=(pos,item["ref"])
            if posref not in cache:
                cache[posref]=reciprocal_allele_check(pos,item["ref"],
                    seq37,seq38,to37,to38)
            pos37,reason=cache[posref]
            if reason:
                source_counts["QTL_reference_liftover_reasons"][reason]+=1
                tally=by_gene[gene]["n_liftover_blockers"]
                tally[reason]=tally.get(reason,0)+1
                continue
            by_gene[gene]["unique_qtl_lifted"]+=1
            if item["alt_af"] is not None:
                by_gene[gene]["n_QTL_ALT_AF_from_observed_ac_an"]+=1
            snps_by_gene.setdefault(gene,[]).append((pos37,item["ref"],
                     item["alt"],item["beta"],item["se"],item["alt_af"]))
    source_counts["unique_liftover_QTL_variants"]=len({
        (k[0],k[1]) for k,v in cache.items() if v[1] is None})
    source_counts["QTL_source_sha256_stream"]=source_counts[
        "QTL_source_sha256_stream"].hexdigest()
    source_counts["QTL_reference_liftover_reasons"]=dict(source_counts[
        "QTL_reference_liftover_reasons"])
    return snps_by_gene,by_gene,source_counts


def parse_original_gwas_file(path,ancestry,wanted_positions):
    """Read publisher original full file. Keep at most target region SNPs."""
    original_header=None
    total=0
    hits={}
    invalid=0
    duplicates=0
    with gzip.open(path,"rt",encoding="utf-8-sig",newline="") as handle:
        reader=csv.DictReader(handle,delimiter="\t")
        original_header=list(reader.fieldnames or [])
        require(set(GWAS.HEADER).issubset(set(original_header)),
                "original GWAS publisher signed beta/se/eaf or allele columns absent")
        for row in reader:
            total+=1
            require(total<=GWAS.MAX_GWAS_ROWS,"original GWAS file unexpectedly large")
            if row["chromosome"].removeprefix("chr")!="10":
                continue
            try:pos=int(row["base_pair_location"])
            except (TypeError,ValueError):continue
            if pos not in wanted_positions:continue
            # Parse the exact extended union, not the prior 2.05Mb-only original GWAS screen.
            try:
                ea=row["effect_allele"].upper()
                oa=row["other_allele"].upper()
                beta=float(row["beta"])
                se=float(row["standard_error"])
                p=float(row["p_value"])
                eaf=float(row["effect_allele_frequency"])
                good=ea in "ACGT" and oa in "ACGT" and ea!=oa and \
                     all(math.isfinite(x) for x in (beta,se,p,eaf)) and \
                     se>0 and 0<=p<=1 and 0<=eaf<=1
            except (KeyError,ValueError,TypeError):
                good=False
            if not good:
                invalid+=1
                continue
            if pos in hits:
                duplicates+=1
                hits[pos]=None
            else:
                hits[pos]=(ea,oa,beta,se,eaf)
    require(total==PINNED.SOURCES[ancestry]["original_rows_in_prior_source_qc"],
            "whole original GWAS parsed count differs from immutable publisher scan")
    return hits,{"GWAS_whole_original_row_count":total,
        "GWAS_original_target_unique_position_count":len(hits),
        "GWAS_target_invalid_numeric_or_allele_rows":invalid,
        "GWAS_target_duplicate_locus_rows":duplicates,
        "original_publisher_accession":PINNED.SOURCES[ancestry]["accession"],
        "original_GWAS_build":"GRCh37",
        "original_GWAS_effect_allele_role":"SOURCE_DEFINED_EFFECT_NOT_ALT",
        "GWAS_N_cases_controls_approved":False}


def summarize_groups(manifest,union,per_dataset,raw_snps,gwas_maps,gwas_counts,chain):
    groups=[]
    baseline={(g["QTL_dataset"],g["gene_id"]):g for g in union["groups"]}
    overall=Counter()
    for record in manifest["groups"]:
        group=record["QTL_dataset"]
        gene=record["gene_id"]
        key=(group,gene)
        entry=per_dataset[group][gene]
        expected=baseline[key]["extended_context_valid_nonpal_SNP_assay_rows"]
        require(entry["n_original_qtl_valid_nonpal"]==expected,
                "source original QTL record count drift: "+str(key))
        snps=raw_snps[group].get(gene,[])
        per_ancestry={}
        for ancestry in GWAS_ANCESTRIES:
            gwas=gwas_maps[ancestry]
            counts=Counter()
            AFabs_sum=0.0
            n_AF=0
            for pos,ref,alt,qbeta,qse,alt_af in snps:
                if pos not in gwas:
                    counts["NO_GWAS_POSITION"]+=1
                    continue
                if gwas[pos] is None:
                    counts["HOLD_DUPLICATE_GWAS_POSITION"]+=1
                    continue
                ea,oa,beta,se,eaf=gwas[pos]
                status,aligned,ALT_eaf=harmonize(ref,alt,ea,oa,beta,eaf)
                counts[status]+=1
                if status.startswith("MATCH_") and alt_af is not None:
                    AFabs_sum+=abs(alt_af-ALT_eaf)
                    n_AF+=1
                    if abs(alt_af-ALT_eaf)>.2:
                        counts["AF_ABS_DELTA_GT_0_2_DESCRIPTIVE"]+=1
            nmatched=counts["MATCH_ALT_EFFECT"]+counts["MATCH_REF_EFFECT_SIGN_FLIPPED"]
            per_ancestry[ancestry]={
                "dataset":PINNED.SOURCES[ancestry]["accession"],
                "counts_by_signed_GWAS_allele_match":dict(counts),
                "n_nonpal_QTL_GWAS_recipient_allele_matches":nmatched,
                "n_valid_QTL_ALT_AF_vs_GWAS_ALT_EAF_pairs":n_AF,
                "mean_absolute_frequency_difference_descriptive":
                    round(AFabs_sum/n_AF,6) if n_AF else None,
                "cannot_use_as_independent_replication":True,
                "sample_N_and_case_fraction_approved":False,
                "formal_coloc_eligible":False,
            }
            overall[(ancestry,"match")]+=nmatched
            overall[(ancestry,"no_GWAS_pos")]+=counts["NO_GWAS_POSITION"]
        groups.append({
          "region_id":record["region_id"],"QTL_dataset":group,
          "gene_id":gene,"gene_symbol":record["gene_symbol"],
          "QTL_source_indexed_nonpal_valid":entry["n_original_qtl_valid_nonpal"],
          "QTL_uniquely_reciprocally_lifted_REF_valid":entry["unique_qtl_lifted"],
          "QTL_liftover_blockers":entry["n_liftover_blockers"],
          "per_ancestry":per_ancestry,
          "exact_original_study_tested_QTL_cis_approved":False,
          "original_gene_specific_effective_N_approved":False,
          "GIGASTROKE_exact_trait_cases_controls_approved":False,
          "palindromic_and_indels_aligned":False,
          "formal_coloc_abf_eligible":False,
          "formal_coloc_susie_eligible":False,
          "science_status":"HOLD_N_AND_FULL_CIS_AND_STRAND_INDELS",
        })
    require(len(groups)==36 and len({x["gene_id"] for x in groups})==20,
            "36 gene-tissue identity records lost during original SNP matching")
    return {
        "schema":"is-p3a-36-real-source-reciprocal-SNP-GWAS-harmonization-v1",
        "technical_status":PASS_STATUS,"scientific_status":"HOLD",
        "n_real_assay_groups":36,"n_unique_genes":20,
        "n_independent_stroke_loci":1,
        "two_way_chain":chain,
        "original_20_gene_reference_union_GRCh38":REFERENCE_UNION,
        "QTL_datasets":{
           d:{"source_id":d,"source_study_n_NOT_gene_effective_N":QTL.EXPECTED_N[d],
              "source_stats":raw_snps[d]["__source_stats__"],
              "gene_n":len(per_dataset[d])}
             for d in EXPECTED_DATASETS},
        "GWAS_ancestries":gwas_counts,
        "matched_nonpal_summary":dict((f"{a}_{k}",v)
                                      for (a,k),v in overall.items()),
        "n_formal_coloc_abf_eligible":0,"n_formal_coloc_abf_executed":0,
        "n_formal_coloc_susie_eligible":0,"n_formal_coloc_susie_executed":0,
        "original_study_QTL_cis_variant_set_approved":False,
        "original_GWAS_phenotype_N_case_control_approved":False,
        "BBJ_independent_replication_proven":False,
        "full_variant_indel_palindromic_harmonization_proven":False,
        "rows":groups,
    }


def main():
    p=argparse.ArgumentParser()
    for field in ("manifest","union","official","liftover","out"):
        p.add_argument("--"+field,type=Path,required=True)
    a=p.parse_args()
    result={"schema":"is-p3a-36-real-source-reciprocal-SNP-GWAS-harmonization-v1",
            "technical_status":"HOLD_NO_FINAL_SOURCE_VALIDATION",
            "scientific_status":"HOLD",
            "n_formal_coloc_abf_eligible":0,"n_formal_coloc_abf_executed":0,
            "n_formal_coloc_susie_eligible":0,"n_formal_coloc_susie_executed":0}
    try:
        manifest=json.loads(a.manifest.read_text(encoding="utf8"))
        union=json.loads(a.union.read_text(encoding="utf8"))
        official=json.loads(a.official.read_text(encoding="utf8"))
        lift=json.loads(a.liftover.read_text(encoding="utf8"))
        require(manifest["schema"]=="is-p3a-36-real-group-keys-from-legacy-393-v1"
                and manifest["scientific_status"]=="HOLD"
                and len(manifest["groups"])==36
                and len({(g["region_id"],g["QTL_dataset"],g["gene_id"])
                         for g in manifest["groups"]})==36,
                "original group manifest missing or duplicated")
        require(union["schema"]=="is-p3a-36-expanded-original-cis-reference-union-v1"
                and union["technical_status"]=="PASS_THREE_REAL_INDEXED_REFERENCE_UNION_SURVEYS"
                and union["scientific_status"]=="HOLD"
                and union["cis_union_reference_Ensembl105_GRCh38_1based"]==REFERENCE_UNION
                and union["n_genuine_assay_groups"]==36
                and union["n_approved_coloc_abf"]==0,
                "source 36 gene cis reference union has not been independently verified")
        require(lift["technical_status"]=="PASS_SOURCE_MD5_AND_ANCHOR_REF_VALIDATION"
                and lift["liftover"]["anchor_status"]=="PASS_SIX_ANCHORS_AND_TWO_FASTAS",
                "original official source 6-point reference liftover unverified")
        require(official["technical_status"]=="PASS_3_SOURCE_IDENTITIES_IN_IMMUTABLE_CATALOGUE",
                "official eQTL source metadata unverified")
        with tempfile.TemporaryDirectory(prefix="IS_P3A_public_liftover_GWAS_") as td:
            folder=Path(td)
            source_files={}
            for n in ("chain","hg19","hg38"):
                source_files[n]=UCSC.fetch_checked(n,folder/(n+".gz"))
            require(all(source_files[k]["observed_SHA256"]==
                        lift["source_provenance"][k]["observed_SHA256"]
                        for k in source_files),"original official reference source SHA256 drift")
            gz=folder/"reverse.chain.gz"
            reverse=fetch_reverse_chain(gz)
            with gzip.open(gz,"rb") as src,(folder/"reverse.chain").open("wb") as dst:
                while (part:=src.read(1024*1024)):dst.write(part)
            with gzip.open(folder/"chain.gz","rb") as src,(folder/"forward.chain").open("wb") as dst:
                while (part:=src.read(1024*1024)):dst.write(part)
            from pyliftover import LiftOver
            to37=LiftOver(str(folder/"reverse.chain"))
            to38=LiftOver(str(folder/"forward.chain"))
            old_seq=UCSC.chromosome_fasta_sequence(folder/"hg19.gz")
            new_seq=UCSC.chromosome_fasta_sequence(folder/"hg38.gz")
            require(hashlib.sha256(old_seq).hexdigest()==
                    lift["chr10_uncompressed_fasta_base_SHA256"]["GRCh37_chr10"]
                    and hashlib.sha256(new_seq).hexdigest()==
                    lift["chr10_uncompressed_fasta_base_SHA256"]["GRCh38_chr10"],
                    "reference raw bases drift")
            raw={}
            per={}
            for index,(dataset,n) in enumerate(EXPECTED_DATASETS.items()):
                if index:time.sleep(9)
                genes={g["gene_id"] for g in manifest["groups"]
                       if g["QTL_dataset"]==dataset}
                url=official["QTL_datasets"][dataset]["source_https_probe"]
                values,summary,counts=gather_qtl(url,genes,old_seq,new_seq,to37,to38)
                require(counts["candidate_unique_nonpal_SNP_rows"]>=1,
                        "no genuine source QTL valid SNPs found")
                raw[dataset]=values
                raw[dataset]["__source_stats__"]=counts
                per[dataset]=summary
            positions={item[0] for dataset in EXPECTED_DATASETS
                       for gene,values in raw[dataset].items()
                       if gene!="__source_stats__"
                       for item in values}
            require(500<len(positions)<500_000,"wrong QTL reciprocal liftOver position range")
            gwas_maps={}
            gwas_counts={}
            for ancestry in GWAS_ANCESTRIES:
                dest=folder/(ancestry+"_original_GWAS.tsv.gz")
                source=GWAS.stream_source(ancestry,dest)
                gmap,counts=parse_original_gwas_file(dest,ancestry,positions)
                counts["original_publisher_full_file_MD5_SHA256_bytes_pass"]=source[
                       "publisher_md5_and_prior_sha256_verified"]
                gwas_maps[ancestry]=gmap
                gwas_counts[ancestry]=counts
            result=summarize_groups(manifest,union,per,raw,gwas_maps,gwas_counts,reverse)
    except (OSError,ValueError,KeyError,TypeError,OverflowError,ImportError,
            UnicodeError,LookupError,TimeoutError) as exc:
        result["technical_status"]="HOLD_SOURCE_OR_ALLELE_LIFTOVER_FAILURE"
        result["error_type"]=type(exc).__name__
        result["error"]=str(exc)[:300]
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:result.get(k) for k in
       ("technical_status","scientific_status","n_real_assay_groups",
        "matched_nonpal_summary","n_formal_coloc_abf_eligible","error")},
        ensure_ascii=False))
    return 0 if result["technical_status"]==PASS_STATUS else 1


if __name__=="__main__":
    raise SystemExit(main())
