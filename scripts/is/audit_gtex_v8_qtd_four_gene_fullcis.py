#!/usr/bin/env python3
"""Read-only bounded full-cis original GTEx V8 versus EBI QTD000161 audit.

Two predeclared GRCh38 windows cover the *four-gene reprocessed QTD tested
universe* from the 2026-10-11 source-lineage audit. Every source is a separate
analytical model. This script does not calculate MR, coloc or causal inference.

A 2.5Mb chr4 and a 5Mb chr10 window are each queried once per source using
tabix; prohibit excessive row-count/memory, retain only four target genes.
No paywalled source or whole-file 3.4GB download.
"""
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import statistics
import time
import urllib.request
import zlib

HERE=Path(__file__).resolve().parent
SPEC=importlib.util.spec_from_file_location("gtex_qtd_bounded",HERE/"audit_gtex_v8_vs_qtd_direct_tabix.py")
BASE=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)

# Exact ranges used in prior independent EBI QTD full-cis reconstruction.
FULL_WINDOWS=[("4",79000000,81500000),("10",101500000,106500000)]
EXPECTED_REPROCESSED_CIS={"FGF5":5511,"C4orf22":5490,"CALHM2":4868,"NEURL":4820}
MAX_ROWS_PER_QUERY=1_000_000
MAX_TARGET_GENE_VARIANTS_PER_QUERY=30000
MAX_INVALID_TARGET_ROWS=100
WAIT_SECONDS=2.5

def variant_type(key):
    _gene,_chr,_pos,ref,alt=key
    return "SNP" if len(ref)==len(alt)==1 else "INDEL_MNV"


def diagnostic_unmatched_qtd(og,reprocessed, *, limit=10):
    # The exact source key is (stable_gene, chr, GRCh38 position, REF, ALT).
    # Nearby or trimmed-equivalent alleles are only candidates, never used
    # as proof of identity without source-specific normalization/ref FASTA.
    original_by_gene={}
    for key in og:
        original_by_gene.setdefault((key[0],key[1]),[]).append(key)
    unmatched=sorted(set(reprocessed)-set(og))
    reports=[]
    for key in unmatched[:limit]:
        gene,chrom,pos,ref,alt=key
        near=[candidate for candidate in original_by_gene.get((gene,chrom),[])
              if abs(candidate[2]-pos)<=5]
        near.sort(key=lambda candidate:(abs(candidate[2]-pos),candidate))
        exact_position=[candidate for candidate in near if candidate[2]==pos]
        reverse_at_position=any(candidate[3]==alt and candidate[4]==ref
                                for candidate in exact_position)
        candidate_rows=[{
            "original_GRCh38_pos":c[2],"original_ref":c[3],"original_alt":c[4],
            "bp_distance_from_qtd":c[2]-pos,
            "reversed_REF_ALT_at_same_position":(c[2]==pos and c[3]==alt and c[4]==ref),
        } for c in near[:10]]
        reports.append({
            "gene_id":gene,"chr":chrom,"qtd_GRCh38_pos":pos,
            "qtd_REF":ref,"qtd_ALT":alt,
            "qtd_variant_type":variant_type(key),
            "qtd_beta":reprocessed[key][0],"qtd_se":reprocessed[key][1],
            "qtd_p":reprocessed[key][2],
            "original_same_gene_same_position_row_count":len(exact_position),
            "original_nearby_plusminus5bp_row_count":len(near),
            "ref_alt_reversed_same_position":reverse_at_position,
            "original_candidate_variants_nearby_bounded":candidate_rows,
            "independent_ref_fasta_indel_normalization_checked":False,
            "matched_after_liftover_or_strand_reversal":False,
            "diagnostic_only_not_harmonized":True,
        })
    return {
        "qtd_only_variant_count":len(unmatched),
        "qtd_only_diagnostics_truncated":len(unmatched)>limit,
        "qtd_only_diagnostics":reports,
    }


def read_source_rows(tabix,header,region, *, max_rows=MAX_ROWS_PER_QUERY):
    chrom,start,end=region
    contig=BASE.known_index_contig(tabix.contigs,chrom)
    if contig is None:
        return {"status":"HOLD_MISSING_CONTIG","rows":{},"scanned":0,"target_rows":0}
    if not BASE.REQUIRED.issubset(header or []):
        return {"status":"HOLD_REQUIRED_COLUMNS","rows":{},"scanned":0,"target_rows":0}
    ix={k:i for i,k in enumerate(header)}
    allowed=set(BASE.TARGET_GENES.values())
    parsed={}
    counts=Counter()
    invalid_source_rows_bounded=[]
    try:
        for line in tabix.fetch(contig,start-1,end):
            counts["region_total"]+=1
            if counts["region_total"]>max_rows:
                return {"status":"HOLD_SCAN_CAP","rows":{},"scanned":counts["region_total"],
                        "target_rows":counts["target_seen"]}
            fields=line.rstrip("\n").split("\t")
            if len(fields)!=len(header):
                return {"status":"HOLD_ROW_WIDTH_MISMATCH","rows":{},"scanned":counts["region_total"],
                        "target_rows":counts["target_seen"]}
            gene=fields[ix["molecular_trait_id"]].split(".")[0]
            if gene not in allowed:
                continue
            counts["target_seen"]+=1
            if counts["target_seen"]>MAX_TARGET_GENE_VARIANTS_PER_QUERY:
                return {"status":"HOLD_TARGET_ROW_CAP","rows":{},"scanned":counts["region_total"],
                        "target_rows":counts["target_seen"]}
            try:
                c=fields[ix["chromosome"]].removeprefix("chr")
                pos=int(fields[ix["position"]])
                ref,alt=(fields[ix[x]].upper() for x in ("ref","alt"))
                p,beta,se=(float(fields[ix[x]]) for x in ("pvalue","beta","se"))
                # Include explicit REF/ALT SNPs, indels and multinucleotide
                # changes up to 200bp, but never symbolic/ambiguous alleles.
                # Keep exact alleles, do not assume strand/normalization equality.
                if (not ref or not alt or ref==alt or len(ref)>200 or len(alt)>200
                    or any(nt not in "ACGT" for nt in ref+alt)):
                    counts["unsupported_allele_encoding_rows"]+=1
                    continue
                counts["single_base_SNP_rows" if len(ref)==len(alt)==1
                       else "complex_indel_mnv_rows"]+=1
                if not (math.isfinite(beta) and math.isfinite(se) and se>0 and
                        math.isfinite(p) and 0<=p<=1 and pos>=1):
                    counts["invalid_statistic_rows"]+=1
                    if len(invalid_source_rows_bounded)<40:
                        invalid_source_rows_bounded.append({
                            "gene_id":gene,"chromosome":c,"position":pos,
                            "ref":ref,"alt":alt,
                            "p_raw":fields[ix["pvalue"]],
                            "beta_raw":fields[ix["beta"]],
                            "se_raw":fields[ix["se"]],
                            "reason":"NONFINITE_OR_NONPOSITIVE_SE_OR_INVALID_P",
                        })
                    continue
                key=(gene,c,pos,ref,alt)
                value=(beta,se,p)
                if key in parsed:
                    counts["duplicate_exact_id_rows"]+=1
                    if parsed[key]!=value:
                        return {"status":"HOLD_CONFLICTING_DUPLICATE_SOURCE_VARIANT",
                                "rows":{},"scanned":counts["region_total"],
                                "target_rows":counts["target_seen"]}
                else:
                    parsed[key]=value
            except (KeyError,TypeError,ValueError,OverflowError):
                counts["invalid_statistic_rows"]+=1
                if len(invalid_source_rows_bounded)<40:
                    invalid_source_rows_bounded.append({
                        "gene_id":gene,
                        "position_raw":fields[ix["position"]],
                        "ref_raw":fields[ix["ref"]],
                        "alt_raw":fields[ix["alt"]],
                        "p_raw":fields[ix["pvalue"]],
                        "beta_raw":fields[ix["beta"]],
                        "se_raw":fields[ix["se"]],
                        "reason":"UNPARSEABLE_ROW",
                    })
        outcome={
            "status":"SCREEN_VALID" if counts["invalid_statistic_rows"]<=MAX_INVALID_TARGET_ROWS
              else "HOLD_TOO_MANY_INVALID_TARGET_ROWS",
            "rows":parsed,"scanned":counts["region_total"],
            "target_rows":counts["target_seen"],"counts":dict(counts),
            "invalid_source_rows_bounded":invalid_source_rows_bounded}
        return outcome
    except (OSError,ValueError,RuntimeError) as exc:
        return {"status":"HOLD_TABIX_REMOTE_FETCH", "rows":{},
                "scanned":counts["region_total"],"target_rows":counts["target_seen"],
                "failure_class":type(exc).__name__}


def summarize_full_cis(original,reprocessed):
    out={}
    for name,gene in BASE.TARGET_GENES.items():
        a={k:v for k,v in original["rows"].items() if k[0]==gene}
        b={k:v for k,v in reprocessed["rows"].items() if k[0]==gene}
        obj=BASE.compare_rows(
            {"status":original["status"],"rows":a,"observed":len(a)},
            {"status":reprocessed["status"],"rows":b,"observed":len(b)})
        orig=set(a)
        reprs=set(b)
        obj.update({
            "gene":name,"gene_id":gene,
            "original_unique_tested_variant":len(a),
            "reprocessed_unique_tested_variant":len(b),
            "original_only_variant":len(orig-reprs),
            "reprocessed_only_variant":len(reprs-orig),
            "original_only_SNP_count":sum(variant_type(k)=="SNP" for k in orig-reprs),
            "original_only_indel_mnv_count":sum(variant_type(k)!="SNP" for k in orig-reprs),
            "reprocessed_only_diagnostic":diagnostic_unmatched_qtd(a,b),
            "is_complete_full_cis_only_if_separate_cis_interval_confirmed":False,
            "reprocessed_expected_previous_QTD_n_variants":EXPECTED_REPROCESSED_CIS[name],
            "reprocessed_variant_denominator_matches_prior_QTD_fullcis_audit":
                len(b)==EXPECTED_REPROCESSED_CIS[name],
        })
        out[name]=obj
    return out


def run(tabix=None,header_fetcher=None,md5_fetcher=None,wait=WAIT_SECONDS):
    result={"schema":"is-GTExv8-original-vs-QTD000161-four-gene-full-cis-v3-unmatched-exact-variant-audit",
       "original":BASE.ORIGINAL,"reprocessed":BASE.REPROCESSED,
       "GRCh38_predeclared_windows":[f"{c}:{s}-{e}" for c,s,e in FULL_WINDOWS],
       "QTL_source_tissue":"Brain_Cerebellar_Hemisphere",
       "original_full_file_checksum_recomputed":False,
       "original_fastqtl_source_equals_historical_inherited_QTL":False,
       "reprocessed_2023_release6_original_bytes_hash_verified":False,
       "all_2225_candidate_genes_tested":False,
       "coloc_performed":False,"MR_performed":False,"scientific_status":"HOLD"}
    try:
        result["official_original_checksum_manifest"]=BASE.read_md5_list(md5_fetcher)
    except (OSError,ValueError) as exc:
        result["official_original_checksum_manifest"]={"gate":"HOLD_UNAVAILABLE",
                "failure_class":type(exc).__name__}
    if tabix is None:
        os.environ["CURL_CA_BUNDLE"]="/etc/ssl/certs/ca-certificates.crt"
        os.environ["SSL_CERT_FILE"]="/etc/ssl/certs/ca-certificates.crt"
        import pysam
        tabix=pysam
    handles=[]
    try:
        for url in (BASE.ORIGINAL,BASE.REPROCESSED):
            handles.append(tabix.TabixFile(url))
        headers=[]
        for source,hdl in zip((BASE.ORIGINAL,BASE.REPROCESSED),handles):
            header=BASE.normalize_header(hdl.header)
            if header is None:
                header=(header_fetcher(source) if header_fetcher else
                        BASE.original_bgzf_first_line(source))["header"]
            headers.append(header)
        result["source_headers"]={"original":headers[0],"reprocessed":headers[1]}
        result["header_schema_pass"]=all(BASE.REQUIRED.issubset(set(h)) for h in headers)
        if not result["header_schema_pass"]:
            result["fullcis_source_gate"]="HOLD_SOURCE_SCHEMA"
            return result
        chrom_results={}
        for region in FULL_WINDOWS:
            chrom=region[0]
            current=[]
            for h,hd in zip(handles,headers):
                time.sleep(wait)
                current.append(read_source_rows(h,hd,region))
            x,y=current
            chrom_results[chrom]={
                "region":f"{region[0]}:{region[1]}-{region[2]}",
                "original_source_scan":{k:v for k,v in x.items() if k!="rows"},
                "reprocessed_source_scan":{k:v for k,v in y.items() if k!="rows"},
                "four_genes":summarize_full_cis(x,y),
            }
        result["chromosomes"]=chrom_results
        genes={}
        for gene in BASE.TARGET_GENES:
            source="4" if gene in ("FGF5","C4orf22") else "10"
            genes[gene]=chrom_results[source]["four_genes"][gene]
        result["four_gene_fullcis_screen"]=genes
        result["total_original_variant_gene_variant_records"]=sum(
            v["original_unique_tested_variant"] for v in genes.values())
        result["total_reprocessed_variant_gene_variant_records"]=sum(
            v["reprocessed_unique_tested_variant"] for v in genes.values())
        result["total_exact_gene_variant_pairs"]=sum(
            v["matched_same_trait_and_exact_GRCh38_REF_ALT_pairs"]
            for v in genes.values())
        result["expected_QTD_fullcis_reprocessed_gene_counts_equal"]=all(
            g["reprocessed_variant_denominator_matches_prior_QTD_fullcis_audit"] for g in genes.values())
        result["all_source_windows_complete"]=all(
            chrom_results[c]["original_source_scan"]["status"]=="SCREEN_VALID"
            and chrom_results[c]["reprocessed_source_scan"]["status"]=="SCREEN_VALID"
            for c in chrom_results)
        result["fullcis_source_gate"]=(
            "FOUR_GENE_FULL_CIS_SOURCE_COMPARISON_SCREEN_PASS_NOT_COLOC"
            if result["total_exact_gene_variant_pairs"]>0 and
               result["expected_QTD_fullcis_reprocessed_gene_counts_equal"] and
               result["all_source_windows_complete"]
            else "HOLD_FULL_CIS_DENOMINATOR_SOURCE_OR_ACCESS"
        )
    except (OSError,ValueError,KeyError,TypeError,zlib.error) as exc:
        result["fullcis_source_gate"]="HOLD_SOURCE_CONNECTION_HEADER_OR_ANALYSIS"
        result["error_class"]=type(exc).__name__
        result["error_message_bounded"]=str(exc)[-160:]
    finally:
        for h in handles:
            h.close()
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    result=run()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("IS_GTEX_FOUR_GENE_FULL_CIS="+json.dumps(result,sort_keys=True,separators=(",",":")))

if __name__=="__main__":
    main()
