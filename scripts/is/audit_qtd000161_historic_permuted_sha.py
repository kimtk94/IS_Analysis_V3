#!/usr/bin/env python3
"""Verify original 2026-10-11 cached QTD000161 permutation gzip SHA vs current EBI.

Prior Google Drive evidence (IS_PRIORITY4_EBI_FULL_CIS_RECONSTRUCTION_MANIFEST.json)
contains original gzip SHA and four p_beta/n_variants values. File is
small (<10MB) and can be downloaded within explicit byte and URL limits.
Even exact gzip identity here does NOT prove all.tsv.gz original-byte
identity or original GTEx FastQTL equivalence.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import urllib.parse
import urllib.request
from pathlib import Path

URL="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000161/QTD000161.permuted.tsv.gz"
HISTORICAL_SHA256="b5de9603b25b7f65f5532442bb5cc90f5d546421387363633967b540b4424bbe"
HISTORICAL_GENES={
 "ENSG00000138675":("FGF5",5511,0.0278074),
 "ENSG00000197826":("C4orf22",5490,0.000591983),
 "ENSG00000138172":("CALHM2",4868,0.0270452),
 "ENSG00000107954":("NEURL",4820,0.40605),
}
MAX_COMPRESSED=10_000_000
MAX_UNCOMPRESSED=70_000_000

def fetch_permutation(fetcher=None):
    if fetcher is not None:
        contents=fetcher(URL)
    else:
        req=urllib.request.Request(URL,headers={"User-Agent":"IS-QTD000161-Permutation-Release-Provenance/1.0"})
        with urllib.request.urlopen(req,timeout=40) as f:
            parsed=urllib.parse.urlsplit(f.url)
            if parsed.scheme!="https" or parsed.hostname!="ftp.ebi.ac.uk":
                raise ValueError("Unapproved official permutation redirect")
            contents=f.read(MAX_COMPRESSED+1)
    if len(contents)>MAX_COMPRESSED or not contents.startswith(b"\x1f\x8b"):
        raise ValueError("Unbounded or unexpected EBI permutation source")
    return contents

def inspect(raw):
    actual_hash=hashlib.sha256(raw).hexdigest()
    text=gzip.decompress(raw)
    if len(text)>MAX_UNCOMPRESSED:raise ValueError("Decompressed source exceeds bounded length")
    reader=csv.DictReader(io.StringIO(text.decode("utf-8-sig")),delimiter="\t")
    required={"molecular_trait_id","n_variants","p_beta"}
    if not required.issubset(reader.fieldnames or []):
        raise ValueError("Source permutation fields missing")
    genes={}
    source_total=0
    for row in reader:
        source_total+=1
        stable=row["molecular_trait_id"].split(".")[0]
        if stable not in HISTORICAL_GENES:continue
        symbol,count,p_beta=HISTORICAL_GENES[stable]
        if symbol in genes:raise ValueError("Duplicated original QTD permutation gene")
        actual=int(row["n_variants"])
        p=float(row["p_beta"])
        if actual<0 or not math.isfinite(p) or not 0<=p<=1:
            raise ValueError("Invalid permutation statistic")
        genes[symbol]={"stable_gene_id":stable,"historical_n_variants":count,
             "live_n_variants":actual,"historical_p_beta":p_beta,"live_p_beta":p,
             "n_variants_equal":count==actual,
             "p_beta_close":math.isclose(p,p_beta,rel_tol=0.00003,abs_tol=1e-8)}
    complete=set(genes)=={v[0] for v in HISTORICAL_GENES.values()}
    stats_match=complete and all(x["n_variants_equal"] and x["p_beta_close"] for x in genes.values())
    whole_sha_pass=actual_hash==HISTORICAL_SHA256
    return {
        "schema":"IS_QTD000161_permuted_gzip_20261011_vs_live_original_v1",
        "official_original_url":URL,
        "actual_compressed_source_bytes":len(raw),
        "current_EBI_source_whole_gzip_SHA256":actual_hash,
        "historic_20261011_cached_permutation_gzip_SHA256":HISTORICAL_SHA256,
        "whole_permuted_gzip_same_bytes_as_20261011_archive":whole_sha_pass,
        "live_source_total_permutation_gene_rows":source_total,
        "four_gene_permutation_comparison":genes,
        "all_four_gene_permutation_pbeta_n_variants_match_history":stats_match,
        "permutation_source_gate":(
            "HISTORIC_PERMUTED_GZIP_EXACT_SHA_MATCH_AND_FOUR_GENE_STATS_PASS"
            if whole_sha_pass and stats_match else
            "HOLD_PERMUTED_GZIP_SHA_DIFFERENT_OR_STATS_CHANGED"),
        "QTD_all_full_nominal_original_2023_gzip_sha256_authenticated":False,
        "is_GTEx_v8_original_FastQTL_statistics":False,
        "GWAS_ancestry_LD_N_approved":False,
        "MR_performed":False,"coloc_performed":False,"scientific_status":"HOLD",
    }

def run(fetcher=None):
    try:
        return inspect(fetch_permutation(fetcher))
    except (OSError,ValueError,UnicodeError,EOFError,gzip.BadGzipFile) as e:
        return {"schema":"IS_QTD000161_permuted_gzip_20261011_vs_live_original_v1",
            "permutation_source_gate":"HOLD_OFFICIAL_PERMUTATION_SOURCE_UNAVAILABLE",
            "failure_class":type(e).__name__,
            "failure_detail_bounded":str(e)[-160:],
            "scientific_status":"HOLD","coloc_performed":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    x=p.parse_args()
    r=run()
    x.output.parent.mkdir(parents=True,exist_ok=True)
    x.output.write_text(json.dumps(r,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("IS_QTD_PERM_HISTORY="+json.dumps(r,sort_keys=True,separators=(",",":")))
if __name__=="__main__":
    main()
