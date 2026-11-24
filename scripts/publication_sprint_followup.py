#!/usr/bin/env python3
"""Publication follow-up: public aggregate analysis only; no raw source uploads."""
import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
import math
import re
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path

from publication_sprint import POOL_SHA, read_gene_table, same_sign, sha

REACTOME = "https://reactome.org/download/current/ReactomePathways.gmt.zip"
SKIN_BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90094001-GCST90095000/"
ORAL_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103322/suppl/GSE103322_HNSCC_all_data.txt.gz"
ORAL_SHA = "51a5a3bb4f3ff8e326655ec613a10c9f895f9d5c943ce01c815303868b6e3a0b"
FIELDS = {
    "chromosome": ("chromosome","chrom","chr"),
    "position": ("base_pair_location","position","pos","bp"),
    "effect_allele": ("effect_allele","effectallele","ea","a1","allele1"),
    "other_allele": ("other_allele","otherallele","nea","a2","allele2","non_effect_allele"),
    "beta": ("beta","effect","effect_size"),
    "se": ("standard_error","se","stderr","standarderror"),
    "p": ("p_value","pvalue","p","pval"),
}
DNA = set("ACGT")

def transfer(url, path, limit):
    request = urllib.request.Request(url, headers={"User-Agent": "publication-sprint-aggregate/1.0", "Accept-Encoding": "identity"})
    digest = hashlib.md5(usedforsecurity=False)
    size = 0
    with urllib.request.urlopen(request, timeout=60) as response, path.open("wb") as out:
        if response.status != 200:
            raise ValueError("Expected complete HTTP response")
        expected = response.headers.get("Content-Length")
        if expected is not None and int(expected) > limit:
            raise ValueError("File exceeds transfer limit")
        while True:
            block = response.read(1024*1024)
            if not block:
                break
            size += len(block)
            if size > limit:
                raise ValueError("File exceeds transfer limit")
            digest.update(block)
            out.write(block)
        if expected is not None and size != int(expected):
            raise ValueError("Incomplete transfer")
    return {"bytes": size, "md5": digest.hexdigest(), "sha256": sha(path)}

def jsonout(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+"\n", encoding="utf8")
    print(json.dumps(obj, indent=2, allow_nan=False), flush=True)

def tsvout(path, rows, fields):
    with path.open("w", newline="", encoding="utf8") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

def bh(pvals):
    if any(not math.isfinite(p) or not 0 <= p <= 1 for p in pvals):
        raise ValueError("Invalid p value")
    n = len(pvals)
    adjusted = [1.0]*n
    running = 1.0
    for rank, i in reversed(list(enumerate(sorted(range(n), key=lambda i: pvals[i]), start=1))):
        running = min(running, n*pvals[i]/rank)
        adjusted[i] = running
    return adjusted

def hypergeom_tail(N, K, n, x):
    if not (0 <= K <= N and 0 <= n <= N and max(0,n-(N-K)) <= x <= min(K,n)):
        raise ValueError("Invalid contingency table")
    if x == max(0,n-(N-K)):
        return 1.0
    def lc(a,b):
        return math.lgamma(a+1)-math.lgamma(b+1)-math.lgamma(a-b+1)
    normal = lc(N,n)
    return min(1.0, math.fsum(math.exp(lc(K,j)+lc(N-K,n-j)-normal) for j in range(x,min(K,n)+1)))

def compare(donor, pool, external):
    shared = sorted(set(donor)&set(pool)&set(external))
    rows = []
    for g in shared:
        a,b,c = donor[g],pool[g],external[g]
        if c["fdr"] is None:
            continue
        stable = a[1] < .05 and b[1] < .05 and same_sign(a[0],b[0])
        coherent = same_sign(a[0],c["log2fc"])
        rows.append({"gene_id":g, "donor_symbol":a[2], "external_symbol":c["symbol"],
                     "donor_log2fc":a[0], "donor_fdr":a[1], "pool_log2fc":b[0],
                     "pool_fdr":b[1], "external_log2fc":c["log2fc"], "external_fdr":c["fdr"],
                     "strong_concordant":stable and min(abs(a[0]),abs(b[0]))>=1 and c["fdr"]<.05 and coherent,
                     "sign_discordant":stable and c["fdr"]<.05 and a[0]*c["log2fc"]<0})
    return len(shared), rows

def parse_gmt(path):
    out = {}
    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if n.endswith(".gmt")]
        if len(names) != 1:
            raise ValueError("Reactome archive must contain exactly one GMT")
        if z.getinfo(names[0]).file_size > 30*1024*1024:
            raise ValueError("Reactome GMT too large")
        for line in z.read(names[0]).decode("utf8").splitlines():
            a = line.split("\t")
            if len(a)<3 or not re.fullmatch(r"R-HSA-\d+",a[1]):
                continue
            if a[1] in out:
                raise ValueError("Duplicate pathway ID")
            out[a[1]] = (a[0], {s.strip().upper() for s in a[2:] if s.strip()})
    if len(out)<500:
        raise ValueError("Unexpected Reactome human pathway count")
    return out

def pathway_enrichment(rows, pathways):
    # Use unique, exactly cross-study-matching symbols; no alias imputation.
    counts = collections.Counter(r["donor_symbol"].strip().upper() for r in rows if r["donor_symbol"].strip())
    external_counts = collections.Counter(r["external_symbol"].strip().upper() for r in rows if r["external_symbol"].strip())
    eligible = {r["gene_id"]:r["donor_symbol"].strip().upper() for r in rows
                if r["donor_symbol"].strip() and r["donor_symbol"].strip().upper()==r["external_symbol"].strip().upper()
                and counts[r["donor_symbol"].strip().upper()]==1 and external_counts[r["external_symbol"].strip().upper()]==1}
    annotated = set().union(*(g for _,g in pathways.values()))
    universe = set(eligible.values()) & annotated
    selected = {eligible[r["gene_id"]] for r in rows if r["strong_concordant"] and r["gene_id"] in eligible} & universe
    if not selected or not universe:
        raise ValueError("No annotated candidate/background")
    result=[]
    for pid,(name,genes) in sorted(pathways.items()):
        genes = genes & universe
        if 10 <= len(genes) <= 500:
            hit = genes & selected
            p = hypergeom_tail(len(universe),len(genes),len(selected),len(hit))
            result.append({"pathway_id":pid,"pathway":name,"background_genes":len(genes),
                           "candidate_overlap":len(hit),"p":p,"overlap_symbols":";".join(sorted(hit))})
    if not result:
        raise ValueError("No eligible Reactome terms")
    for r,fdr in zip(result,bh([r["p"] for r in result])):
        r["fdr_bh"] = fdr
    result.sort(key=lambda r:(r["fdr_bh"],r["p"],r["pathway_id"]))
    return {"unique_matched_symbols":len(eligible),"annotated_background":len(universe),
            "annotated_candidates":len(selected),"tested_pathways":len(result),
            "fdr05_pathways":sum(r["fdr_bh"]<.05 for r in result)},result

def ocular(args):
    if sha(args.pool_zip)!=POOL_SHA:
        raise ValueError("P5d ZIP hash mismatch")
    sys.path.insert(0,str(args.ocular_source/"scripts/ocular"))
    from ocular_p5j_official_TM_overlap import download_zip,parse_official,TM_ALL
    with zipfile.ZipFile(args.pool_zip) as z:
        donor=read_gene_table(z,"p5_all_gene_results.tsv.gz","logFC","FDR",16050)
        pool=read_gene_table(z,"p5d_poolunit_all_gene_results.tsv.gz","pool_logFC","pool_FDR",21081)
    with tempfile.TemporaryDirectory() as d:
        d=Path(d)
        supp=d/"tm.zip"
        download_zip(supp)
        external=parse_official(supp)[TM_ALL]
        shared,rows=compare(donor,pool,external)
        strong=[r for r in rows if r["strong_concordant"]]
        discordant=[r for r in rows if r["sign_discordant"]]
        if shared!=12251 or len(strong)!=241 or len(discordant)!=64:
            raise ValueError("Published comparison differs from previous run; review source version")
        out=args.out.parent
        out.mkdir(parents=True,exist_ok=True)
        tsvout(out/"strong_241.tsv",strong,list(rows[0]))
        tsvout(out/"discordant_64.tsv",discordant,list(rows[0]))
        summary={"status":"DESCRIPTIVE_GENE_EVIDENCE_AND_EXPLORATORY_PATHWAYS",
                 "common_gene_ids":shared,"external_fdr_available":len(rows),
                 "strong_concordant":len(strong),"discordant":len(discordant),
                 "pool_zip_sha256":POOL_SHA,"external_supplement_sha256":sha(supp),
                 "pathway_interpretation":"One exploratory overrepresentation family, 10-500 annotated background genes per term; BH includes zero-overlap terms. Gene correlation, expression/selection bias and Reactome pathway overlap prevent confirmatory causal conclusions. Donor and pool are not independent studies."}
        try:
            gmt=d/"reactome.zip"
            source=transfer(REACTOME,gmt,5*1024*1024)
            meta,result=pathway_enrichment(rows,parse_gmt(gmt))
            tsvout(out/"reactome_all_terms.tsv",result,list(result[0]))
            summary.update({"reactome_source":REACTOME,"reactome_archive":source,**meta,
                            "top_pathways":result[:8]})
        except (OSError,ValueError,zipfile.BadZipFile) as e:
            summary.update({"status":"GENE_EVIDENCE_COMPLETE_PATHWAYS_UNRESOLVED","pathway_error_type":type(e).__name__,"pathway_error":str(e)[:180]})
        jsonout(args.out,summary)

def normalize(s):
    return re.sub(r"[^a-z0-9]+","_",s.lower().strip()).strip("_")

def skin_schema(header):
    index={}
    names=[normalize(x) for x in header]
    for role,aliases in FIELDS.items():
        found=[i for i,name in enumerate(names) if name in aliases]
        if len(found)!=1:
            raise ValueError("Missing/ambiguous required role: "+role)
        index[role]=found[0]
    return index

def scan_gwas(handle):
    reader=csv.reader(handle,delimiter="\t")
    header=next(reader)
    index=skin_schema(header)
    n=0
    valid=0
    invalid=collections.Counter()
    significant=0
    zero_p=0
    chrom_counts=collections.Counter()
    minima={}
    for row in reader:
        n+=1
        if n>15_000_000:
            raise ValueError("Row limit")
        if len(row)!=len(header):
            invalid["width"]+=1
            continue
        try:
            chrom=int(re.sub(r"^chr","",row[index["chromosome"]],flags=re.I))
            pos=int(row[index["position"]])
            if not 1<=chrom<=22 or pos<1:
                invalid["non_autosomal_or_position"]+=1
                continue
            a,b=(row[index[k]].upper() for k in ("effect_allele","other_allele"))
            if a not in DNA or b not in DNA or a==b:
                invalid["not_biallelic_snv"]+=1
                continue
            beta,se,p=(float(row[index[k]]) for k in ("beta","se","p"))
            if not all(math.isfinite(x) for x in (beta,se,p)) or se<=0 or not 0<=p<=1:
                invalid["numeric"]+=1
                continue
        except (ValueError,IndexError):
            invalid["parse"]+=1
            continue
        valid+=1
        chrom_counts[chrom]+=1
        if p==0:
            zero_p+=1
            continue
        if p<5e-8:
            significant+=1
        if chrom not in minima or p<minima[chrom]["p"]:
            minima[chrom]={"chromosome":chrom,"position":pos,"p":p}
    return {"header":header,"schema":{k:header[i] for k,i in index.items()},
            "data_rows":n,"valid_autosomal_snv_rows":valid,"excluded_rows":dict(invalid),
            "positive_p_lt_5e8_rows":significant,"zero_p_rows_not_counted_significant":zero_p,
            "autosomes_with_valid_rows":len(chrom_counts),"minimum_positive_p_by_chromosome":[minima[k] for k in sorted(minima)],
            "interpretation":"Rows, not unique independent loci. No LD clumping, MR, effect-allele harmonization across sources or phenotype equivalence test."}

def skin(args):
    accession=args.accession
    if accession not in ("GCST90094903","GCST90094904"):
        raise ValueError("Unapproved accession")
    name=accession+"_buildGRCh37.tsv"
    base=SKIN_BASE+accession+"/"
    result={"status":"HOLD_SOURCE_NOT_VERIFIED","accession":accession,"url":base+name}
    try:
        with tempfile.TemporaryDirectory() as d:
            d=Path(d)
            manifest=d/"md5sum.txt"
            transfer(base+"md5sum.txt",manifest,65536)
            hashes=[]
            for line in manifest.read_text(encoding="utf8").splitlines():
                m=re.fullmatch(r"\s*([a-fA-F0-9]{32})\s+\*?(.+?)\s*",line)
                if m and m[2].rsplit("/",1)[-1]==name:
                    hashes.append(m[1].lower())
            if len(hashes)!=1:
                raise ValueError("Exact provider MD5 not found")
            result["expected_md5"]=hashes[0]
            data=d/name
            result["transfer"]=transfer(base+name,data,1500*1024*1024)
            if result["transfer"]["md5"]!=hashes[0]:
                raise ValueError("Provider MD5 mismatch")
            result["status"]="FULL_SOURCE_MD5_VERIFIED_SCHEMA_PENDING"
            with data.open("rt",encoding="utf-8-sig",newline="") as f:
                first=f.readline()
                result["observed_header"]=first.strip().split("\t")[:80]
                f.seek(0)
                result.update(scan_gwas(f))
            result["status"]="FULL_SOURCE_SCANNED_NO_MR"
    except (OSError,ValueError,StopIteration,UnicodeError,csv.Error) as e:
        result["error_type"]=type(e).__name__
        result["error"]=str(e)[:200]
    jsonout(args.out,result)

def oral_metadata(handle):
    reader=csv.reader(handle,delimiter="\t")
    first=next(reader)
    if len(first)>1 and first[-1]=="":
        first=first[:-1]
    n=len(first)-1
    if n!=5902 or len(set(first[1:]))!=n:
        raise ValueError("Unexpected public cell columns")
    required={"processed_by_maxima_enzyme","lymph_node","classified_as_cancer_cell","classified_as_non_cancer_cells","non_cancer_cell_type"}
    rows={}
    for _ in range(5):
        r=next(reader)
        if len(r)==n+2 and r[-1]=="":
            r=r[:-1]
        if len(r)!=n+1 or normalize(r[0]) not in required or normalize(r[0]) in rows:
            raise ValueError("Unexpected official metadata layout")
        rows[normalize(r[0])]=r[1:]
    if set(rows)!=required:
        raise ValueError("Missing official metadata")
    for field in required-{"non_cancer_cell_type"}:
        if set(rows[field])-{"0","1"}:
            raise ValueError("Unexpected classification values")
    site=rows["lymph_node"]
    cancer=rows["classified_as_cancer_cell"]
    non=rows["classified_as_non_cancer_cells"]
    enzyme=rows["processed_by_maxima_enzyme"]
    cross=collections.Counter((s,c,t) for s,c,t in zip(site,cancer,non))
    return {"cell_columns":n,"site_classification_cell_counts":[{"lymph_node":s,"cancer":c,"noncancer":t,"cells":v} for (s,c,t),v in sorted(cross.items())],
            "enzyme_by_site_cell_counts":[{"lymph_node":s,"maxima":m,"cells":v} for (s,m),v in sorted(collections.Counter(zip(site,enzyme)).items())],
            "noncancer_type_cell_counts":dict(collections.Counter(rows["non_cancer_cell_type"])),
            "patient_mapping_verified":False,"patient_level_inference_performed":False,
            "interpretation":"Public labels only. These are cell counts, not independent patient observations. No inference of patient identities from column names; no hypothesis p value."}

def oral(args):
    result={"status":"HOLD_SOURCE_OR_METADATA"}
    try:
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"oral.gz"
            result["transfer"]=transfer(ORAL_URL,path,120*1024*1024)
            if result["transfer"]["sha256"]!=ORAL_SHA:
                raise ValueError("Public matrix differs from previous hash")
            with gzip.open(path,"rt",encoding="utf8",newline="") as h:
                result.update(oral_metadata(h))
            result["status"]="EXPLICIT_SITE_AND_CLASSIFICATION_CELL_COUNTS"
    except (OSError,ValueError,StopIteration,UnicodeError,csv.Error) as e:
        result.update({"error_type":type(e).__name__,"error":str(e)[:180]})
    jsonout(args.out,result)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("analysis",choices=["ocular","skin","oral"])
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--pool-zip",type=Path)
    p.add_argument("--ocular-source",type=Path)
    p.add_argument("--accession")
    args=p.parse_args()
    globals()[args.analysis](args)

if __name__=="__main__":
    main()
