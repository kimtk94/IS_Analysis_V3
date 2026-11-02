#!/usr/bin/env python3
"""Eight bounded, distinct scientific-evidence tasks, executed ONLY in GitHub Actions.

Public-source observations are not scientific GO. Do not export individual patient
identifiers, raw counts, private URLs or complete expression matrices.
"""
from __future__ import annotations
import argparse
import collections
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import math
import os
import pathlib
import re
import struct
import urllib.parse
import urllib.request
import zipfile
import zlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "results/research_operations/20261103"
TASKS = {
    "is": "Inspect a real GIGASTROKE EAS AIS summary-statistics source header and allele/effect schema",
    "glaucoma": "Independently recalculate BH-FDR for the archived 57-pair edgeR gene-level output",
    "hcc": "Independently verify official PDC tumor/normal case pairing in both proteomics studies",
    "radiomics": "Check real public NIfTI voxel-spacing headers or verify 14-case CI result accessibility",
    "ocular": "Audit all 46 ENA runs for two FASTQ paths, byte sizes and published MD5 syntax",
    "muscle": "Resolve nine independently archived Yang 2024 RT rsIDs against GTEx v10",
    "oral": "Verify GSE103322 expression resource accessibility and clinical donor linkage prerequisites",
    "skin": "Validate the bounded actual perceived-age GWAS header/allele schema against official readme",
}
NEXT = {
    "is": "Complete source file SHA256, SNP effects and EUR/EAS allele/LD harmonization",
    "glaucoma": "Re-fit leave-one-pool-out QL on donor-level integer counts and independently replicate",
    "hcc": "Prove 330 quant columns uniquely map to official aliquots and RNA cases",
    "radiomics": "Measure registered PV NIfTI spacing for 14 cases and MRI-to-RNA patient join",
    "ocular": "Resolve 1h/6h/48h, strain and donor/timepoint pairing before raw-count analysis",
    "muscle": "Acquire all 20 RT lead loci and ancestry-matched LD/eQTL/coloc",
    "oral": "Verify donor identities in real expression/clinical matrices and independent cohort",
    "skin": "Whole-file checksum, build, full SNP QC and MR pleiotropy/coloc",
}
ALLOWED = {"ftp.ebi.ac.uk", "www.ebi.ac.uk", "pdc.cancer.gov", "zenodo.org",
           "www.ncbi.nlm.nih.gov", "bioconductor.org", "www.bioconductor.org",
           "gtexportal.org", "data.bris.ac.uk"}
class Hold(Exception): pass
class IntegrityFail(Exception): pass

class SafeRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        host = urllib.parse.urlsplit(newurl).hostname
        if host not in ALLOWED or urllib.parse.urlsplit(newurl).scheme != "https":
            raise Hold("SOURCE_REDIRECT_NOT_ALLOWED")
        return super().redirect_request(req, fp, code, msg, headers, newurl)

OPENER = urllib.request.build_opener(SafeRedirect)

def fetch(url, cap=1000000, headers=None, method="GET", data=None):
    u = urllib.parse.urlsplit(url)
    if u.scheme != "https" or u.hostname not in ALLOWED:
        raise Hold("UNVERIFIED_SOURCE_HOST")
    req = urllib.request.Request(url, data=data, method=method,
        headers={"User-Agent":"MASTER-DEGREE-20261103-CI-evidence/1.0",
                 "Accept":"*/*", **(headers or {})})
    with OPENER.open(req, timeout=25) as response:
        b = response.read(cap+1) if method != "HEAD" else b""
        if len(b) > cap: raise Hold("SOURCE_EXCEEDS_BOUNDED_LIMIT")
        return b, {"http":response.status,"bytes":len(b),
                   "sha256":hashlib.sha256(b).hexdigest(),
                   "content_range":response.headers.get("Content-Range"),
                   "content_length":response.headers.get("Content-Length")}

def is_task():
    # Directory path is derived from the accession block convention, not guessed file names.
    url = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
           "GCST90104001-GCST90105000/GCST90104545/")
    page, p = fetch(url, 500000)
    hrefs = re.findall(r'href=["\x27]([^"\x27]+)["\x27]', page.decode("utf-8","replace"), re.I)
    files = [x for x in hrefs if re.fullmatch(r"[\w.+-]+\.(?:tsv|txt)(?:\.gz)?",x,re.I)]
    if not files: raise Hold("EAS_AIS_OFFICIAL_SUMSTATS_FILE_NOT_DISCOVERED")
    # Prefer harmonised files if actually present.
    files.sort(key=lambda x:(not ("harmon" in x.lower()),x))
    file_url = urllib.parse.urljoin(url,files[0])
    part, m = fetch(file_url, 65536, {"Range":"bytes=0-65535","Accept-Encoding":"identity"})
    if m["http"] != 206 and m["content_range"] is None:
        raise Hold("SUMSTATS_RANGE_NOT_SUPPORTED")
    if files[0].endswith(".gz"):
        try: head=zlib.decompressobj(16+zlib.MAX_WBITS).decompress(part, 200000).decode("utf-8","replace")
        except zlib.error: raise Hold("SUMSTATS_GZIP_PREFIX_INVALID")
    else: head=part.decode("utf-8","replace")
    lines=[x for x in head.splitlines() if x and not x.startswith("##")]
    if not lines: raise Hold("SUMSTATS_HEADER_NOT_READABLE")
    cols=lines[0].lstrip("#").split("\t")
    if len(cols)<5 or len(cols)!=len(set(cols)): raise IntegrityFail("SUMSTATS_DUPLICATE_OR_SHORT_HEADER")
    lc=[c.lower() for c in cols]
    allele=any("effect_allele" in c for c in lc) and any("other_allele" in c for c in lc)
    effect=any(c in lc for c in ("beta","hm_beta","odds_ratio","hm_odds_ratio","effect_size"))
    if not allele or not effect: raise Hold("SUMSTATS_EFFECT_ALLELE_COLUMNS_NOT_IDENTIFIED")
    return "PASS_BOUNDED_ORIGINAL_SUMSTATS_HEADER", {
        "study":"GCST90104545","official_listing_sha256":p["sha256"],
        "files_discovered":len(files),"selected_source_file":files[0],
        "compressed_prefix_sha256":m["sha256"],"prefix_bytes":len(part),
        "header_column_count":len(cols),"allele_fields_present":allele,
        "effect_fields_present":effect,"whole_file_checksum_verified":False,
    }

def _bh(pvals):
    n=len(pvals); order=sorted(range(n),key=lambda i:pvals[i])
    q=[0.]*n; m=1.
    for k in range(n-1,-1,-1):
        i=order[k]; m=min(m,pvals[i]*n/(k+1)); q[i]=m
    return q

def glaucoma_task():
    folder=pathlib.Path("/tmp/research-p5-source")
    if not folder.exists(): raise Hold("ARCHIVED_P5_ARTIFACT_NOT_ACCESSIBLE_IN_CI")
    candidates=[]
    for p in folder.rglob("*"):
        if not p.is_file() or p.stat().st_size>25000000 or p.suffix.lower() not in (".tsv",".csv"):
            continue
        try:
            with p.open(encoding="utf-8-sig") as f:
                first=f.readline()
            sep="\t" if "\t" in first else ","
            cols=[x.strip().lower() for x in next(csv.reader([first],delimiter=sep))]
            if any(x in cols for x in ("pvalue","p_value","p.value","pval","p")) and any(
                    x in cols for x in ("fdr","padj","adj.p.val","qvalue")):
                candidates.append((p,sep,cols))
        except (OSError,UnicodeError): continue
    if not candidates:
        return "HOLD_P5_GENE_LEVEL_SCHEMA",{"archive_access":True,
            "candidate_tables":0,"scanned_files":sum(p.is_file() for p in folder.rglob("*"))}
    p,sep,cols=max(candidates,key=lambda x:x[0].stat().st_size)
    pi=next(cols.index(x) for x in ("pvalue","p_value","p.value","pval","p") if x in cols)
    qi=next(cols.index(x) for x in ("fdr","padj","adj.p.val","qvalue") if x in cols)
    values=[]
    with p.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.reader(f,delimiter=sep); next(reader)
        for row in reader:
            if len(row)<=max(pi,qi): continue
            try: a,b=float(row[pi]),float(row[qi])
            except (ValueError,TypeError): continue
            if not all(math.isfinite(v) and 0<=v<=1 for v in (a,b)):
                raise IntegrityFail("P5_NONFINITE_OR_INVALID_PVALUES")
            values.append((a,b))
    if len(values)<1000: raise Hold("P5_INSUFFICIENT_VALID_GENE_ROWS")
    recomputed=_bh([a for a,b in values])
    diffs=[abs(a-b) for a,(_,b) in zip(recomputed,values)]
    bad=sum(d>max(1e-5,1e-3*max(a,b)) for d,a,(_,b) in zip(diffs,recomputed,values))
    if bad>0: raise IntegrityFail("P5_BH_FDR_DISCREPANCY_"+str(bad))
    return "PASS_ARCHIVED_P5_BH_RECALCULATION",{
        "tested_genes":len(values),"bh_recalculated":len(values),
        "fdr_below_0_05":sum(q<0.05 for q in recomputed),
        "max_abs_fdr_delta":max(diffs),"source_table_sha256":hashlib.sha256(p.read_bytes()).hexdigest(),
        "donor_level_refit":False}

def _pdc(study):
    query='{ biospecimenPerStudy(pdc_study_id:"'+study+'",acceptDUA:true){aliquot_id sample_id case_id sample_type} }'
    body,_=fetch("https://pdc.cancer.gov/graphql",5000000,
      {"Content-Type":"application/json","Accept":"application/json"},
      method="POST",data=json.dumps({"query":query}).encode())
    obj=json.loads(body)
    if obj.get("errors"):raise Hold("PDC_GRAPHQL_SCHEMA_OR_ACCESS_ERROR")
    rows=obj.get("data",{}).get("biospecimenPerStudy")
    if not isinstance(rows,list) or not rows:raise Hold("PDC_EMPTY_BIOSPECIMEN")
    seen={};types=collections.defaultdict(set);sample_cases={}
    for r in rows:
        a,s,c,t=(r.get(x) for x in ("aliquot_id","sample_id","case_id","sample_type"))
        if not all(isinstance(v,str) and v for v in (a,s,c)):
            raise IntegrityFail("PDC_MISSING_PARENT_IDS")
        if a in seen and seen[a]!=(s,c):raise IntegrityFail("PDC_AMBIGUOUS_ALIQUOT_PARENT")
        if a in seen:raise IntegrityFail("PDC_DUPLICATE_ALIQUOT")
        if s in sample_cases and sample_cases[s]!=c:raise IntegrityFail("PDC_SAMPLE_MULTIPLE_CASES")
        seen[a]=(s,c);sample_cases[s]=c
        typ=str(t or "").lower()
        if "normal" in typ: types[c].add("N")
        elif "tumor" in typ or "tumour" in typ:types[c].add("T")
    canonical=sorted((a,s,c) for a,(s,c) in seen.items())
    return {"aliquots":len(seen),"cases":len(set(c for _,c in seen.values())),
      "paired_tumor_normal_cases":sum(v.issuperset({"T","N"}) for v in types.values()),
      "typed_cases":sum(bool(v) for v in types.values()),
      "mapping_sha256":hashlib.sha256(json.dumps(canonical,separators=(",",":")).encode()).hexdigest()}

def hcc_task():
    a,b=_pdc("PDC000198"),_pdc("PDC000199")
    if a["mapping_sha256"]!=b["mapping_sha256"]:
        raise IntegrityFail("PDC_CROSS_STUDY_PARENT_MAPPING_CONFLICT")
    if a["aliquots"]!=331 or a["cases"]!=171:
        raise Hold("PDC_DENOMINATOR_CHANGED_REQUIRES_REVIEW")
    if a["paired_tumor_normal_cases"]<1:raise Hold("PDC_SAMPLE_TYPE_CANNOT_CLASSIFY")
    return "PASS_OFFICIAL_CASE_PAIR_COHORT_AUDIT",{
        "proteome":a,"phosphoproteome":b,
        "quant_column_to_aliquot_verified":False,
        "rna_patient_join_verified":False}

def _nifti_header(raw,compressed):
    if compressed:
        try: raw=zlib.decompressobj(16+zlib.MAX_WBITS).decompress(raw,4096)
        except zlib.error:raise Hold("NIFTI_GZIP_PREFIX_INVALID")
    if len(raw)<348:raise Hold("NIFTI_HEADER_INCOMPLETE")
    endian="<" if struct.unpack("<I",raw[:4])[0]==348 else ">"
    if struct.unpack(endian+"I",raw[:4])[0]!=348:raise Hold("NIFTI_1_HEADER_NOT_FOUND")
    pixdim=struct.unpack(endian+"8f",raw[76:108])
    spacing=list(pixdim[1:4])
    if any(not math.isfinite(x) or x<=0 or x>100 for x in spacing):
        raise IntegrityFail("NIFTI_INVALID_VOXEL_SPACING")
    return spacing

def radiomics_task():
    folder=pathlib.Path("/tmp/research-radiomics-source")
    local=[p for p in folder.rglob("*") if p.is_file()] if folder.exists() else []
    # Source artifact is preferred but its schema cannot be assumed.
    source_artifact={"available":bool(local),"file_count":len(local),
                     "contains_tabular_result":any(p.suffix.lower() in (".csv",".tsv") for p in local)}
    obj=json.loads(fetch("https://zenodo.org/api/records/8179129",3000000)[0])
    files=obj.get("files") or []
    nifti=[f for f in files if str(f.get("key","")).lower().endswith((".nii",".nii.gz"))]
    if not nifti:
        return "HOLD_ZENODO_NO_DIRECT_NIFTI",{
            "zenodo_record_id":8179129,"files_listed":len(files),
            "direct_nifti_files":0,"prior_14_case_artifact":source_artifact}
    f=nifti[0]; url=(f.get("links") or {}).get("self")
    if not url:raise Hold("ZENODO_NIFTI_DOWNLOAD_LINK_MISSING")
    part,meta=fetch(url,65536,{"Range":"bytes=0-65535","Accept-Encoding":"identity"})
    if meta["http"]!=206 and not meta["content_range"]:raise Hold("ZENODO_NIFTI_RANGE_UNSUPPORTED")
    spacing=_nifti_header(part,str(f["key"]).endswith(".gz"))
    return "PASS_ONE_PUBLIC_NIFTI_SPACING_ONLY",{
        "zenodo_record_id":8179129,"nifti_file_count":len(nifti),
        "one_nifti_spacing_mm":spacing,"prefix_sha256":meta["sha256"],
        "prior_14_case_artifact":source_artifact,"registered_PV_14_verified":False}

def ocular_task():
    url=("https://www.ebi.ac.uk/ena/portal/api/filereport?"
         "accession=PRJNA1243143&result=read_run&format=tsv"
         "&fields=run_accession,sample_accession,fastq_ftp,fastq_md5,fastq_bytes")
    body,meta=fetch(url,2000000)
    rows=list(csv.DictReader(io.StringIO(body.decode("utf-8-sig")),delimiter="\t"))
    if len(rows)!=46:raise Hold("ENA_RUN_COUNT_NOT_46")
    runids=[r.get("run_accession") for r in rows]
    if len(set(runids))!=46:raise IntegrityFail("ENA_DUPLICATE_RUN")
    bad=0;unique=set();sizes=[];md5_count=0
    for r in rows:
        paths=(r.get("fastq_ftp") or "").split(";")
        md5s=(r.get("fastq_md5") or "").split(";")
        bts=(r.get("fastq_bytes") or "").split(";")
        if not (len(paths)==len(md5s)==len(bts)==2):bad+=1;continue
        for path,m,size in zip(paths,md5s,bts):
            if not re.fullmatch(r"[a-fA-F0-9]{32}",m) or not path.startswith("ftp.sra.ebi.ac.uk/"):
                bad+=1;continue
            try: sz=int(size)
            except ValueError:bad+=1;continue
            if sz<=0:bad+=1;continue
            if path in unique:raise IntegrityFail("ENA_REUSED_FASTQ_PATH")
            unique.add(path);sizes.append(sz);md5_count+=1
    return ("PASS_46_RUN_FASTQ_MANIFEST" if bad==0 else "HOLD_FASTQ_MANIFEST_INCOMPLETE"),{
        "runs":46,"fastq_paths_unique":len(unique),"md5_syntax_valid":md5_count,
        "invalid_or_incomplete_records":bad,"declared_total_bytes":sum(sizes),
        "ena_response_sha256":meta["sha256"],"actual_fastq_bytes_downloaded":False}

def muscle_task():
    p=ROOT/"config/muscle/stage1_seed_loci.tsv"
    if not p.exists():raise Hold("SOURCE_SEED_LOCI_MISSING")
    with p.open(encoding="utf-8-sig",newline="") as f:rows=list(csv.DictReader(f,delimiter="\t"))
    rt=[r for r in rows if r.get("cohort_type")=="RT"]
    rsids=sorted(set(r.get("rsid","") for r in rt))
    if len(rsids)!=9 or any(not re.fullmatch(r"rs\d+",x) for x in rsids):
        raise Hold("YANG_SEED_9_LOCI_NOT_REPRODUCIBLE")
    resolved=0;not_found=0;errors=0
    for rsid in rsids:
        url="https://gtexportal.org/api/v2/dataset/variant?"+urllib.parse.urlencode(
            {"snpId":rsid,"datasetId":"gtex_v10","itemsPerPage":2000})
        try:
            data=json.loads(fetch(url,1000000)[0])
            rec=data.get("data")
            if not isinstance(rec,list):errors+=1
            elif rec:resolved+=1
            else:not_found+=1
        except Exception:errors+=1
    return ("PASS_9_GTEX_VARIANT_RESOLUTION" if errors==0 else "HOLD_GTEX_PARTIAL_ACCESS"),{
        "source_seed_sha256":hashlib.sha256(p.read_bytes()).hexdigest(),
        "yang_RT_loci":9,"full_target_loci":20,"gtex_v10_resolved":resolved,
        "not_found":not_found,"api_errors":errors,"eQTL_verified":False}

def oral_task():
    url="https://bioconductor.org/packages/release/data/experiment/html/GSE103322.html"
    body,meta=fetch(url,2000000)
    text=body.decode("utf-8","replace")
    if "GSE103322" not in text:raise Hold("BIOCONDUCTOR_ACCESSION_NOT_PRESENT")
    hrefs=re.findall(r'href=["\x27]([^"\x27]+)["\x27]',text,re.I)
    src=[x for x in hrefs if "GSE103322_" in x and x.endswith(".tar.gz")]
    return ("PASS_BIOCONDUCTOR_EXPERIMENT_RESOURCE_LINK" if src else
            "HOLD_NO_DIRECT_EXPERIMENT_RESOURCE_LINK"),{
        "accession":"GSE103322","official_package_sha256":meta["sha256"],
        "direct_package_links":len(src),"raw_expression_downloaded":False,
        "donor_ids_verified":False,"patient_level_proteome_join":False}

def skin_task():
    readme,rm=fetch("https://data.bris.ac.uk/datasets/21crwsnj4xwjm2g4qi8chathha/readme.txt",500000)
    required=("ALLELE1","ALLELE0","A1FREQ","LOG_OR","LOG_SE")
    missing=[x for x in required if x not in readme.decode("utf-8","replace").upper()]
    if missing:raise IntegrityFail("GWAS_README_FIELDS_CHANGED")
    url=("https://data.bris.ac.uk/datasets/21crwsnj4xwjm2g4qi8chathha/"
         "perceived_age_summary_results.txt")
    part,meta=fetch(url,65536,{"Range":"bytes=0-65535","Accept-Encoding":"identity"})
    if meta["http"]!=206 and not meta["content_range"]:
        raise Hold("BRISTOL_RANGE_NOT_SUPPORTED")
    lines=part.decode("utf-8","replace").splitlines()
    if not lines:raise Hold("BRISTOL_HEADER_NOT_AVAILABLE")
    cols=re.split(r"\s+",lines[0].strip())
    if len(cols)!=len(set(cols)):raise IntegrityFail("BRISTOL_DUPLICATE_HEADER_COLUMNS")
    if not all(x in (c.upper() for c in cols) for x in required):
        raise Hold("BRISTOL_ACTUAL_FILE_HEADER_NOT_MATCH_README")
    return "PASS_BOUNDED_GWAS_HEADER_ALLELE_SCHEMA",{
        "readme_sha256":rm["sha256"],"file_prefix_sha256":meta["sha256"],
        "header_columns":len(cols),"partial_lines":len(lines),
        "whole_gwas_sha256_verified":False,"genome_build_verified":False}

RUNNERS={"is":is_task,"glaucoma":glaucoma_task,"hcc":hcc_task,
         "radiomics":radiomics_task,"ocular":ocular_task,"muscle":muscle_task,
         "oral":oral_task,"skin":skin_task}

def execute(track):
    if os.environ.get("GITHUB_ACTIONS")!="true":
        raise SystemExit("Research execution prohibited outside GitHub Actions")
    r={"date_kst":"2026-11-03","track":track,"selected_task":TASKS[track],
       "execution":"GITHUB_ACTIONS_ONLY","scientific_gate":"HOLD",
       "next_gate":NEXT[track],"ckd":"HOLD_NO_EXECUTION",
       "observed_utc":dt.datetime.now(dt.timezone.utc).isoformat()}
    try:
        status,evidence=RUNNERS[track]()
        r.update(probe_status=status,evidence=evidence)
    except IntegrityFail as exc:
        r.update(probe_status="FAIL_INTEGRITY",blocker=str(exc),evidence={})
    except Exception as exc:
        r.update(probe_status="HOLD_SOURCE_OR_INPUT",blocker=type(exc).__name__+": "+str(exc)[:130],evidence={})
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/(track+".json")).write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({k:r.get(k) for k in ("track","probe_status","scientific_gate","blocker")},ensure_ascii=False))
    return 1 if r["probe_status"].startswith("FAIL") else 0

if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--track",choices=TASKS,required=True)
    args=ap.parse_args();raise SystemExit(execute(args.track))
