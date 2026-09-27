#!/usr/bin/env python3
from pathlib import Path
import csv
import json
import os
import re
import time
import urllib.request

ROOT=Path(os.environ.get("IS_ANALYSIS_ROOT","/srv/is-analysis"))
INST=ROOT/"results/metabolic_resilience/stage3_full_pgwas/instruments/STAGE3C2_FROZEN_INSTRUMENTS_ALL.tsv"
AUDIT=ROOT/"results/metabolic_resilience/stage3_full_pgwas/audit"
OUT=AUDIT/"STAGE3C4_DBSNP_RSID_ANNOTATION.tsv"
SUMMARY=AUDIT/"STAGE3C4_DBSNP_RSID_ANNOTATION.json"
RAW=AUDIT/"STAGE3C4_DBSNP_RAW"
RAW.mkdir(parents=True,exist_ok=True)

ASSEMBLY="GCF_000001405.25"
API=("https://api.ncbi.nlm.nih.gov/variation/v0/vcf/file/set_rsids"
     f"?assembly={ASSEMBLY}")
BATCH_SIZE=50
RS_RE=re.compile(r"^rs[0-9]+$")

if not INST.exists():
    print("[FAIL] frozen instrument file missing")
    raise SystemExit(2)

with INST.open("r",encoding="utf-8",newline="") as f:
    rows=list(csv.DictReader(f,delimiter="\t"))

unique={}
for r in rows:
    vid=r["variant_id"]
    if vid in unique:
        continue
    unique[vid]={
        "variant_id":vid,
        "gene":r["gene_symbol"],
        "chrom":str(r["chrom_hg19"]).replace("chr",""),
        "pos":int(r["pos_hg19"]),
        "ref":(r.get("ld_ref") or r.get("allele0") or r.get("other_allele")).upper(),
        "alt":(r.get("ld_alt") or r.get("allele1") or r.get("effect_allele")).upper(),
    }

variants=sorted(unique.values(),key=lambda x:(int(x["chrom"]),x["pos"],x["ref"],x["alt"]))
results={(x["chrom"],x["pos"],x["ref"],x["alt"]):{**x,"rsid":".","ncbi_id_raw":".","status":"NOT_QUERIED"} for x in variants}

def annotate_batch(batch,batch_no):
    body="\n".join("\t".join([x["chrom"],str(x["pos"]),".",x["ref"],x["alt"],".",".","."]) for x in batch)+"\n"
    for attempt in range(1,6):
        try:
            req=urllib.request.Request(
                API,data=body.encode("utf-8"),method="POST",
                headers={"Content-Type":"text/plain","Accept":"text/plain",
                         "User-Agent":"metabolic-resilience-master-thesis/1.0"}
            )
            with urllib.request.urlopen(req,timeout=120) as resp:
                text=resp.read().decode("utf-8",errors="replace")
            (RAW/f"batch_{batch_no:03d}.txt").write_text(text,encoding="utf-8")
            print(f"[PASS] batch={batch_no} n={len(batch)} bytes={len(text)}")
            return text
        except Exception as exc:
            print(f"[WARN] batch={batch_no} attempt={attempt} {exc!r}")
            time.sleep(min(2**attempt,20))
    return None

api_fail=0
for start in range(0,len(variants),BATCH_SIZE):
    batch=variants[start:start+BATCH_SIZE]
    batch_no=start//BATCH_SIZE+1
    text=annotate_batch(batch,batch_no)
    if text is None:
        api_fail+=1
        continue
    returned=0
    for line in text.splitlines():
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        p=line.split("\t")
        if len(p)<5:
            p=line.split()
        if len(p)<5:
            continue
        chrom=str(p[0]).replace("chr","")
        try: pos=int(p[1])
        except Exception: continue
        rid,ref,alt=p[2],p[3].upper(),p[4].upper()
        target=results.get((chrom,pos,ref,alt))
        if target is None:
            target=results.get((chrom,pos,alt,ref))
        if target is None:
            continue
        returned+=1
        rsids=sorted({x.strip() for x in rid.replace(";",",").split(",") if RS_RE.match(x.strip())})
        if len(rsids)==1:
            target["rsid"]=rsids[0]; target["status"]="PASS_RSID"
        elif len(rsids)>1:
            target["rsid"]=";".join(rsids); target["status"]="REVIEW_MULTIPLE_RSID"
        elif rid not in ("","."):
            target["status"]="NON_RSID_ID"
        else:
            target["status"]="NO_RSID_DBSNP"
        target["ncbi_id_raw"]=rid
    print(f"batch={batch_no} returned_rows={returned}")

out_rows=sorted(results.values(),key=lambda x:(int(x["chrom"]),x["pos"],x["ref"],x["alt"]))
fields=["variant_id","gene","chrom","pos","ref","alt","rsid","ncbi_id_raw","status"]
with OUT.open("w",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
    w.writeheader(); w.writerows(out_rows)

status_counts={}
for r in out_rows:
    status_counts[r["status"]]=status_counts.get(r["status"],0)+1
n=len(out_rows); n_rsid=sum(r["status"]=="PASS_RSID" for r in out_rows)
payload={
    "assembly":ASSEMBLY,
    "source":"NCBI Variation Services / dbSNP",
    "instrument_rows_original":len(rows),
    "unique_variants":n,
    "rsid_recovered":n_rsid,
    "rsid_recovery_pct":round(100*n_rsid/n,3) if n else 0,
    "multiple_rsid":sum(r["status"]=="REVIEW_MULTIPLE_RSID" for r in out_rows),
    "no_rsid_dbsnp":sum(r["status"]=="NO_RSID_DBSNP" for r in out_rows),
    "not_queried":sum(r["status"]=="NOT_QUERIED" for r in out_rows),
    "api_failed_batches":api_fail,
    "status_counts":status_counts,
    "canonical_variant_key":"GRCh37 chr:pos:REF:ALT",
    "policy":"C2 instrument definition remains frozen. rsID is secondary annotation only."
}
SUMMARY.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
print(json.dumps(payload,indent=2))
raise SystemExit(0 if n>0 and n_rsid==n and api_fail==0 else 2)
