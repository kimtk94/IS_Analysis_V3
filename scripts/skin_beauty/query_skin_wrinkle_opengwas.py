#!/usr/bin/env python3
from __future__ import annotations
import csv, gzip, json, os, re, time, urllib.error, urllib.request
from pathlib import Path

BASE="https://api.opengwas.io/api"
OUTCOMES={
 "ebi-a-GCST90094903":"facial_wrinkles_under_eye",
 "ebi-a-GCST90094904":"facial_wrinkles_crows_feet",
}
ROOT=Path(os.environ.get("ROOT","/srv/is-analysis/IS_Analysis_V3"))
OUT=ROOT/"results/skin_beauty/stage0_audit"
AUDIT=OUT/"STAGE0_SKIN_GWAS_AUDIT.json"
JWT=os.environ.get("OPENGWAS_JWT","").strip()

if not JWT:
    raise SystemExit("OPENGWAS_JWT is missing")
OUT.mkdir(parents=True, exist_ok=True)

def api(path, payload=None, retries=4):
    url=BASE+path
    data=None if payload is None else json.dumps(payload).encode()
    method="GET" if payload is None else "POST"
    headers={"Authorization":"Bearer "+JWT,"Accept":"application/json","X-Api-Source":"skin-beauty-mr/1.0"}
    if data is not None: headers["Content-Type"]="application/json"
    last=None
    for attempt in range(retries):
        req=urllib.request.Request(url,data=data,headers=headers,method=method)
        try:
            with urllib.request.urlopen(req,timeout=120) as r:
                raw=r.read().decode("utf-8",errors="replace")
                return json.loads(raw) if raw.strip() else {}, dict(r.headers)
        except urllib.error.HTTPError as e:
            body=e.read().decode("utf-8",errors="replace")
            last=RuntimeError(f"HTTP {e.code} {path}: {body[:1000]}")
            if e.code not in (429,500,502,503,504): break
            retry=e.headers.get("Retry-After")
            delay=min(15, int(retry) if retry and retry.isdigit() else 2*(attempt+1))
            time.sleep(delay)
        except Exception as e:
            last=e; time.sleep(2*(attempt+1))
    raise last

def rows_from(obj):
    if isinstance(obj,list): return obj
    if isinstance(obj,dict):
        for k in ("data","results","result"):
            if isinstance(obj.get(k),list): return obj[k]
        # gwasinfo may return keyed dict
        vals=[v for v in obj.values() if isinstance(v,dict)]
        if vals: return vals
    return []

def write_tsv(path, rows):
    rows=[r for r in rows if isinstance(r,dict)]
    cols=[]
    for r in rows:
        for k in r:
            if k not in cols: cols.append(k)
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=cols or ["status"],delimiter="\t",extrasaction="ignore")
        w.writeheader()
        for r in rows: w.writerow({k:(json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v) for k,v in r.items()})

def get_rsids(exposure, n=8):
    op=gzip.open if str(exposure).endswith(".gz") else open
    with op(exposure,"rt",encoding="utf-8",errors="replace") as f:
        first=f.readline(); delim="\t" if first.count("\t")>=first.count(",") else ","
        cols=first.rstrip("\r\n").split(delim); norm=[re.sub(r"[^a-z0-9]+","_",c.lower()).strip("_") for c in cols]
        choices=["rsid","rs_id","snp","variant","variant_id","variant_id_hg19"]
        idx=None
        for c in choices:
            if c in norm: idx=norm.index(c); break
        if idx is None: raise RuntimeError("No rsID/variant column found in exposure")
        out=[]
        for line in f:
            v=line.rstrip("\r\n").split(delim)[idx]
            if re.fullmatch(r"rs\d+",v) and v not in out:
                out.append(v)
                if len(out)>=n: break
        return out

print("===== OPEN GWAS JWT VALIDATION =====")
user,_=api("/user"); print("JWT valid; /user returned successfully")

print("===== OPENGWAS METADATA =====")
meta,_=api("/gwasinfo",{"id":list(OUTCOMES)})
meta_rows=rows_from(meta)
for r in meta_rows:
    if r.get("id") in OUTCOMES: r["analysis_label"]=OUTCOMES[r["id"]]
write_tsv(OUT/"STAGE0_OPENGWAS_METADATA.tsv",meta_rows)
print("metadata rows =",len(meta_rows))

print("===== OPENGWAS FILE AVAILABILITY =====")
try:
    files,_=api("/gwasinfo/files",{"id":list(OUTCOMES)})
    file_rows=rows_from(files)
except Exception as e:
    file_rows=[{"status":"warning","error":str(e)}]
write_tsv(OUT/"STAGE0_OPENGWAS_FILES.tsv",file_rows)
print("file rows =",len(file_rows))

print("===== ASSOCIATION SMOKE TEST (NO PROXIES) =====")
report=json.loads(AUDIT.read_text())
exposure=report.get("preferred_exposure")
if not exposure or not Path(exposure).exists(): raise SystemExit("Preferred exposure missing")
rsids=get_rsids(Path(exposure),8)
print("exposure =",exposure); print("test rsids =",",".join(rsids))

assoc_rows=[]; errors=[]
# Keep each request deliberately small; official ieugwasr chunks associations to avoid timeouts.
for oid,label in OUTCOMES.items():
    try:
        obj,h=api("/associations",{
            "variant":rsids[:4], "id":[oid], "proxies":0,
            "r2":0.8, "align_alleles":1, "palindromes":1, "maf_threshold":0.3
        },retries=4)
        rr=rows_from(obj)
        for r in rr:
            r["analysis_label"]=label; r["smoketest_requested_rsids"]=",".join(rsids[:4])
        assoc_rows.extend(rr)
        print(oid,"association rows =",len(rr))
    except Exception as e:
        errors.append({"id":oid,"analysis_label":label,"error":str(e)})
        print("WARNING",oid,"association smoke test failed:",e)

write_tsv(OUT/"STAGE0_ASSOCIATION_SMOKETEST.tsv",assoc_rows if assoc_rows else [{"status":"no_rows_returned"}])
write_tsv(OUT/"STAGE0_ASSOCIATION_SMOKETEST_ERRORS.tsv",errors if errors else [{"status":"none"}])

summary={
 "jwt_ok":True, "metadata_rows":len(meta_rows), "file_rows":len(file_rows),
 "preferred_exposure":exposure, "test_rsids":rsids,
 "association_rows":len(assoc_rows), "association_errors":errors,
 "note":"/tophits intentionally not used; Stage 0 validates exact exposure rsIDs via /associations with proxies=0."
}
(OUT/"STAGE0_OPENGWAS_API_SUMMARY.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding="utf-8")

# Metadata is the hard Stage-0 gate. A transient association 5xx is recorded as warning, not fatal.
if len(meta_rows) < 2: raise SystemExit(2)
print("OpenGWAS Stage 0 API gate = PASS")
