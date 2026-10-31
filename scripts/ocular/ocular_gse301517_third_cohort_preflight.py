#!/usr/bin/env python3
"""GSE301517 official GEO metadata + normalized CPM P0 access preflight.

Fetch bounded, checksum the public Series SOFT and CPM table; no raw FASTQ,
individual-level human genotypes, or donor profiles are published by this audit.
GEO sample replicates must never be labeled 14 independent donor samples.
"""
import argparse,csv,gzip,hashlib,io,json,re,time,urllib.request
from collections import Counter,defaultdict
from pathlib import Path

GSE="GSE301517"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE301nnn/GSE301517"
CPM=BASE+"/suppl/GSE301517_Control_vs_Dex_Gene_CPM_MSTR.tsv.gz"
SOFT=BASE+"/soft/GSE301517_family.soft.gz"
EXPECTED={
 "TM62":{"Control":4,"Dex":4},
 "TM70":{"Control":2,"Dex":2},
 "TM72":{"Control":1,"Dex":1},
}
REFERENCE={"TM62":3,"TM71":2,"TM83":2}
CORE=["MYOC","TGFBI","MMP2","ACTA2","SDC1","EGR3","LPAR1","TFPI"]
def fetch(url, cap):
    error=None
    for attempt in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"ocular-GSE301517-P0-manifest/1.0"})
            with urllib.request.urlopen(req,timeout=110) as res:
                data=res.read(cap+1)
            if len(data)>cap:raise ValueError("REMOTE_FILE_EXCEEDS_SIZE_CAP")
            if data[:2]!=b"\x1f\x8b":raise ValueError("NON_GZIP_SOURCE")
            return data
        except (TimeoutError, OSError, ValueError) as exc:
            error=type(exc).__name__+": "+str(exc)[:180]
            if attempt<2:time.sleep(2*(attempt+1))
    raise RuntimeError("GEO_DOWNLOAD_FAILED "+str(error))
def read_soft(raw):
    result={}
    current=None
    with gzip.open(io.BytesIO(raw),"rt",encoding="utf8",errors="replace") as f:
        for line in f:
            if line.startswith("^SAMPLE = "):
                current=line.split("=",1)[1].strip()
                result[current]={"gsm":current,"title":""}
            if current is None:continue
            if line.startswith("!Sample_title = "):
                result[current]["title"]=line.split("=",1)[1].strip()
            elif line.startswith("!Sample_treatment_protocol_ch1 = "):
                result[current].setdefault("treatment_protocol",[]).append(line.split("=",1)[1].strip()[:400])
            elif line.startswith("!Sample_growth_protocol_ch1 = "):
                result[current].setdefault("growth_protocol",[]).append(line.split("=",1)[1].strip()[:250])
            elif line.startswith("!Sample_characteristics_ch1 = "):
                result[current].setdefault("characteristics",[]).append(line.split("=",1)[1].strip()[:250])
    return result
def classify(title):
    t=title.strip()
    arm="Resource" if re.search(r"\bResource\b",t,re.I) else ("Control" if re.search(r"\bControl\b",t,re.I) else "Dex" if re.search(r"\bDex\b",t,re.I) else "UNKNOWN")
    m=re.search(r"\bTM\s*([0-9]+)\b",t,re.I)
    strain="TM"+m.group(1) if m else "UNKNOWN"
    rep=re.search(r"\bRep\s*,?\s*(\d+)\b",t,re.I)
    return arm,strain,int(rep.group(1)) if rep else None
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    report={"status":"HOLD_UNVERIFIED","dataset":GSE,
      "normalization":"CPM (published table); NOT RAW INTEGER COUNTS",
      "independent_biological_strain_count_expected":3,
      "technical_replicates_are_not_donors":True,
      "source_urls":{"cpm":CPM,"soft":SOFT},
      "clinical_steroid_delta_iop":"NOT_AVAILABLE"}
    raw=fetch(SOFT,30_000_000)
    entries=read_soft(raw)
    expected_all=[r for r in entries.values() if r["gsm"].startswith("GSM")]
    report["soft_sha256"]=hashlib.sha256(raw).hexdigest()
    report["GSM_count"]=len(expected_all)
    counts=Counter()
    per_strain_reps=defaultdict(set)
    unexpected=[]
    meta=[]
    for row in expected_all:
        arm,strain,rep=classify(row["title"])
        counts[(strain,arm)]+=1
        if rep is not None:per_strain_reps[(strain,arm)].add(rep)
        if arm=="UNKNOWN" or strain=="UNKNOWN":unexpected.append(row["gsm"])
        # Published GSM identifiers only, no inferred human donor IDs.
        meta.append({"gsm":row["gsm"],"arm":arm,"strain":strain,"replicate_id":rep,
          "is_treatment_comparison":arm in ("Dex","Control")})
    checks=all(counts[(s,a)]==n for s,d in EXPECTED.items() for a,n in d.items())
    checks=checks and all(counts[(s,"Resource")]==n for s,n in REFERENCE.items())
    checks=checks and len(expected_all)==21 and not unexpected
    checks=checks and all(len(per_strain_reps[(s,a)])==n for s,d in EXPECTED.items() for a,n in d.items())
    report["group_counts"]={s:{a:counts[(s,a)] for a in ("Control","Dex")} for s in EXPECTED}
    report["resource_group_counts"]={s:counts[(s,"Resource")] for s in REFERENCE}
    report["sample_design_PASS"]=checks
    report["GSM_accessions_expected_only"]=sorted(x["gsm"] for x in expected_all)
    report["treatment_protocol_text_unique"]=sorted(set(v for r in expected_all for v in r.get("treatment_protocol",[])))[:12]
    if not checks:report["status"]="HOLD_GEO_SAMPLE_DESIGN_DRIFT"
    with (args.out/"GSE301517_21_published_sample_manifest.tsv").open("w",encoding="utf8",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=list(meta[0]),delimiter="\t")
        wr.writeheader();wr.writerows(sorted(meta,key=lambda x:x["gsm"]))
    cpm=fetch(CPM,26_000_000)
    report["cpm_sha256"]=hashlib.sha256(cpm).hexdigest()
    report["cpm_compressed_bytes"]=len(cpm)
    # Complete gzip/CRC through text-mode streaming. Inspect every row and header,
    # but keep only aggregate schema; do not publish raw sample-expression vectors.
    with gzip.open(io.BytesIO(cpm),"rt",encoding="utf8",newline="") as f:
        header=f.readline().rstrip("\r\n").split("\t")
        if len(header)<4:raise ValueError("GEO_CPM_TOO_FEW_COLUMNS")
        examples=[];rows=0;wrong_columns=0
        first_col=set()
        for line in f:
            rows+=1
            fields=line.rstrip("\r\n").split("\t")
            if len(fields)!=len(header):wrong_columns+=1
            if len(examples)<4:examples.append({"identifier_preview":fields[:min(3,len(fields))],
              "field_count":len(fields)})
            if fields: first_col.add(fields[0])
            if rows>150000:raise ValueError("CPM_ROW_LIMIT_REACHED")
    report["cpm_complete_gzip_crc"]="PASS"
    report["cpm_rows"]=rows
    report["cpm_columns"]=len(header)
    report["cpm_header"]=header[:90]
    report["cpm_initial_rows_preview"]=examples
    report["first_column_unique_values"]=len(first_col)
    report["row_width_errors"]=wrong_columns
    report["first_column_unique"]=len(first_col)==rows
    report["header_contains_GSM_count"]=sum(any(g==h for g in report["GSM_accessions_expected_only"]) for h in header)
    report["status"]="GO_SOURCE_DESIGN_CPM_FORMAT_PREFLIGHT" if checks and rows>10000 and not wrong_columns and len(first_col)==rows else "HOLD_CPM_SCHEMA_OR_DESIGN"
    (args.out/"GSE301517_P0_official_access_and_schema.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("GSE301517_P0",json.dumps({k:report[k] for k in ("status","GSM_count","group_counts","resource_group_counts","cpm_rows","cpm_columns","cpm_sha256","cpm_header","cpm_initial_rows_preview")},ensure_ascii=False),flush=True)
if __name__=="__main__":main()
