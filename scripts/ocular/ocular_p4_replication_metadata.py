#!/usr/bin/env python3
"""Independent replication accession PRJNA1243143: metadata-only treatment mapping."""
import argparse,csv,io,json,re,urllib.parse,urllib.request
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

FIELDS="run_accession,sample_accession,experiment_title,sample_title,library_layout,fastq_ftp"
def get_registry():
    q={"accession":"PRJNA1243143","result":"read_run","fields":FIELDS,"format":"tsv"}
    url="https://www.ebi.ac.uk/ena/portal/api/filereport?"+urllib.parse.urlencode(q)
    req=urllib.request.Request(url,headers={"User-Agent":"ocular-replication-metadata/1.0"})
    with urllib.request.urlopen(req,timeout=55) as r:
        text=r.read(4*1024*1024+1)
    if len(text)>4*1024*1024:raise ValueError("REGISTRY_TOO_LARGE")
    rd=csv.DictReader(io.StringIO(text.decode("utf8","replace")),delimiter="\t")
    if not {"run_accession","sample_accession"}.issubset(rd.fieldnames or []):
        raise ValueError("UNEXPECTED_ENA_COLUMNS "+str(rd.fieldnames))
    return list(rd),rd.fieldnames
def classify(s):
    x=s.lower()
    is_dex=bool(re.search(r"\b(dex(?:amethasone)?|glucocorticoid)\b",x))
    is_ctrl=bool(re.search(r"\b(vehicle|control|untreated|ethanol)\b",x))
    treat=("AMBIGUOUS" if is_dex and is_ctrl else
           "DEX" if is_dex else "CONTROL" if is_ctrl else "UNKNOWN")
    times=[]
    for key,pattern in [
      ("1h",r"\b1\s*(?:h|hr|hour)s?\b"),
      ("6h",r"\b6\s*(?:h|hr|hour)s?\b"),
      ("48h",r"\b(?:48\s*(?:h|hr|hours?)|2\s*days?)\b")]:
        if re.search(pattern,x):times.append(key)
    return treat,(times[0] if len(times)==1 else "UNKNOWN")
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--output",required=True);args=ap.parse_args()
    out={"accession":"PRJNA1243143","run_utc":datetime.now(timezone.utc).isoformat(),
         "status":"HOLD_REPLICATION_DESIGN"}
    try:
        rows,fields=get_registry()
        if not rows:raise ValueError("NO_RUNS_RETURNED")
        assignments=Counter()
        for row in rows:
            title=" ".join(row.get(k,"") for k in ("experiment_title","sample_title"))
            tr,t=classify(title)
            assignments[(tr,t)]+=1
        samples={r["sample_accession"] for r in rows}
        out.update({"run_count":len(rows),"sample_count":len(samples),
            "fields":fields,"layout_counts":dict(Counter(r.get("library_layout") for r in rows)),
            "has_fastq":sum(bool(r.get("fastq_ftp")) for r in rows),
            "treatment_time_counts":{"/".join(k):n for k,n in assignments.items()},
            "treatment_classified":sum(n for (t,_),n in assignments.items() if t in ("DEX","CONTROL")),
            "timepoint_classified":sum(n for (_,t),n in assignments.items() if t!="UNKNOWN"),
            "donor_pairing":"NOT_VERIFIED_FROM_ENA_FILE_REPORT"})
        if out["treatment_classified"]>=0.8*len(rows) and out["timepoint_classified"]>=0.8*len(rows):
            out["status"]="CONDITIONAL_METADATA_GO_ONLY"
    except Exception as e:
        out["error"]=type(e).__name__+": "+str(e)[:180]
    Path(args.output).write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print("REPLICATION_METADATA",out["status"],"runs",out.get("run_count"),
          "treatment_classified",out.get("treatment_classified"),
          "timepoint_classified",out.get("timepoint_classified"),
          "error",out.get("error"))
    # Registry exploration may validly yield HOLD. A successful workflow is not a GO claim.
if __name__=="__main__":main()
