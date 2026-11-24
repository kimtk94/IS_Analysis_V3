#!/usr/bin/env python3
"""P5c: sample-specific ENA metadata audit for PRJNA1243143, no read downloads."""
import argparse,csv,io,json,re,urllib.parse,urllib.request
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
BASE="https://www.ebi.ac.uk/ena/portal/api/"
WANT=["run_accession","sample_accession","sample_alias","sample_title",
      "sample_description","experiment_alias","experiment_title",
      "library_name","experimental_factor","library_layout","fastq_ftp"]
def fetch(url,cap=4000000):
    request=urllib.request.Request(url,headers={"User-Agent":"ocular-external-replication-audit/1.1"})
    with urllib.request.urlopen(request,timeout=55) as resp:
        b=resp.read(cap+1)
    if len(b)>cap:raise ValueError("TOO_MANY_METADATA_BYTES")
    return b.decode("utf8","replace")
def fields():
    s=fetch(BASE+"returnFields?result=read_run",cap=500000)
    result=set()
    rows=list(csv.DictReader(io.StringIO(s),delimiter="\t"))
    for r in rows:
        result.add(r.get("columnId",""))
    # ReturnFields response differs by API version; fallback to known safe fields.
    return result
def classify(text):
    x=" "+text.lower()+" "
    dex=bool(re.search(r'\b(?:dex|dexamethasone)\b',x))
    ctl=bool(re.search(r'\b(?:vehicle|ethanol|untreated|control)\b',x))
    treatment="AMBIGUOUS" if dex and ctl else "DEX" if dex else "CONTROL" if ctl else "UNKNOWN"
    times=[k for k,p in (("1h",r"\b1\s*(?:h|hr|hour)s?\b"),
                         ("6h",r"\b6\s*(?:h|hr|hour)s?\b"),
                         ("2d",r"\b(?:2\s*day[s]?|48\s*(?:h|hr|hour)s?)\b"))
           if re.search(p,x)]
    time=times[0] if len(times)==1 else "UNKNOWN"
    strains=set(re.findall(r"\b(?:TM|SC)\s?(\d{2,3})\b",text.upper()))
    # distinguish TM from SCE so that strain # alone cannot pair across tissues
    proper=set(re.findall(r"\b(?:TM|SC)\s?\d{2,3}\b",text.upper()))
    return treatment,time,tuple(sorted(x.replace(" ","") for x in proper))
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True);a=p.parse_args()
    result={"study":"PRJNA1243143","run_utc":datetime.now(timezone.utc).isoformat(),
            "status":"HOLD_EXTERNAL_SAMPLE_DESIGN","classifier_policy":"Sample-specific fields only; study-wide DEX titles not treatment labels"}
    try:
        supported=fields()
        requested=[c for c in WANT if c in supported] if {"run_accession","sample_accession"}<=supported else [
            "run_accession","sample_accession","sample_title","sample_alias","experiment_title","library_layout","fastq_ftp"]
        url=BASE+"filereport?"+urllib.parse.urlencode({
          "accession":"PRJNA1243143","result":"read_run",
          "fields":",".join(requested),"format":"tsv"})
        rd=csv.DictReader(io.StringIO(fetch(url)),delimiter="\t")
        rows=list(rd)
        if not rows:raise ValueError("NO_SAMPLE_RECORDS")
        if not {"run_accession","sample_accession"}<=set(rd.fieldnames or []):
            raise ValueError("MISSING_KEY_COLUMNS")
        labels={}
        for r in rows:
            text=" ".join(r.get(col,"") or "" for col in
                ("sample_alias","sample_title","sample_description","experimental_factor","library_name"))
            # Exclude experiment_title because it frequently repeats the parent study DEX title.
            labels[r["sample_accession"]]=classify(text)
        groups=Counter(labels.values())
        treat=Counter(v[0] for v in labels.values())
        time=Counter(v[1] for v in labels.values())
        mapped=defaultdict(set)
        for tr,t,groups in labels.values():
            for strain in groups: mapped[(strain,t)].add(tr)
        confirmed=[(strain,t) for (strain,t),arms in mapped.items()
                   if {"CONTROL","DEX"}<=arms and t!="UNKNOWN"]
        tm=[x for x in confirmed if x[0].startswith("TM")]
        result.update({"runs":len(rows),"unique_biosamples":len(labels),
           "available_field_names":requested,
           "sample_treatment_counts":dict(treat),"sample_time_counts":dict(time),
           "distinct_TM_strains_with_time_matched_arms":len(set(x[0] for x in tm)),
           "TM_strain_timepoints_with_both_arms":len(tm),
           "one_to_one_sample_condition":"NOT_INFERRED_WITHOUT_EXPLICIT_ATTRIBUTES",
           "fastq_links_listed":sum(bool(r.get("fastq_ftp")) for r in rows)})
        if len(tm)>=3:
            result["status"]="CONDITIONAL_GO_EXTERNAL_STRAIN_TIME_METADATA"
        else:
            result["warning"]="Sample-level treatment/time/strain pairing unverified, cannot label 46 run accessions as 46 donors"
    except Exception as e:
        result["error"]=type(e).__name__+": "+str(e)[:230]
    dst=Path(a.out);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("P5C",result["status"],"runs",result.get("runs"),
          "tm_paired",result.get("TM_strain_timepoints_with_both_arms"),
          "warning",result.get("warning"),"error",result.get("error"))
if __name__=="__main__":main()
