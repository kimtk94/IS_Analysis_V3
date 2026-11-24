#!/usr/bin/env python3
"""Find verifiable per-run timepoint and already-provided counts in official article/ENA.

No raw FASTQ download. No inferred treatment timepoint from sample suffix.
Result HOLD is a scientific conclusion, independent of CI SUCCESS.
"""
import argparse,urllib.request,urllib.parse,time,xml.etree.ElementTree as ET
import hashlib,json,re,io,csv,zipfile,os
from pathlib import Path
from collections import Counter

PMC="PMC12882097"
XML_URL="https://www.ebi.ac.uk/europepmc/webservices/rest/"+PMC+"/fullTextXML"
ENA_URL="https://www.ebi.ac.uk/ena/browser/api/xml/"
SHARE_MAX=35_000_000
def get(url,limit=3_000_000):
    last=""
    for retry in range(2):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"ocular-master-thesis-PMC-ENA-source-audit/1.0",
              "Accept":"application/xml,text/xml,application/zip,*/*"})
            with urllib.request.urlopen(req,timeout=36) as f:
                body=f.read(limit+1)
                final=f.url;ctype=f.headers.get("Content-Type","")
            if not body or len(body)>limit:raise ValueError("RESPONSE_BAD_SIZE")
            return body,final,ctype
        except (TimeoutError,ValueError,OSError) as e:
            last=type(e).__name__+" "+str(e)[:150]
            if retry==0:time.sleep(2)
    raise RuntimeError(last)
def node_text(node):
    return re.sub(r"\s+"," "," ".join(node.itertext())).strip()
def scan_article():
    b,final,ctype=get(XML_URL,limit=9_000_000)
    root=ET.fromstring(b)
    ns="{http://www.w3.org/1999/xlink}href"
    refs=[]
    for x in root.findall(".//supplementary-material"):
        attr={"id":x.get("id"),"supplementary_type":x.get("content-type"),
            "label":node_text(x)[:160],"hrefs":[]}
        for child in x.iter():
            h=child.get(ns) or child.get("href")
            if h:attr["hrefs"].append(h)
        refs.append(attr)
    links=[]
    for x in root.iter():
        href=x.get(ns) or x.get("href")
        if href and re.search(r"\.(xlsx?|tsv|csv|txt|gz|zip|pdf|docx?)($|\?)",href,re.I):
            links.append({"tag":x.tag.rsplit("}",1)[-1],"href":href,
                          "label":node_text(x)[:100]})
    # Fetch no supplementary data yet: report exact official URLs and filenames.
    body=root.find(".//body")
    if body is None:raise ValueError("OFFICIAL_PAPER_BODY_NOT_FOUND")
    methods=[node_text(x) for x in body.findall(".//sec") if
        re.search(r"(rna sequencing|cell culture|treatment|bioinformatic|patient)",node_text(x)[:220],re.I)]
    checks={
      "journal_fullTextXML_source":XML_URL,"xml_sha256":hashlib.sha256(b).hexdigest(),
      "article_supplementary_records":refs,"direct_data_file_links":links[:60],
      "publication_method_excerpts":[s[:950] for s in methods[:10]],
      "article_text_contains_any_raw_count_link":any(
          re.search(r"count|matrix|featurecounts",x["href"],re.I) for x in links),
      "timepoint_list_from_publication_methods":"1-hour, 6-hour, 2-days",
      "publication_48h_TM_strains_initial":10,"publication_48h_TM_strains_QC_passed":9,
      "publication_48h_SCE_strains_initial":5,"publication_48h_SCE_strains_QC_passed":4,
      "publication_1h_TM_strains":3,"publication_1h_SCE_strains":3,
      "publication_6h_TM_strains":4,"publication_6h_SCE_strains":4,
      "publication_count_not_sample_ID_crosswalk":True}
    return checks

def read_samples(fp):
    with fp.open(encoding="utf8",newline="") as h:
        rows=list(csv.DictReader(h,delimiter="\t"))
    if len(rows)!=46 or len({r["run_accession"] for r in rows})!=46:raise ValueError("SOURCE_46_RUN_DRIFT")
    return rows

def inspect_ena_xml(rows):
    # Read one representative BioSample per nominally paired library, not all
    # 46, to bound API usage. Follow up only when actual time annotations seen.
    pairs={}
    for r in rows:
        alias=r["sample_alias"]
        stem=re.sub(r"(?i)_(Con-EtOH|Dex)_Trypsin(?:_\d+)?$","",alias)
        stem=re.sub(r"(?i)_(vehicle|dex)$","",stem)
        stem=re.sub(r"-[12]$","",stem)
        pairs.setdefault(stem,[]).append(r)
    if len(pairs)!=23 or any(len(v)!=2 for v in pairs.values()):
        raise ValueError("PAIR_GROUP_MISMATCH")
    output=[]
    for stem,rs in sorted(pairs.items()):
        r=rs[0]
        sample=r["sample_accession"]
        url=ENA_URL+sample
        entry={"stem":stem,"representative_biosample":sample,
               "two_runs_in_source_pair":2,
               "timepoint_verified":False,"timepoint_value":"",
               "explicit_time_evidence":[],"other_source_attrs":[]}
        try:
            data,final,ctype=get(url,limit=900_000)
            root=ET.fromstring(data)
            for att in root.findall(".//SAMPLE_ATTRIBUTE"):
                name=att.findtext("TAG","").strip()
                val=att.findtext("VALUE","").strip()
                if re.search(r"time|duration|hour|exposure|treatment",name,re.I):
                    frag={"tag":name[:90],"value":val[:140]}
                    entry["explicit_time_evidence"].append(frag)
                    if re.search(r"time|duration|hour|exposure",name,re.I):
                        m=re.search(r"(?i)\b(1|6|48)\s*(?:h|hour)",val)
                        if m:entry["timepoint_verified"]=True;entry["timepoint_value"]=m.group(1)+"h"
                elif re.search(r"culture|strain|source|tissue",name,re.I):
                    entry["other_source_attrs"].append({"tag":name[:90],"value":val[:140]})
            entry["biosample_xml_sha256"]=hashlib.sha256(data).hexdigest()
            entry["status"]="FETCHED"
        except Exception as e:
            entry["status"]="HOLD_XML_FETCH"
            entry["reason"]=type(e).__name__+": "+str(e)[:160]
        output.append(entry)
    return output

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    info={}
    try:
        info=scan_article()
        info["status"]="PASS_OFFICIAL_PAPER_XML_SUPPLEMENT_INVENTORY"
    except Exception as e:
        info={"status":"HOLD_PAPER_FULLTEXT_UNAVAILABLE",
              "error":type(e).__name__+": "+str(e)[:180]}
    (a.out/"MEHROTRA_PMC_PAPER_SUPPLEMENT_AND_COUNT_FILES.json").write_text(
        json.dumps(info,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    rows=read_samples(a.manifest)
    items=inspect_ena_xml(rows)
    total=sum(r["timepoint_verified"] for r in items)
    enasummary={"status":"GO_23_PAIR_BIOSAMPLE_META_ONLY" if len(items)==23 else "HOLD",
        "source_project":"PRJNA1243143",
        "source_pairs":len(items),
        "pair_representative_biosamples_fetched":sum(r["status"]=="FETCHED" for r in items),
        "pair_timepoint_metadata_verified_count":total,
        "timepoint_distribution_only_for_explicit_source_records":
           dict(Counter(r["timepoint_value"] for r in items if r["timepoint_verified"])),
        "timepoint_by_alias_suffix_inferred":False,
        "48h_TM_analysis_raw_recount_ready":False,
        "raw_FASTQ_files_retrieved":0,
        "biological_strain_count_requires_original_sample_sheet":True}
    (a.out/"MEHROTRA_23_PAIR_BIOSAMPLE_TIMEPOINT_QC.json").write_text(
        json.dumps(enasummary,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    (a.out/"MEHROTRA_23_PAIR_BIOSAMPLE_TIME_EVIDENCE.json").write_text(
        json.dumps(items,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("PAPER",info.get("status"),"supplement",len(info.get("article_supplementary_records",[])),
          "file_links",len(info.get("direct_data_file_links",[])),flush=True)
    print("BIOSAMPLE",json.dumps(enasummary,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
