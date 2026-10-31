#!/usr/bin/env python3
"""P5h: discover original IOVS full-gene DEG supplement via Europe PMC XML.

Metadata only. Never invent supplementary URLs or claim replication from names.
"""
import argparse,json,re,time,urllib.request,xml.etree.ElementTree as ET
from datetime import datetime,timezone
from pathlib import Path
SOURCE="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12882097/fullTextXML"
def get(url):
    errors=[]
    for i in range(3):
        try:
            req=urllib.request.Request(url,headers={
                "User-Agent":"ocular-publication-supplement-verification/1.0",
                "Accept":"application/xml"})
            with urllib.request.urlopen(req,timeout=45) as x:
                raw=x.read(12*1024*1024+1)
                if len(raw)>12*1024*1024:raise ValueError("XML_EXCEEDS_12MB_LIMIT")
                if not raw.lstrip().startswith(b"<?xml") and not raw.lstrip().startswith(b"<"):
                    raise ValueError("RESPONSE_IS_NOT_XML")
                return raw
        except Exception as e:
            errors.append(type(e).__name__+":"+str(e)[:100])
            if i<2:time.sleep(i+1)
    raise RuntimeError(";".join(errors))
def flatten(x):
    return " ".join(" ".join(x.itertext()).split())
def tags(root):
    items=[]
    for el in root.iter():
        name=el.tag.rsplit("}",1)[-1]
        if name not in ("supplementary-material","supplementary-materials"):
            continue
        hrefs=[]
        for sub in el.iter():
            for k,v in sub.attrib.items():
                if k.rsplit("}",1)[-1] in ("href","filename"):
                    hrefs.append(v)
        title=" ".join((flatten(y) for y in el if y.tag.rsplit("}",1)[-1] in ("label","caption","title")))
        if not title:title=flatten(el)[:120]
        items.append({"label":title[:250],"references":list(dict.fromkeys(hrefs))[:15]})
    return items
def main():
    a=argparse.ArgumentParser();a.add_argument("--output",required=True);p=a.parse_args()
    out={"source":SOURCE,"doi":"10.1167/iovs.67.2.13",
      "pmcid":"PMC12882097","status":"HOLD","utc":datetime.now(timezone.utc).isoformat()}
    try:
        doc=ET.fromstring(get(SOURCE))
        supplements=tags(doc)
        out["n_supplement_nodes"]=len(supplements)
        out["supplements"]=supplements
        candidates=[v for v in supplements if re.search(
           r"\bsupp(?:lementary)?\s*(?:table|tbl)\s*4\b",v["label"],re.I)]
        out["candidate_full_TM_all_gene_table4"]=candidates
        out["status"]="PASS_SUPPLEMENT_TABLE4_REFERENCE" if candidates and any(x["references"] for x in candidates) else "HOLD_TABLE4_FILE_REFERENCE_MISSING"
        out["next_gate"]="Download only verified referenced Table 4 and parse actual per-gene log2FC/FDR; then compare on matching identifiers"
    except Exception as e:
        out["error"]=type(e).__name__+": "+str(e)[:250]
    dest=Path(p.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print("P5H",out["status"],"supplement_nodes",out.get("n_supplement_nodes"),
          "table4",len(out.get("candidate_full_TM_all_gene_table4",[])),
          "error",out.get("error"))
    if out["status"]!="PASS_SUPPLEMENT_TABLE4_REFERENCE":raise SystemExit(2)
if __name__=="__main__":main()
