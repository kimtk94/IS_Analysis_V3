#!/usr/bin/env python3
"""Bounded official PMC Cloud version and supplementary-table inventory.
No paywall access, no private data, no automatic candidate gene reconstruction.
"""
import argparse,hashlib,io,json,re,time,zipfile
from pathlib import Path
from urllib.parse import urlparse,urlencode,parse_qs,quote,unquote
from urllib.request import Request,urlopen
from xml.etree import ElementTree as ET
import openpyxl

BASE="https://pmc-oa-opendata.s3.amazonaws.com/"
HOST="pmc-oa-opendata.s3.amazonaws.com"
PAPERS={"PNAS":{"pmcid":"PMC13778663","pmid":"43013538","doi":"10.1073/pnas.2609070123"},"IOVS_ERRATUM":{"pmcid":"PMC13712379","doi":None}}
def req(ok,msg):
    if not ok:raise ValueError(msg)
def fetch(url,cap):
    req(urlparse(url).scheme=="https" and urlparse(url).hostname==HOST,"UNAPPROVED_CLOUD_HOST")
    with urlopen(Request(url,headers={"Accept-Encoding":"identity","User-Agent":"Ocular-research-PMC-source-inventory/1.0"}),timeout=40) as r:
        req(r.status==200 and urlparse(r.geturl()).hostname==HOST,"HTTP_OR_HOST")
        n=r.headers.get("Content-Length")
        if n:req(0<int(n)<=cap,"SOURCE_SIZE_CAP")
        raw=r.read(cap+1)
    req(0<len(raw)<=cap and (not n or len(raw)==int(n)),"SOURCE_SIZE_OR_TRUNCATION")
    return raw,{"url":url,"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest()}
def cloud_url(value,prefix):
    u=urlparse(value)
    if u.scheme=="s3" and u.netloc=="pmc-oa-opendata":
        key=u.path.lstrip("/")
    elif u.scheme=="https" and u.hostname==HOST:
        key=u.path.lstrip("/")
    else:raise ValueError("UNAPPROVED_METADATA_OBJECT_URL")
    decoded=unquote(key)
    req(decoded.startswith(prefix) and ".." not in decoded.split("/") and "/" not in decoded[len(prefix):],"WRONG_ARTICLE_OBJECT")
    md5=parse_qs(u.query).get("md5",[])
    req(len(md5)==1 and re.fullmatch("[a-fA-F0-9]{32}",md5[0]),"NO_PUBLISHED_MD5")
    return BASE+quote(decoded,safe="/"),decoded,md5[0].lower()
def object_read(value,prefix,cap):
    url,key,md5=cloud_url(value,prefix)
    raw,src=fetch(url,cap);req(hashlib.md5(raw).hexdigest()==md5,"PUBLISHED_MD5_MISMATCH")
    return raw,{**src,"key":key,"published_md5":md5}
def inventory(p):
    ns={"s":"http://s3.amazonaws.com/doc/2006-03-01/"}
    raw,source=fetch(BASE+"?"+urlencode({"list-type":"2","prefix":p["pmcid"]+".","delimiter":"/","max-keys":20}),200000)
    xml=ET.fromstring(raw)
    req(xml.findtext("s:IsTruncated",namespaces=ns)=="false","VERSION_LIST_INCOMPLETE")
    prefixes=[n.text for n in xml.findall("s:CommonPrefixes/s:Prefix",ns)]
    req(len(prefixes)<=10 and all(re.fullmatch(p["pmcid"]+r"\.[1-9][0-9]*/",v) for v in prefixes),"VERSION_PREFIX_IDENTITY")
    out={"status":"VERSION_INVENTORY_READ","version_listing":source,"versions":prefixes,"articles":[]}
    if not prefixes:return {**out,"status":"HOLD_NO_PUBLIC_CLOUD_VERSION"}
    for prefix in sorted(prefixes,key=lambda s:int(s.split(".")[-1][:-1])):
        article={"prefix":prefix,"status":"HOLD"}
        try:
            ver=prefix[:-1]
            raw,meta_source=fetch(BASE+"metadata/"+ver+".json",500000)
            m=json.loads(raw)
            req(m["pmcid"]==p["pmcid"] and str(m["version"])==ver.split(".")[-1] and m.get("is_retracted") is False,"ARTICLE_IDENTITY_OR_RETRACTION")
            req(p.get("doi") is None or str(m.get("doi")).lower()==p["doi"],"DOI_MISMATCH")
            req(p.get("pmid") is None or str(m.get("pmid"))==p["pmid"],"PMID_MISMATCH")
            req(p.get("doi") is not None or ("erratum" in str(m.get("title","")).lower() or "correction" in str(m.get("title","")).lower()) and ("trabecular" in str(m.get("title","")).lower() or "glucocorticoid" in str(m.get("title","")).lower()),"ERRATUM_TITLE_IDENTITY")
            article.update(metadata_source=meta_source,metadata={k:m.get(k) for k in ("pmcid","version","pmid","doi","title","citation","is_pmc_openaccess","is_manuscript","is_retracted","license_code")})
            # Access only objects explicitly exposed by the official PMC metadata.
            if m.get("is_pmc_openaccess") is not True:
                article["status"]="HOLD_NOT_CONFIRMED_OPEN_ACCESS";out["articles"].append(article);continue
            req(isinstance(m.get("license_code"),str) and bool(m["license_code"]),"LICENSE_UNVERIFIED")
            raw,src=object_read(m["xml_url"],prefix,5000000)
            root=ET.fromstring(raw)
            ids={a.get("pub-id-type"):a.text for a in root.findall(".//article-meta/article-id")}
            if p.get("doi"):req(ids.get("doi")==p["doi"],"XML_DOI_MISMATCH")
            article["xml_source"]=src
            refs=[]
            for tag in root.iter():
                if tag.tag.split("}")[-1] not in ("supplementary-material","media"):continue
                href=tag.get("{http://www.w3.org/1999/xlink}href")
                if href:
                    label=" ".join("".join(tag.find("label").itertext()).split()) if tag.find("label") is not None else ""
                    # Retain filenames/labels only, not full publication prose.
                    refs.append({"href":href[:250],"label":label[:100]})
            article["supplement_references"]=refs[:150]
            # Extract accessions only; do not reproduce full paper.
            article["public_accession_tokens"]=sorted(set(re.findall(r"\b(?:GSE[0-9]+|PRJNA[0-9]+|SRP[0-9]+)\b"," ".join(root.itertext()))))
            media=m.get("media_urls") or []
            req(isinstance(media,list) and len(media)<=500 and all(isinstance(v,str) for v in media),"MEDIA_SCHEMA")
            listed=[]
            for v in media:
                url,key,md5=cloud_url(v,prefix)
                listed.append({"key":key,"url":url,"md5":md5})
            article["supplement_file_inventory"]=[v for v in listed if re.search(r"\.(xlsx|xls|tsv|csv|txt|zip|pdf)$",v["key"],re.I)]
            referenced={unquote(v["href"].rsplit("/",1)[-1]) for v in refs}
            excel=[v for v in media if unquote(urlparse(v).path).lower().endswith(".xlsx") and unquote(urlparse(v).path.rsplit("/",1)[-1]) in referenced]
            article["referenced_xlsx_count"]=len(excel)
            article["workbooks"]=[]
            if len(excel)>8:
                article["status"]="HOLD_REVIEW_XLSX_INVENTORY";out["articles"].append(article);continue
            for u in excel:
                try:
                    time.sleep(1)
                    data,ms=object_read(u,prefix,20000000)
                    with zipfile.ZipFile(io.BytesIO(data)) as z:
                        req(len(z.infolist())<=5000 and sum(x.file_size for x in z.infolist())<=120000000,"XLSX_EXPANSION_CAP")
                    book=openpyxl.load_workbook(io.BytesIO(data),read_only=True,data_only=True)
                    try:
                        req(len(book.sheetnames)<=80,"SHEET_COUNT_CAP")
                        tabs=[]
                        for sheet in book:
                            # Schema only: first 6 rows and up to 24 cols, no full gene list.
                            first=[]
                            for ix,row in enumerate(sheet.iter_rows(min_row=1,max_row=6,max_col=24,values_only=True),1):
                                # Retain only first two non-empty rows, capped at 160 chars each cell.
                                if any(v is not None for v in row):
                                    first.append({"row":ix,"values":[str(v)[:160] if v is not None else None for v in row]})
                                if len(first)==2:break
                            tabs.append({"name":sheet.title,"reported_max_row":sheet.max_row,"reported_max_column":sheet.max_column,"schema_preview":first})
                        article["workbooks"].append({"status":"REFERENCED_WORKBOOK_SCHEMA_READ","source":ms,"sheets":tabs})
                    finally:book.close()
                except (OSError,ValueError,KeyError,TypeError,zipfile.BadZipFile) as ex:
                    article["workbooks"].append({"status":"HOLD_WORKBOOK","source":u,"reason":str(ex)[:200]})
            article["status"]="PUBLIC_SOURCE_INVENTORY_COMPLETE"
        except (OSError,ValueError,TypeError,KeyError,ET.ParseError) as ex:article["reason"]=str(ex)[:250]
        out["articles"].append(article)
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    out={"status":"PMC_OFFICIAL_PUBLIC_SOURCE_INVENTORY_NOT_OVERLAP","source":"https://registry.opendata.aws/ncbi-pmc","documentation":"https://pmc.ncbi.nlm.nih.gov/tools/pmcaws/","papers":{},"limits":["Official public cloud only. No full-text web scraping, paywall bypass, restricted records, or redistribution of complete publisher supplementary files.","Verify article version, DOI/PMID and published object MD5. Workbook schema previews are not full gene sets. Do not reconstruct 67 genes from counts or select sheets by matching expected count.","No panel overlap or novelty result until exact published table identities and complete rows are verified."]}
    for name,p in PAPERS.items():
        try:out["papers"][name]=inventory(p)
        except (OSError,ValueError,TypeError,KeyError,ET.ParseError) as ex:out["papers"][name]={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(ex)[:250]}
        time.sleep(1)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n");print(json.dumps(out,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
