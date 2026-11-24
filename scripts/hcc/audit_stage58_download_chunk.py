#!/usr/bin/env python3
"""Stage58: inspect only source-proven HCCDB repository/download route JS.

Paths AND SHA256 below were observed directly in Stage57 clean Chromium
Actions #41438110070. Fetching is strictly bounded and SHA-pinned. Inspect
anonymous code *structure*, not patient data, private APIs, or inferred download
URLs. A public static chunk does not prove downloadable RNA matrices.
"""
import argparse
import datetime as dt
import hashlib
import json
import pathlib
import re
import urllib.parse
import urllib.request

ORIGIN="http://lifeome.net:809"
CHUNKS={
 "repository":{"path":"/static/js/repository.509a9848.js","bytes":108379,
      "sha256":"cf02b87d3424e0f1605f0d7b772948df51f6e9c52820445a2850e13a3b568893"},
 "download":{"path":"/static/js/download.b22d03a9.js","bytes":11335,
      "sha256":"973cddbabef72e012a6be3dd2c418f71c5e4ac27c851bd4d70ebd7a28b42c777"},
}
MAX_JS=512*1024
TARGET_FILES=("HCCDB25_mRNA_level3.txt","HCCDB25.sample.txt","HCCDB25.patient.txt")
KEYWORDS=("download","repository","buke","bulk","HCCDB25","api","axios","fetch(",
          "window.open","location.href","createElement","POST","GET","zip","tsv","rsem")
RAW_LITERAL=re.compile(r"""["']([A-Za-z0-9_./:$?=&{}-]{2,140})["']""")
BLOCK_ID=re.compile(r"(?i)(?:OEP|OED|OES|OEZ)[0-9]{4,}|(?:token|secret|password|bearer|patient)")

def allowed_path(path):
    """Whitelist only two browser-observed first-party static JS paths."""
    if path not in {v["path"] for v in CHUNKS.values()}:
        raise ValueError("JS_PATH_NOT_OBSERVED_BY_PUBLIC_BROWSER")
    return ORIGIN+path

def exact_public_js(chunk,opener=urllib.request.urlopen):
    target=CHUNKS[chunk]
    url=allowed_path(target["path"])
    req=urllib.request.Request(url,headers={
        "User-Agent":"HCC-Stage58-source-pinned-public-JS/1.0",
        "Range":"bytes=0-"+str(MAX_JS-1),
        "Accept":"application/javascript,text/javascript,*/*"})
    with opener(req,timeout=35) as r:
        if r.status not in (200,206) or r.geturl()!=url:
            raise ValueError("SOURCE_ORIGIN_OR_STATUS_UNVERIFIED")
        data=r.read(MAX_JS+1)
        content_type=r.headers.get("Content-Type","").split(";")[0].strip().lower()
        content_range=r.headers.get("Content-Range","")
        length=r.headers.get("Content-Length")
        status=r.status
    if len(data)!=target["bytes"] or len(data)>MAX_JS:
        raise ValueError("SOURCE_JS_EXPECTED_BYTES_MISMATCH")
    if status==206:
        m=re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",content_range.strip())
        complete=bool(m and int(m[1])==0 and int(m[2])+1==len(data)
                      and int(m[3])==len(data))
    else:
        complete=length is not None and int(length)==len(data)
    if not complete:
        raise ValueError("PUBLIC_JS_COMPLETE_SOURCE_NOT_VERIFIED")
    if content_type not in ("application/javascript","text/javascript",
                            "application/x-javascript","text/plain","application/octet-stream"):
        raise ValueError("SOURCE_JS_MIME_UNEXPECTED")
    digest=hashlib.sha256(data).hexdigest()
    if digest!=target["sha256"]:
        raise ValueError("PUBLIC_JS_PINNED_SHA256_MISMATCH")
    return data,{"bytes":len(data),"sha256":digest,"pinned_hash_matches":True,
                 "http_status":status,"full_source":True}

def safe_literal(value):
    """Report only code-owned, privacy-safe path hints, never query params or IDs."""
    if not isinstance(value,str) or not value or len(value)>135:
        return None
    if BLOCK_ID.search(value) or "//" in value or "@" in value:
        return None
    try:
        url=urllib.parse.urlsplit(value)
        if url.scheme or url.netloc: return None
        value=url.path
        if not value.startswith("/"):value="/"+value
        if not re.fullmatch(r"/[A-Za-z0-9_./:$-]{2,125}",value):
            return None
        if any(x in (".","..") for x in value.split("/")):
            return None
        if not any(x in value.lower() for x in ("api","download","repository",
               "bulk","buke","file","data","hccdb")):return None
        return value
    except (ValueError,UnicodeError):
        return None

def inspect_chunk(data,chunk):
    text=data.decode("utf-8","replace")
    count={kw:len(re.findall(re.escape(kw),text,re.I)) for kw in KEYWORDS}
    literal=set()
    for match in RAW_LITERAL.finditer(text):
        v=safe_literal(match.group(1))
        if v:literal.add(v)
    # Only structural proximity counts: no JS snippets, original URLs or
    # potentially identifying metadata can be emitted in logs or artifacts.
    api_near_download=0
    for m in re.finditer("download",text,re.I):
        nearby=text[max(0,m.start()-400):m.end()+400]
        if re.search(r"axios|fetch\(|XMLHttpRequest|\.post\(|\.get\(",nearby,re.I):
            api_near_download+=1
    return {"chunk":chunk,"source_keyword_counts":count,
       "source_literal_safe_path_hints":sorted(literal)[:35],
       "source_literal_safe_path_hints_n":len(literal),
       "download_word_near_request_call_n":api_near_download,
       "expected_RNA_triplet_file_references":{name:text.count(name) for name in TARGET_FILES},
       "verified_RNA_original_matrix_content":False,
       "verified_official_NODE_or_PDC_sample_identity":False}

def audit(chunk,opener=urllib.request.urlopen):
    out={"chunk":chunk,"status":"PUBLIC_ROUTE_CHUNK_SOURCE_HOLD",
         "rna_matrix_downloaded":False,
         "official_RNA_PDC_biopsy_identity":False}
    try:
        raw,source=exact_public_js(chunk,opener)
        out["source"]=source
        out["analysis"]=inspect_chunk(raw,chunk)
        out["status"]="PUBLIC_PINNED_ROUTE_CHUNK_ANALYZED"
    except Exception as exc:
        out["error_class"]=type(exc).__name__
        out["error_code"]=str(exc) if isinstance(exc,ValueError) else "SOURCE_TRANSPORT_HOLD"
    return out

def review(out):
    if out.get("rna_matrix_downloaded") is not False or (
        out.get("official_RNA_PDC_biopsy_identity") is not False):
        raise ValueError("BIOLOGICAL_GO_UNSUPPORTED")
    return {"schema":"hcc-stage58-public-route-source-postflight-v1",
            "official_source_js_verified":out.get("status")=="PUBLIC_PINNED_ROUTE_CHUNK_ANALYZED",
            "RNA_processed_expression_file_verified":False,
            "G30_B":"HOLD_PDC_OFFICIAL_ALIQUOT_ID",
            "G30_C":"HOLD_HCCDB25_RNA_MATRIX_AND_PATIENT_IDS",
            "G30_D":"STOP_NO_SAME_BIOPSY_RNA_PROTEIN",
            "biological_claims_allowed":False}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--chunk",required=True,choices=sorted(CHUNKS))
    ap.add_argument("--out",required=True)
    args=ap.parse_args()
    source=audit(args.chunk)
    report={"schema":"hcc-stage58-source-pinned-download-module-analysis-v1",
            "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "source":source,"review":review(source)}
    path=pathlib.Path(args.out)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print("HCC_STAGE58_ROUTE_CHUNK",json.dumps(report))
if __name__=="__main__":main()
