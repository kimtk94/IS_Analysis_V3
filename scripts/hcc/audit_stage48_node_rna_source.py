#!/usr/bin/env python3
"""Stage48: bounded, unauthenticated NODE OEP000321 source discovery.
Do not assume public project HTML proves accessible OED/OEZ records, processed
RNA matrices, or a patient/protein join. No raw identifiers, filenames, links,
tokens, patient metadata, or expression measurements are emitted.
"""
import argparse
import csv
import datetime as dt
import hashlib
from html.parser import HTMLParser
import json
import pathlib
import re
import urllib.error
import urllib.parse
import urllib.request

ACCESSION="OEP000321"
PAGE="https://www.biosino.org/node/project/detail/"+ACCESSION
HELP="https://www.biosino.org/node/help"
LIMIT=262144
OES=re.compile(r"^OES[0-9]{5,}$",re.I)
OED=re.compile(r"^OED[0-9]{5,}$",re.I)
OEZ=re.compile(r"^OEZ[0-9]{5,}$",re.I)
PUBLIC_DOWNLOAD=re.compile(r"^/download/node/data/public/(OED[0-9]+)$",re.I)
PROCESSED=re.compile(r"(?i)(rsem|expression|gene.?counts|fpkm|tpm|gene.?matrix)")

class Anchors(HTMLParser):
    def __init__(self):
        super().__init__();self.links=[];self.scripts=0
    def handle_starttag(self,tag,attrs):
        d=dict(attrs)
        if tag=="a" and d.get("href"):self.links.append(d["href"])
        elif tag=="script" and d.get("src"):self.scripts+=1

def public_href(href):
    """Accept only an explicit public NODE OED download URI, never construct one."""
    try:
        u=urllib.parse.urlsplit(href)
        if u.scheme and (u.scheme!="https" or u.hostname!="www.biosino.org"):return False
        if not u.scheme and href.startswith("//"):return False
        return bool(PUBLIC_DOWNLOAD.fullmatch(u.path) and not u.query and not u.fragment)
    except ValueError:
        return False

def inspect_html(html):
    parser=Anchors();parser.feed(html)
    ids=set(re.findall(r"\b(?:OED|OEZ|OES)[0-9]{5,}\b",html,re.I))
    counted={prefix:sum(s.upper().startswith(prefix) for s in ids)
             for prefix in ("OES","OED","OEZ")}
    links=[h for h in parser.links if public_href(h)]
    return {"landing_embeds_project_accession":ACCESSION in html,
            "html_anchor_links_n":len(parser.links),
            "javascript_src_n":parser.scripts,
            "embedded_accession_token_counts":counted,
            "literal_public_download_anchor_n":len(links),
            "linked_file_ids_proven_to_belong_to_project":False,
            "matrix_file_or_download_verified":False,
            "raw_identifiers_or_links_exported":False}

def bounded_page(url,opener=urllib.request.urlopen):
    if url not in (PAGE,HELP):raise ValueError("non-allowlisted public route")
    req=urllib.request.Request(url,headers={
        "User-Agent":"HCC-OEP000321-public-metadata-audit/1.0",
        "Accept":"text/html","Cache-Control":"no-cache"})
    with opener(req,timeout=20) as resp:
        actual=urllib.parse.urlsplit(resp.geturl())
        if actual.scheme!="https" or actual.hostname!="www.biosino.org":
            raise ValueError("unexpected redirect host")
        mime=resp.headers.get("Content-Type","").split(";")[0].lower()
        if mime not in ("text/html","application/xhtml+xml"):
            raise ValueError("unexpected public route content type")
        payload=resp.read(LIMIT+1)
        if len(payload)>LIMIT:raise ValueError("public HTML inspection limit")
        status=resp.status
    return {"http_status":status,"content_type":mime,"html_sha256":hashlib.sha256(payload).hexdigest(),
            "html_bytes":len(payload)},payload.decode("utf-8","replace")

def key(x):
    return re.sub(r"[^a-z0-9]","",str(x).lower())

def inspect_manifest_rows(rows):
    """Optional NODE official export; do not equate an OEZ/OED record with a matrix."""
    allowed={"projectid","sampleid","dataid","analysisid","dataname","datatype",
             "filename","filetype","datasize","downloadlink"}
    records=list(rows)
    if len(records)>100000:raise ValueError("manifest row cap")
    if not records:return {"status":"EMPTY_SOURCE_MANIFEST","processed_candidate_rows_n":0}
    parsed=[{key(k):str(v or "").strip() for k,v in row.items()} for row in records]
    if not any("projectid" in r for r in parsed):
        return {"status":"PROJECT_COLUMN_MISSING_HOLD","processed_candidate_rows_n":0}
    for r in parsed:
        if any(k not in allowed for k in r):continue
        if r.get("projectid") and r["projectid"]!=ACCESSION:
            raise ValueError("unrelated cohort in source manifest")
    matched=[r for r in parsed if r.get("projectid")==ACCESSION]
    sample={r["sampleid"] for r in matched if OES.fullmatch(r.get("sampleid",""))}
    data={r["dataid"] for r in matched if OED.fullmatch(r.get("dataid",""))}
    analysis={r["analysisid"] for r in matched if OEZ.fullmatch(r.get("analysisid",""))}
    candidates=[r for r in matched if PROCESSED.search(
        " ".join(r.get(k,"") for k in ("dataname","filename","datatype")))]
    explicitly_public=sum(public_href(r.get("downloadlink","")) for r in candidates)
    return {"status":"OFFICIAL_EXPORT_METADATA_CANDIDATES_ONLY",
            "node_source_rows_n":len(matched),
            "node_valid_unique_OES_sample_ids_n":len(sample),
            "node_valid_unique_OED_data_ids_n":len(data),
            "node_valid_unique_OEZ_analysis_ids_n":len(analysis),
            "processed_candidate_rows_n":len(candidates),
            "processed_candidate_rows_with_explicit_public_OED_href_n":explicitly_public,
            "public_download_file_content_inspected":False,
            "processed_RNA_matrix_dimensions_verified":False,
            "PDC_same_biopsy_sample_mapping_verified":False,
            "raw_identifiers_links_and_filenames_exported":False}

def audit(manifest_path=None):
    result={"accession":ACCESSION,"project_url":PAGE,"status":"PUBLIC_SOURCE_UNVERIFIED",
            "raw_file_or_expression_matrix_downloaded":False,
            "NODE_RNA_sample_identity_verified":False,
            "PDC_aliquot_to_NODE_RNA_sample_join_verified":False}
    for label,url in (("project",PAGE),("help",HELP)):
        try:
            transport,html=bounded_page(url)
            entry={"status":"PUBLIC_HTML_READ",**transport}
            if label=="project":entry["html_source_indicators"]=inspect_html(html)
            result[label]=entry
        except urllib.error.HTTPError as exc:
            result[label]={"status":"HTTP_ERROR","http_status":exc.code}
        except Exception as exc:
            result[label]={"status":"ACCESS_HOLD","error_class":type(exc).__name__}
    if manifest_path:
        path=pathlib.Path(manifest_path)
        if path.stat().st_size>10*1024*1024:raise ValueError("manifest bytes cap")
        with path.open("r",encoding="utf-8-sig",newline="") as fh:
            delimiter="\t" if path.suffix.lower()==".tsv" else ","
            result["source_export_manifest"]=inspect_manifest_rows(csv.DictReader(fh,delimiter=delimiter))
    else:
        result["source_export_manifest"]={"status":"OFFICIAL_NODE_EXPORT_NOT_ACQUIRED"}
    result["status"]=("PUBLIC_PROJECT_HTML_ONLY" if
        result.get("project",{}).get("status")=="PUBLIC_HTML_READ" else "SOURCE_ACCESS_HOLD")
    return result

def review(source):
    if source.get("raw_file_or_expression_matrix_downloaded") is not False or (
        source.get("PDC_aliquot_to_NODE_RNA_sample_join_verified") is not False):
        raise ValueError("unsupported biological data access assertion")
    return {"schema":"hcc-stage48-node-original-rna-postflight-v1",
            "official_project_public_page_reachable":source.get("project",{}).get("http_status")==200,
            "NODE_processed_RNA_expression_matrix_verified":False,
            "NODE_sample_to_PDC_official_aliquot_bridge_verified":False,
            "G30_B":"HOLD_PDC_QUANT_COLUMN_OFFICIAL_ALIQUOT_ID",
            "G30_C":"HOLD_NODE_OEP000321_ORIGINAL_RNA_MATRIX_AND_IDS",
            "G30_D":"STOP_NO_INDEPENDENT_MATCHED_PATIENT_RNA_PROTEIN",
            "biological_or_clinical_claims_allowed":False}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    p.add_argument("--manifest-path",default=None)
    args=p.parse_args()
    obj={"schema":"hcc-stage48-node-original-source-inventory-v1",
         "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
         "source":audit(args.manifest_path)}
    rev=review(obj["source"])
    for dest,content in ((args.out,obj),(args.review_out,rev)):
        path=pathlib.Path(dest);path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(content,indent=2)+"\n",encoding="utf-8")
    print("HCC_STAGE48_NODE_SOURCE",json.dumps(obj))
    print("HCC_STAGE48_REVIEW",json.dumps(rev))
if __name__=="__main__":main()
