#!/usr/bin/env python3
"""Stage36 read-only PDC Gao GCT heatmap path discovery and GCT 1.3 header gate.

Original GCT is not an ordinary PDC Protein Assembly file. Portal page visibility
or a study heatmap link does not prove same 330 quantDataMatrix columns. This
module does not invent GCT URLs. Only exact source aliases can form candidates.
No source URLs, individual IDs, per-column hashes or matrix values are emitted.
"""
import argparse
import collections
import datetime as dt
import hashlib
from html.parser import HTMLParser
import io
import json
import pathlib
import re
import urllib.error
import urllib.parse
import urllib.request

from probe_stage30_gao_catalog import PDC
from probe_stage30_gao_quant_api_header_to_biospecimen import score

MAX_HTML_BYTES=8*1024*1024
MAX_GCT_HEADER_BYTES=512*1024
PORTAL="https://pdc.cancer.gov/pdc/explore-quantitation-data"
DOCUMENTATION="https://pdc-docs.cancer.gov/pdc/explore-quantitation-data"
STUDY_PREFIX="https://pdc.cancer.gov/pdc/study/"
SAFE_PAGE_HOSTS={"pdc.cancer.gov","pdc-docs.cancer.gov"}

class LinkFinder(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.links=[]
        self.rows=[]
        self.current=None
        self.current_row=None
    def handle_starttag(self,tag,attrs):
        if tag=="tr":
            self.current_row={"txt":[],"links":[]}
        if tag=="a":
            self.current={"href":dict(attrs).get("href",""),"txt":[]}
    def handle_data(self,data):
        if self.current is not None:
            self.current["txt"].append(data)
        if self.current_row is not None:
            self.current_row["txt"].append(data)
    def handle_endtag(self,tag):
        if tag=="a" and self.current is not None:
            link=(self.current["href"],"".join(self.current["txt"]).strip())
            self.links.append(link)
            if self.current_row is not None:
                self.current_row["links"].append(link)
            self.current=None
        if tag=="tr" and self.current_row is not None:
            self.rows.append(self.current_row)
            self.current_row=None

def is_safe_page(url):
    try:
        p=urllib.parse.urlparse(url)
        return p.scheme=="https" and p.hostname in SAFE_PAGE_HOSTS and not p.username and not p.password
    except ValueError:
        return False

def classify_anchor(href,text,base):
    """No URL returned; classify canonical path as explicit GCT or heatmap."""
    if not href or href.lstrip().lower().startswith(("javascript:","data:","blob:","#")):
        return None
    url=urllib.parse.urljoin(base,href)
    p=urllib.parse.urlparse(url)
    if p.scheme!="https" or not p.hostname:
        return None
    path=urllib.parse.unquote(p.path).lower()
    explicit_gct=bool(re.search(r"\.gct(?:\.gz)?$",path))
    heatmap=any(tok in (p.path+" "+text).lower() for tok in ("heatmap","morpheus","log2 ratio"))
    if not (explicit_gct or heatmap):
        return None
    if explicit_gct:
        return "explicit_gct"
    return "heatmap_or_viewer"

def page_inventory(raw,base,pid):
    txt=raw.decode("utf-8","replace")
    parser=LinkFinder()
    parser.feed(txt)
    # A link is study-associated ONLY if the anchor text, relative href, or
    # its enclosing local document window contains the explicit accession.
    direct_matching=[(url,label) for url,label in parser.links if pid in url or pid in label]
    row_matching=[link for row in parser.rows if pid in "".join(row["txt"]) for link in row["links"]]
    matching=list(set(direct_matching+row_matching))
    counts=collections.Counter(classify_anchor(href,label,base) for href,label in matching)
    # Requiring accession in the anchor prevents unrelated study links from
    # being misrepresented as an accessible source for the requested study.
    return {
        "page_mentions_target_study":pid in txt,
        "target_study_associated_links_n":len(matching),
        "target_study_explicit_GCT_links_n":counts["explicit_gct"],
        "target_study_heatmap_or_viewer_links_n":counts["heatmap_or_viewer"],
        "all_page_explicit_GCT_links_unassigned_n":
            sum(classify_anchor(h,t,base)=="explicit_gct" for h,t in parser.links),
        "full_original_GCT_retrieved":False,
        "source_url_or_study_anchor_exported":False,
    }

def http_get_page(url,opener=urllib.request.urlopen):
    if not is_safe_page(url):raise ValueError("untrusted page URL")
    req=urllib.request.Request(url,headers={"User-Agent":"HCC-stage36-public-GCT-source-audit/1",
                                            "Accept":"text/html"})
    with opener(req,timeout=25) as resp:
        raw=resp.read(MAX_HTML_BYTES+1)
        ctype=resp.headers.get("Content-Type","").lower()
        final=resp.geturl()
    if len(raw)>MAX_HTML_BYTES:raise ValueError("HTML source exceeds 8MiB")
    if "html" not in ctype or not is_safe_page(final):raise ValueError("not official HTML")
    return raw

def discover(pid,fetch=http_get_page):
    if pid not in PDC:raise ValueError("unsupported source study")
    out={"pdc_study_id":pid,"status":"SOURCE_GCT_ACCESS_HOLD",
         "GCT_original_file_verified":False,"quant_header_to_biological_aliquot_verified":False,
         "any_urls_or_sample_identifiers_exported":False}
    audits=[]
    # The PDC docs and portal are source-verified, not guessed file URLs.
    for url in (STUDY_PREFIX+pid,PORTAL,DOCUMENTATION):
        entry={"official_page_checked":True,"status":"ACCESS_HOLD"}
        try:
            raw=fetch(url)
            entry.update(page_inventory(raw,url,pid))
            entry["response_SHA256"]=hashlib.sha256(raw).hexdigest()
            entry["status"]="OFFICIAL_HTML_CHECKED"
        except Exception as exc:
            entry["error_class"]=type(exc).__name__
        audits.append(entry)
    out["pages"]=audits
    out["pages_successful_n"]=sum(x["status"]=="OFFICIAL_HTML_CHECKED" for x in audits)
    out["target_study_direct_GCT_href_candidates_n"]=sum(
        x.get("target_study_explicit_GCT_links_n",0) for x in audits)
    out["source_url_only_candidate_not_download_verified"]=True
    if out["pages_successful_n"]:
        out["status"]="PORTAL_HTML_ACCESS_VERIFIED_GCT_NOT_DOWNLOADED"
    return out

def parse_gct13_header(blob):
    """Parse full GCT header and ALL column metadata (not numerical matrix)."""
    if len(blob)>MAX_GCT_HEADER_BYTES:raise ValueError("header bytes exceed safety bound")
    data=blob.decode("utf-8-sig")
    lines=data.splitlines()
    if len(lines)<3 or lines[0].strip()!="#1.3":
        raise ValueError("GCT is not v1.3 or header incomplete")
    dim=lines[1].split("\t")
    if len(dim)!=4 or any(not re.fullmatch(r"\d+",v.strip()) for v in dim):
        raise ValueError("not GCT1.3 four dimensions")
    nr,nc,nrowmeta,ncolmeta=map(int,dim)
    if not 1<=nr<=500000 or not 1<=nc<=5000 or not 0<=nrowmeta<=40 or not 0<=ncolmeta<=150:
        raise ValueError("GCT dimensions outside allowed source scope")
    if len(lines)<3+ncolmeta:
        raise ValueError("column annotation section truncated")
    head=lines[2].split("\t")
    offset=1+nrowmeta
    if head[0].strip().lower()!="id" or len(head)!=offset+nc:
        raise ValueError("GCT header width mismatch")
    aliases=[s.strip() for s in head[offset:]]
    if any(not s for s in aliases) or len(set(aliases))!=nc:
        raise ValueError("nonunique GCT sample IDs")
    metadata={}
    for line in lines[3:3+ncolmeta]:
        row=line.split("\t")
        if len(row)!=len(head):
            raise ValueError("GCT column metadata width mismatch")
        name=row[0].strip()
        if not name or name in metadata or name.lower()=="id":
            raise ValueError("GCT metadata row name duplicate")
        metadata[name]=[x.strip() for x in row[offset:]]
    return {"sample_ids":aliases,"metadata":metadata,
            "declared_feature_rows_n":nr,"declared_samples_n":nc,
            "row_metadata_count_n":nrowmeta,"column_metadata_count_n":ncolmeta,
            "source_GCT_numeric_body_verified":False}

def alias_map(pdc_gct,official):
    """Strict per-column bijection candidate, never position-only mapping."""
    ids=pdc_gct["sample_ids"]
    metadata=pdc_gct["metadata"]
    required=("aliquot_id","aliquot_submitter_id","sample_type","case_id")
    missing=[key for key in required if key not in metadata]
    all_records=[]
    distinct_uuid=set()
    ambiguous=0
    missing_n=0
    record_by_column=[]
    for j,h in enumerate(ids):
        candidate_sets=[]
        for key in ("aliquot_id","aliquot_submitter_id"):
            a=metadata.get(key,[])
            if a and a[j] in official.get(key,{}):
                candidate_sets.append({rec["aliquot"] for rec in official[key][a[j]]})
        # GCT header itself might be the submitter identifier.
        if h in official.get("aliquot_submitter_id",{}):
            candidate_sets.append({rec["aliquot"] for rec in official["aliquot_submitter_id"][h]})
        if h in official.get("aliquot_id",{}):
            candidate_sets.append({rec["aliquot"] for rec in official["aliquot_id"][h]})
        if not candidate_sets:
            missing_n+=1
            record_by_column.append(None)
            continue
        union=set().union(*candidate_sets)
        if len(union)!=1 or any(x!=union for x in candidate_sets):
            ambiguous+=1
            record_by_column.append(None)
            continue
        uuid=next(iter(union))
        distinct_uuid.add(uuid)
        record_by_column.append(uuid)
        all_records.append(uuid)
    return {
        "gct_declared_sample_columns_n":len(ids),
        "mandatory_source_annotation_fields_present_n":sum(x in metadata for x in required),
        "mandatory_source_annotation_fields_missing_n":len(missing),
        "columns_with_no_exact_official_aliquot_id_n":missing_n,
        "columns_with_ambiguous_official_aliquot_n":ambiguous,
        "columns_with_one_unambiguous_official_aliquot_n":len(all_records),
        "distinct_official_aliquot_uuids_n":len(distinct_uuid),
        "duplicate_aliquot_columns_n":len(all_records)-len(distinct_uuid),
        "GCT_to_official_aliquot_330_candidate":
            len(ids)==330 and len(distinct_uuid)==330 and missing_n==ambiguous==0 and not missing,
        "GCT_to_quant_matrix_330_bridge_proven":False,
        "source_full_numeric_match_proven":False,
        "case_or_aliquot_identifiers_exported":False
    }

def compare_quant_column_names(parsed,quant_ids):
    """Exact alias-only evidence; this does not validate numerical identity."""
    ids=parsed["sample_ids"]
    aliases={s:i for i,s in enumerate(ids)}
    meta=parsed["metadata"]
    candidate_columns=collections.defaultdict(set)
    for h in ("aliquot_alias","aliquot_id","aliquot_submitter_id"):
        for j,v in enumerate(meta.get(h,[])):
            if v:
                candidate_columns[v].add(j)
    for j,v in enumerate(ids):candidate_columns[v].add(j)
    matched=[candidate_columns.get(a,set()) for a in quant_ids]
    unambiguous=[next(iter(x)) for x in matched if len(x)==1]
    return {
        "quant_header_columns_n":len(quant_ids),
        "quant_to_GCT_exact_alias_hits_n":sum(bool(x) for x in matched),
        "quant_to_GCT_unique_column_candidates_n":len(unambiguous),
        "quant_to_GCT_ambiguous_alias_candidates_n":sum(len(x)>1 for x in matched),
        "distinct_GCT_columns_for_unique_header_candidates_n":len(set(unambiguous)),
        "quant_to_GCT_330_candidate":
            len(quant_ids)==len(ids)==330 and len(unambiguous)==len(set(unambiguous))==330,
        "same_original_numeric_assay_certified":False,
        "raw_header_values_exported":False
    }

def review(payload):
    if payload.get("schema")!="hcc-stage36-original-pdc-GCT-discovery-v1":
        raise ValueError("source schema")
    vals=payload.get("studies")
    if not isinstance(vals,list) or {x.get("pdc_study_id") for x in vals}!=set(PDC):
        raise ValueError("unapproved source studies")
    for x in vals:
        for flag in ("GCT_original_file_verified","quant_header_to_biological_aliquot_verified",
                     "any_urls_or_sample_identifiers_exported"):
            if x.get(flag) is not False:raise ValueError("unverified scientific source assertion")
        if x.get("status") not in ("SOURCE_GCT_ACCESS_HOLD","PORTAL_HTML_ACCESS_VERIFIED_GCT_NOT_DOWNLOADED"):
            raise ValueError("unsupported source status")
        if len(x.get("pages",[]))!=3:raise ValueError("source page coverage changed")
        if any(q.get("source_url_or_study_anchor_exported") is True for q in x["pages"]):
            raise ValueError("source identity leak")
    return {
        "schema":"hcc-stage36-GCT-access-postflight-v1",
        "official_PDC_pages_accessible_study_n":sum(x["pages_successful_n"]>0 for x in vals),
        "GCT_original_file_download_verified_study_n":0,
        "GCT_original_to_quant_matrix_330_verified_study_n":0,
        "G30_B":"HOLD_GCT_ORIGINAL_SOURCE_AND_330_COLUMN_ALIQUOT_JOIN",
        "G30_C":"HOLD_NODE_RNA_SOURCE_QUANT",
        "G30_D":"STOP_NO_INDEPENDENT_MATCHED_RNA_PROTEIN",
        "novel_HCC_biology_claims_allowed":False,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    args=p.parse_args()
    payload={"schema":"hcc-stage36-original-pdc-GCT-discovery-v1",
        "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "studies":[discover(pid) for pid in PDC]}
    verdict=review(payload)
    for f,obj in ((args.out,payload),(args.review_out,verdict)):
        path=pathlib.Path(f);path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(obj,indent=2)+"\n",encoding="utf-8")
    print("STAGE36_GCT_ACCESS",json.dumps(payload))
    print("STAGE36_GCT_REVIEW",json.dumps(verdict))

if __name__=="__main__":main()
