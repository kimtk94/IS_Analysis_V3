#!/usr/bin/env python3
"""2026 NCBI PMC Cloud article-media manifest and bounded XLSX verification.

Per NCBI 2026 official PMC Cloud README, article versions are stored at
s3://pmc-oa-opendata/PMC<id>.<version> and metadata at
metadata/PMC<id>.<version>.json. Read metadata first; resolve media_urls
*from NCBI metadata*, never guess an S3 object or assume old PMC HTTP works.

No patient identifiers, gene values, raw/signed URLs or full workbooks
are persisted. Public 2023 HCC supplementary mmc3.xlsx remains only a
candidate until a validated OOXML workbook is actually downloaded.
"""
import argparse
import collections
import datetime as dt
import hashlib
import io
import math
import statistics
import json
import pathlib
import re
import urllib.error
import urllib.parse
import urllib.request
import zipfile
import xml.etree.ElementTree as ET

BASE="https://pmc-oa-opendata.s3.amazonaws.com/"
PMCID="PMC10783603"
DOI="10.1016/j.xcrm.2023.101315"
TARGET="mmc3.xlsx"
MAX_META=2*1024*1024
MAX_XLSX=24*1024*1024
VERSIONS=(1,2,3)
NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

def official(url):
    p=urllib.parse.urlsplit(url)
    return (p.scheme=="https" and p.hostname=="pmc-oa-opendata.s3.amazonaws.com"
            and p.path.startswith(("/metadata/"+PMCID+".","/"+PMCID+"."))
            and not p.username and not p.password)

def resolve_cloud_uri(raw):
    """Translate NCBI-issued s3:// bucket URIs to documented anonymous HTTPS."""
    p=urllib.parse.urlsplit(raw)
    if p.scheme=="s3" and p.hostname=="pmc-oa-opendata":
        if not p.path.startswith("/"+PMCID+".") or "/../" in p.path:
            raise ValueError("unapproved S3 media key")
        # md5 is publisher object provenance, not a required S3 GET parameter.
        value=BASE+p.path.lstrip("/")
        if not official(value):raise ValueError("cloud key not in allowed PMC prefix")
        return value
    if official(raw):
        return raw
    raise ValueError("media URL not within official PMC Cloud")

def bounded(url,cap,timeout=12):
    if not official(url):raise ValueError("non-official cloud URL")
    req=urllib.request.Request(url,headers={
        "User-Agent":"HCC-PMC-2026-OpenData-Verification/1.0",
        "Accept":"application/json,application/octet-stream,*/*"})
    with urllib.request.urlopen(req,timeout=timeout) as res:
        if not official(res.url):raise ValueError("redirect outside official cloud")
        b=res.read(cap+1)
        if len(b)>cap:raise ValueError("size limit exceeded")
        return b,res.status,res.headers.get("Content-Type","").split(";")[0]

def parse_article_metadata(data):
    meta=json.loads(data)
    if str(meta.get("pmcid","")).upper()!=PMCID:
        raise ValueError("PMCID mismatch")
    if meta.get("doi") and str(meta["doi"]).lower()!=DOI:
        raise ValueError("DOI mismatch")
    raw=meta.get("media_urls") or []
    if not isinstance(raw,list):raise ValueError("media_urls type invalid")
    candidate_urls=[]
    xlsx_filenames=[]
    entry_shapes=collections.Counter()
    link_hosts=collections.Counter()
    public_suffixes=collections.Counter()
    dictionary_fields=set()
    for entry in raw:
        entry_shapes[type(entry).__name__]+=1
        if isinstance(entry,str):url=entry
        elif isinstance(entry,dict):
            dictionary_fields.update(str(x)[:60] for x in entry.keys())
            url=entry.get("url") or entry.get("s3_url") or entry.get("media_url") or ""
        else:continue
        link_hosts[urllib.parse.urlsplit(url).hostname or "NO_HOST"]+=1
        suffix=urllib.parse.unquote(urllib.parse.urlsplit(url).path).lower().rsplit(".",1)
        if len(suffix)==2 and len(suffix[-1])<=12:
            public_suffixes[suffix[-1]]+=1
        try:resolved=resolve_cloud_uri(url)
        except ValueError:continue
        filename=urllib.parse.unquote(urllib.parse.urlsplit(url).path).rsplit("/",1)[-1]
        if re.fullmatch(r"[a-zA-Z0-9_.-]{1,120}\.xlsx", filename, re.I):
            xlsx_filenames.append(filename)
        if re.search(r"(?:^|[-_.])mmc3\.xlsx$",filename,re.I):
            candidate_urls.append(resolved)
    return {
        "license_code":str(meta.get("license_code") or "UNKNOWN")[:32],
        "open_access":meta.get("is_pmc_openaccess"),
        "article_version":meta.get("version"),
        "media_url_count":len(raw),
        "media_entry_type_counts":dict(entry_shapes),
        "media_entry_public_keys":sorted(dictionary_fields),
        "media_url_host_counts":dict(link_hosts),
        "media_file_extension_counts":dict(public_suffixes),
        "official_xlsx_media_filenames":sorted(set(xlsx_filenames)),
        "matching_file_count":len(candidate_urls),
        "candidate_urls_in_memory":candidate_urls,
    }

def sheet_shape(z,part,first_rows=8):
    """Bounded first rows and dimensionality; no cell contents are exported."""
    row_summaries=[]
    dimension=None
    try:
        with z.open(part) as stream:
            for event,elem in ET.iterparse(stream,events=("end",)):
                if elem.tag==NS+"dimension":
                    dimension=elem.get("ref")
                if elem.tag==NS+"row":
                    cells=elem.findall(NS+"c")
                    nc=sum(1 for c in cells if c.find(NS+"v") is not None)
                    numeric=sum(1 for c in cells if c.get("t") not in ("s","str","inlineStr","b","e") and c.find(NS+"v") is not None)
                    row_summaries.append({"number_of_cells":len(cells),
                                          "numeric_cells":numeric,
                                          "nonempty_value_elements":nc})
                    elem.clear()
                    if len(row_summaries)>=first_rows:break
                elif elem.tag not in (NS+"c",NS+"v",NS+"is",NS+"t"):
                    # Avoid excessive XML memory while keeping child cell data.
                    if elem.tag!=NS+"worksheet":elem.clear()
    except Exception as exc:
        return {"status":"UNVERIFIED","error_type":type(exc).__name__}
    # Dimension may not appear in all producer XML. A row prefix is not all rows.
    return {"status":"FIRST_ROWS_SCANNED","dimension_ref":dimension,
            "rows_inspected":len(row_summaries),
            "first_row_width":row_summaries[0]["number_of_cells"] if row_summaries else None,
            "row_summaries":row_summaries,
            "header_identifiers_and_expression_values_logged":False}

def read_shared_strings(z,limit=125000):
    """Keep publisher cell text in memory only; never store or print identifiers."""
    path="xl/sharedStrings.xml"
    if path not in z.namelist():return []
    if z.getinfo(path).file_size>36*1024*1024:
        raise ValueError("shared strings XML too large for privacy-safe sampling")
    result=[]
    with z.open(path) as stream:
        for _,el in ET.iterparse(stream,events=("end",)):
            if el.tag==NS+"si":
                val="".join((t.text or "") for t in el.iter(NS+"t"))
                result.append(val)
                el.clear()
                if len(result)>limit:raise ValueError("too many shared string values")
    return result

def classify_patient_label(x):
    """Conservative candidate T/N matching; unclassified is not a patient."""
    val=str(x).strip()
    if len(val)>90 or len(val)<2:return None
    # Common forms: T_123, N_123, 123_T, 123_N, PT123-Tumor, PT123-Normal
    patterns=(
       r"^(tumou?r|normal|t|n)[-_ .]+([a-zA-Z0-9][a-zA-Z0-9_-]{1,60})$",
       r"^([a-zA-Z0-9][a-zA-Z0-9_-]{1,60})[-_ .]+(tumou?r|normal|t|n)$",
       r"^(t|n)([0-9]{1,10})$",
       r"^([0-9]{1,10})(t|n)$",
    )
    for i,pattern in enumerate(patterns):
        m=re.match(pattern,val,re.I)
        if m:
            kind=m.group(1) if i in (0,2) else m.group(2)
            key=m.group(2) if i in (0,2) else m.group(1)
            tag="T" if kind.lower() in ("t","tumor","tumour") else "N"
            return (key.lower(),tag,i)
    return None

def sample_header_audit(z,sheet_path,expected_sample_columns=304):
    """Read row-2 labels; aggregate T/N pairs without outputting patient IDs."""
    shared=read_shared_strings(z)
    headers=[]
    with z.open(sheet_path) as stream:
        for _,row in ET.iterparse(stream,events=("end",)):
            if row.tag!=NS+"row":continue
            if str(row.get("r"))!="2":
                row.clear()
                continue
            for c in row.findall(NS+"c"):
                typ=c.get("t","n")
                v=c.find(NS+"v")
                x=None if v is None else v.text
                if typ=="s" and x is not None:
                    i=int(x)
                    x=shared[i] if 0<=i<len(shared) else ""
                elif typ=="inlineStr":
                    x="".join(y.text or "" for y in c.iter(NS+"t"))
                headers.append("" if x is None else str(x))
            break
    if not headers:raise ValueError("no sample header row 2")
    # Three feature annotation columns are expected for 6512 proteins.
    # Treat that as a provisional hypothesis, not a verified patient count.
    candidate_headers=headers[3:]
    out={"header_row_excel":2,"observed_row2_header_cells":len(headers),
         "annotation_columns_excluded_hypothesis":3,
         "candidate_sample_columns":len(candidate_headers),
         "matches_expected_304_columns":len(candidate_headers)==expected_sample_columns,
         "shared_string_records_loaded":len(shared),
         "patient_ids_or_labels_logged":False,
         "tumor_normal_pairing_verified":False}
    paired={}
    by_pattern={}
    for label in candidate_headers:
        c=classify_patient_label(label)
        if c is None:continue
        key,kind,method=c
        paired.setdefault(key,set()).add(kind)
        by_pattern[method]=by_pattern.get(method,0)+1
    out["candidate_id_pattern_matches"]=sum(by_pattern.values())
    out["candidate_header_match_patterns"]=by_pattern
    out["candidate_keys_with_T_and_N"]=sum(types=={"T","N"} for types in paired.values())
    out["candidate_keys_T_only"]=sum(types=={"T"} for types in paired.values())
    out["candidate_keys_N_only"]=sum(types=={"N"} for types in paired.values())
    out["unclassified_headers"]=len(candidate_headers)-sum(by_pattern.values())
    out["interpretation"]="Heuristic T/N label grouping; requires author metadata/manual verification before scientific PASS."
    return out

def _excel_column_number(ref):
    letters=re.match(r"^[A-Z]+",ref.upper())
    if not letters:return 0
    n=0
    for ch in letters.group():
        n=n*26+ord(ch)-64
    return n

def protein_numeric_coverage(z,part,expected_sample_columns=304):
    """Process actual numeric protein values, emit only aggregated coverage."""
    observed_sample_counts=[0]*expected_sample_columns
    features=[]
    row_count=0
    with z.open(part) as stream:
        for _,row in ET.iterparse(stream,events=("end",)):
            if row.tag!=NS+"row":continue
            row_num=int(row.get("r") or "0")
            if row_num<=2:
                row.clear()
                continue
            count=0
            for cell in row.findall(NS+"c"):
                idx=_excel_column_number(cell.get("r",""))
                if not 4<=idx<4+expected_sample_columns:continue
                if cell.get("t") in ("s","str","inlineStr","e","b"):continue
                v=cell.find(NS+"v")
                if v is None or v.text is None:continue
                try:
                    x=float(v.text)
                except ValueError:continue
                if not math.isfinite(x):continue
                count+=1
                observed_sample_counts[idx-4]+=1
            features.append(count)
            row_count+=1
            row.clear()
    if not row_count:raise ValueError("No protein feature rows")
    values=sum(features)
    available=row_count*expected_sample_columns
    ordered=sorted(observed_sample_counts)
    return {
        "status":"NUMERIC_VALUE_COVERAGE_AUDITED",
        "number_of_protein_feature_rows":row_count,
        "assayed_sample_columns":expected_sample_columns,
        "total_possible_expression_cells":available,
        "actual_finite_numeric_expression_cells":values,
        "numeric_fraction":round(values/available,6),
        "sample_valid_protein_count_min":min(ordered),
        "sample_valid_protein_count_median":statistics.median(ordered),
        "sample_valid_protein_count_p10":ordered[(len(ordered)-1)//10],
        "sample_valid_protein_count_max":max(ordered),
        "protein_feature_valid_sample_count_median":statistics.median(features),
        "protein_features_with_at_least_152_samples":sum(x>=152 for x in features),
        "protein_features_with_all_304_samples":sum(x==expected_sample_columns for x in features),
        "patient_ids_or_numeric_values_logged":False,
        "rna_protein_join_confirmed":False,
        "note":"Quantification completeness in one protein-only cohort; clinical matching and RNA values remain unverified."
    }

def xlsx_only_manifest(raw):
    if not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("not XLSX")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if "xl/workbook.xml" not in z.namelist():raise ValueError("no workbook")
        if z.getinfo("xl/workbook.xml").file_size>MAX_META:raise ValueError("oversize xml")
        t=ET.fromstring(z.read("xl/workbook.xml"))
        sh=t.find(NS+"sheets")
        names=[s.get("name") for s in sh] if sh is not None else []
        xml_parts=sorted(
            (k for k in z.namelist() if re.fullmatch(r"xl/worksheets/sheet\d+\.xml",k)),
            key=lambda k:int(re.search(r"sheet(\d+)",k).group(1)))
        profiles=[]
        for i,part in enumerate(xml_parts):
            info=z.getinfo(part)
            profiles.append({
                "sheet_index":i+1,
                "sheet_title":names[i] if i<len(names) else "UNKNOWN",
                "uncompressed_xml_bytes":info.file_size,
                **sheet_shape(z,part)
            })
        coverage={"status":"UNVERIFIED"}
        if len(xml_parts)>=2 and len(names)>=2 and "6512 proteins" in names[1]:
            try:coverage=protein_numeric_coverage(z,xml_parts[1])
            except Exception as exc:coverage={"status":"UNVERIFIED","error_type":type(exc).__name__}
        patient_header={"status":"UNVERIFIED"}
        selected=[part for idx,part in enumerate(xml_parts) if idx<len(names) and "6512 proteins" in names[idx]]
        if len(selected)==1:
            try:
                patient_header=sample_header_audit(z,selected[0])
                patient_header["status"]="HEADER_ROW_STRUCTURE_AUDITED"
            except Exception as exc:
                patient_header={"status":"UNVERIFIED","error_type":type(exc).__name__}
        return {"xlsx_verified":True,"sheet_count":len(names),"sheet_names":names,
                "worksheet_structures":profiles,
                "protein_sample_header_audit":patient_header,
                "protein_numeric_coverage":coverage,
                "patient_expression_values_or_raw_ids_logged":False}


def audit():
    r={"pmcid":PMCID,"doi":DOI,
       "target":TARGET,"status":"UNVERIFIED",
       "official_source":"https://pmc-oa-opendata.s3.amazonaws.com/README.txt",
       "metadata_attempts":[],"download_verified":False,
       "same_patient_rna_protein_join":False}
    for version in VERSIONS:
        meta_url=BASE+"metadata/"+PMCID+"."+str(version)+".json"
        attempt={"version":version}
        try:
            b,status,kind=bounded(meta_url,MAX_META)
            info=parse_article_metadata(b)
            candidates=info.pop("candidate_urls_in_memory")
            attempt.update({"http":status,"content_type":kind,**info})
            r["metadata_attempts"].append(attempt)
            if not candidates:continue
            for url in candidates:
                try:
                    raw,code,ctype=bounded(url,MAX_XLSX,timeout=22)
                    out=xlsx_only_manifest(raw)
                    r.update({"status":"S3_XLSX_MANIFEST_VERIFIED",
                              "download_verified":True,
                              "source_article_version":version,
                              "download_bytes":len(raw),
                              "sha256":hashlib.sha256(raw).hexdigest(),
                              "content_type":ctype,
                              "manifest":out})
                    return r
                except urllib.error.HTTPError as exc:
                    attempt["media_http_error_code"]=exc.code
                except Exception as exc:
                    attempt["media_error_type"]=type(exc).__name__
        except urllib.error.HTTPError as exc:
            attempt.update({"http":exc.code,"metadata_status":"HTTP_ERROR"})
            r["metadata_attempts"].append(attempt)
        except Exception as exc:
            attempt.update({"metadata_status":"UNVERIFIED","error_type":type(exc).__name__})
            r["metadata_attempts"].append(attempt)
    r["status"]="PMC_2026_CLOUD_NOT_VERIFIED"
    return r

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="results/hcc/feasibility/pmc2026_cloud_mmc3_audit.json")
    args=p.parse_args()
    report={"checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "audit":audit(),"V02":"BLOCKED"}
    path=pathlib.Path(args.out)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False))

if __name__=="__main__":main()
