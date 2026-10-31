#!/usr/bin/env python3
"""Public-source-only study provenance audit for ocular manuscript.

No arbitrary LLM inferences, private data, or whole copyrighted article text.
Attempts official DOI/PubMed/EuropePMC/PMC for 2026 Mehrotra erratum,
and official EuropePMC search for GSE313736 original paper/gene list.
CI success is NOT scientific GO: inspect status, provenance and excerpts.
"""
import argparse,html,json,re,time,urllib.parse,urllib.request,xml.etree.ElementTree as ET
from pathlib import Path

ROOT="https://www.ebi.ac.uk/europepmc/webservices/rest"
ERRATUM_PMID="42895266"
ERRATUM_DOI="10.1167/iovs.67.12.34"
ORIGINAL_PMID="41636427"
GSE="GSE313736"
CAP=3_500_000

def fetch(url,cap=CAP):
    last=""
    for attempt in range(2):
        try:
            request=urllib.request.Request(url,headers={
                "User-Agent":"ocular-master-degree-public-provenance/1.1 (noncommercial research)",
                "Accept":"application/json, application/xml, text/xml, text/html;q=0.9, */*;q=0.7"})
            with urllib.request.urlopen(request,timeout=32) as f:
                payload=f.read(cap+1)
                url_final=f.url;ctype=f.headers.get("Content-Type","")
            if len(payload)>cap:raise ValueError("RESPONSE_ABOVE_BOUND")
            if not payload:raise ValueError("EMPTY_RESPONSE")
            return {"bytes":payload,"final_url":url_final,"content_type":ctype}
        except (OSError,TimeoutError,ValueError) as e:
            last=type(e).__name__+": "+str(e)[:200]
            if attempt==0:time.sleep(2)
    raise RuntimeError(last)

def safeget(url,label,store,cap=CAP):
    try:
        f=fetch(url,cap=cap)
        store.append({"label":label,"url":url,"status":"FETCHED",
           "response_bytes":len(f["bytes"]),"final_url":f["final_url"],
           "content_type":f["content_type"]})
        return f["bytes"]
    except Exception as e:
        store.append({"label":label,"url":url,"status":"UNAVAILABLE",
                     "reason":str(e)[:220]})
        return None

def parse_xml(data):
    try:return ET.fromstring(data)
    except ET.ParseError:return None

def clean(s):
    return re.sub(r"\s+"," ",html.unescape(re.sub("<[^>]*>"," ",s))).strip()

def snippets(s,terms,max_results=8):
    out=[]
    # Preserve only short limited excerpts, no wholesale redistribution.
    normalized=clean(s)
    for keyword in terms:
        for m in list(re.finditer(re.escape(keyword),normalized,re.I))[:2]:
            low=max(0,m.start()-160);high=min(len(normalized),m.end()+270)
            value=normalized[low:high]
            if value not in out:out.append(value)
            if len(out)>=max_results:return out
    return out

def query(q,attempts):
    url=ROOT+"/search?"+urllib.parse.urlencode(
        {"query":q,"format":"json","pageSize":15})
    raw=safeget(url,"europepmc-search-"+q[:45],attempts,cap=1_500_000)
    if raw is None:return []
    try:
        ob=json.loads(raw)
        return ob.get("resultList",{}).get("result",[])
    except (json.JSONDecodeError,TypeError):return []

def eutils_xml(pmid,attempts):
    u="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"+
      urllib.parse.urlencode({"db":"pubmed","id":pmid,"rettype":"xml","retmode":"xml"})
    return safeget(u,"NCBI-PubMed-EFetch-"+pmid,attempts,cap=1_500_000)

def audit(out):
    attempts=[]
    result={"status":"HOLD_ERRATUM_CONTENT_NOT_AUTHENTICATED",
       "pmid":ERRATUM_PMID,"doi":ERRATUM_DOI,
       "study":"Mehrotra et al., 2026 IOVS 100nM 48h TM",
       "publication_link_verified":False,"erratum_fulltext_recovered":False,
       "correction_categorization":"UNKNOWN",
       "data_values_affected":"UNKNOWN",
       "analysis_should_not_assume_no_data_change":True,
       "attempts":attempts}
    err=query("EXT_ID:"+ERRATUM_PMID+" AND SRC:MED",attempts)
    if not err:err=query("DOI:"+ERRATUM_DOI,attempts)
    result["europepmc_erratum_entries"]=[{
        k:z.get(k) for k in ("id","pmid","pmcid","doi","title","journalTitle",
                            "firstPublicationDate","isOpenAccess","hasPDF","hasTextMinedTerms")}
        for z in err[:4]]
    matching=[x for x in err if x.get("pmid")==ERRATUM_PMID or x.get("doi")==ERRATUM_DOI]
    result["publication_link_verified"]=bool(matching)
    raw=eutils_xml(ERRATUM_PMID,attempts)
    if raw:
        tree=parse_xml(raw)
        if tree is not None:
            # Explicit source must pair erratum to original, not merely same title.
            reftexts=[" ".join(x.itertext()) for x in tree.findall(".//CommentsCorrections")]
            result["pubmed_comments_corrections"]=[clean(x)[:400] for x in reftexts[:5]]
            result["pubmed_records_link_to_original"]=any(
              ORIGINAL_PMID in x or "10.1167/iovs.67.2.13" in x for x in reftexts)
    pmcids=[x.get("pmcid") for x in matching if x.get("pmcid")]
    # PubMed/EuropePMC can contain no PMCID for publisher-only corrections.
    sources=[
       ("IOVS-publisher-articleid","https://iovs.arvojournals.org/article.aspx?articleid=2818715"),
       ("IOVS-publisher-doi","https://iovs.arvojournals.org/doi/10.1167/iovs.67.12.34"),
       ("IOVS-corrected-DOI","https://doi.org/10.1167/iovs.67.12.34")]
    for pmcid in pmcids:
        sources.insert(0,("EuropePMC-official-fullTextXML",ROOT+"/"+pmcid+"/fullTextXML"))
    accessed=[]
    for label,u in sources:
        body=safeget(u,label,attempts,cap=2_500_000)
        if body is None:continue
        text=body.decode("utf8","replace")
        norm=clean(text)
        # Some web servers return a captcha/consent page with HTTP 200;
        # do NOT turn that into "correction does not affect data."
        if any(x in norm.lower() for x in ("captcha","checking your browser","access denied")):
            continue
        is_suitable=((ERRATUM_DOI.lower() in norm.lower() or
            ("erratum" in norm.lower() and "trabecular meshwork" in norm.lower()))
           and len(norm)>150)
        if not is_suitable:continue
        excerpt=snippets(norm,["When the article","When this article",
             "originally published","inadvertently","error in","corrected",
             "correction","supplementary","figure","gene"],max_results=7)
        accessed.append({"source":label,"url":u,
          "normalized_length":len(norm),"excerpt":excerpt})
    result["publisher_or_fulltext_records"]=accessed
    # A full erratum can be available but its meaning must be read and checked.
    result["erratum_fulltext_recovered"]=any(
        x["source"] in ("EuropePMC-official-fullTextXML","IOVS-publisher-articleid",
                        "IOVS-publisher-doi","IOVS-corrected-DOI") and
        bool(x["excerpt"]) for x in accessed)
    if result["erratum_fulltext_recovered"]:
        result["status"]="REVIEW_ERRATUM_FULLTEXT_FOUND_IMPACT_NEEDS_MANUAL_COMPARISON"
    else:
        result["status"]="HOLD_ERRATUM_FULLTEXT_UNAVAILABLE_OR_UNVERIFIABLE"
    study_queries=[
        'GSE313736',
        '"Single Cell Profiling of Human Trabecular Meshwork Cells"',
        '"trabecular meshwork" "74 donors" "dexamethasone"']
    study=[]
    for q in study_queries:
        arr=query(q,attempts)
        study.append({"query":q,"matches":[{k:r.get(k) for k in
         ("title","authorString","firstPublicationDate","doi","pmid","pmcid",
          "journalTitle","pubType","isPreprint")}
          for r in arr[:10]]})
    result2={"status":"HOLD_AUTHOR_80_GENE_IDS_NOT_FETCHED",
       "dataset":GSE,"study_searches":study,
       "original_abstract_degs_up":63,"original_abstract_degs_down":17,
       "individual_80_author_gene_ids_verified":False,
       "no_guessing_gene_list_from_83_or_88_project_genes":True,
       "official_series_supplement_inventory_previously_checked":True,
       "search_queries_are_not_proof_of_nonpublication":True,
       "attempts":attempts}
    source_targets=[r for item in study for r in item["matches"] if
      r.get("title") and any(term in r["title"].lower() for term in
       ("single cell profiling of human trabecular",
        "single-cell profiling of human trabecular"))]
    result2["possible_original_paper_records"]=source_targets
    if source_targets:
        result2["status"]="REVIEW_POSSIBLE_AUTHOR_PUBLICATION_FOUND_NEEDS_SUPPLEMENT"
    out.mkdir(parents=True,exist_ok=True)
    (out/"MEHROTRA_ERRATUM_OFFICIAL_SOURCE_AUDIT.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    (out/"GSE313736_AUTHOR80_PUBLICATION_SEARCH.json").write_text(
        json.dumps(result2,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("ERRATUM",result["status"],"publication_record",result["publication_link_verified"],
      "fulltext_found",result["erratum_fulltext_recovered"])
    print("GSE313736",result2["status"],"possible_papers",len(source_targets))
    print("SOURCE_ATTEMPTS",json.dumps(attempts,ensure_ascii=False))
    # Authenticated failure to retrieve is a scientific HOLD, not CI false success.
    if not matching and not raw:raise SystemExit("MISSING_BOTH_OFFICIAL_ERRATUM_METADATA_SOURCES")
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    audit(a.out)
