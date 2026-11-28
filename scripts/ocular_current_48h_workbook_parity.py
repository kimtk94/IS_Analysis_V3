#!/usr/bin/env python3
"""Compare official current 48h workbook bytes and frozen selected-row values.
Never replace the historical archive or infer full historical equivalence."""
import argparse,collections,csv,hashlib,io,json,math,re,zipfile
from pathlib import Path
from urllib.parse import urlencode,unquote,urlparse
from xml.etree import ElementTree as ET
from publication_sprint import POOL_SHA,read_gene_table,sha
from publication_sprint_followup import compare
from ocular_p5j_official_TM_overlap import parse_official,TM_ALL,TM_HIGH,OFFICIAL_NAME
from ocular_frozen_241_external import PIN,read_panel
import ocular_final_publication_inventory as cloud

OLD="308ffad592cbe1eb5262c5bec1df467257f8859b4905198ad3db9d5513b1c818"
def req(ok,msg):
    if not ok:raise ValueError(msg)
def member(path):
    with zipfile.ZipFile(path) as z:
        names=[n for n in z.namelist() if n.rsplit("/",1)[-1]==OFFICIAL_NAME]
        req(len(names)==1 and z.getinfo(names[0]).file_size<=20000000,"Exact workbook member")
        return z.read(names[0])
def digest(raw):return hashlib.sha256(raw).hexdigest()
def cloud_workbooks():
    ns={"s":"http://s3.amazonaws.com/doc/2006-03-01/"}
    raw,src=cloud.fetch(cloud.BASE+"?"+urlencode({"list-type":"2","prefix":"PMC12882097.","delimiter":"/","max-keys":20}),200000)
    xml=ET.fromstring(raw)
    req(xml.findtext("s:IsTruncated",namespaces=ns)=="false","Incomplete versions")
    prefixes=[n.text for n in xml.findall("s:CommonPrefixes/s:Prefix",ns)]
    req(0<len(prefixes)<=4 and all(re.fullmatch(r"PMC12882097\.[1-9][0-9]*/",p) for p in prefixes),"Version inventory")
    out={"listing":src,"versions":[]}
    for prefix in prefixes:
        try:
            raw,meta=cloud.fetch(cloud.BASE+"metadata/"+prefix[:-1]+".json",500000)
            d=json.loads(raw)
            req(d["pmcid"]=="PMC12882097" and str(d["version"])==prefix[:-1].split(".")[-1] and d["doi"]=="10.1167/iovs.67.2.13" and d.get("is_retracted") is False and d.get("is_pmc_openaccess") is True and d.get("license_code"),"Article identity/license")
            media=[v for v in d.get("media_urls",[]) if unquote(urlparse(v).path.rsplit("/",1)[-1])==OFFICIAL_NAME]
            req(len(media)==1,"Exact public workbook")
            raw,source=cloud.object_read(media[0],prefix,20000000)
            out["versions"].append({"status":"PUBLIC_WORKBOOK_VERIFIED","prefix":prefix,"metadata":meta,"license":d["license_code"],"workbook":source})
        except (OSError,ValueError,KeyError,TypeError) as ex:out["versions"].append({"status":"HOLD_SOURCE","prefix":prefix,"reason":str(ex)[:180]})
    return out
def archived_rows(path,name):
    with zipfile.ZipFile(path) as z:
        names=[n for n in z.namelist() if n.rsplit("/",1)[-1]==name]
        req(len(names)==1,"Archived table identity")
        rows=list(csv.DictReader(io.StringIO(z.read(names[0]).decode()),delimiter="\t"))
        req(len(rows)==len({v["gene_id"] for v in rows}),"Archived gene duplicates")
        return {v["gene_id"]:v for v in rows}
def equivalent(old,new):
    if set(old)!=set(new):return {"keys_match":False,"old_only":len(set(old)-set(new)),"new_only":len(set(new)-set(old)),"changed_values":None}
    diff=[]
    for k in sorted(old):
        a,b=old[k],new[k]
        req(set(a)==set(b),"Selected table schema")
        for f in a:
            if isinstance(b[f],float):same=float(a[f])==b[f]
            else:same=a[f]==str(b[f])
            if not same:diff.append({"gene_id":k,"field":f})
    return {"keys_match":True,"rows":len(old),"changed_values":len(diff),"first_differences":diff[:20]}
def run(a):
    previous=json.loads(a.hold.read_text())
    req(previous["result"]["selection_48h"]["status"]=="HOLD_SOURCE_HASH_CHANGED" and previous["provenance"]["run_id"]==45578789122,"Prior HOLD")
    x=json.loads(a.targets.read_text())["result"]
    req(x["status"]=="PUBLISHED_COUNTS_RECONCILED_EXACT_ID_OVERLAP" and len(set(x["source_gene_ids"]))==67,"Exact 67 targets")
    ids=x["source_gene_ids"]
    req(sha(a.panel)==PIN and sha(a.pool)==POOL_SHA,"Frozen archives")
    panel=read_panel(a.panel)
    src={"zip_sha256":sha(a.external),"historical_zip_sha256":OLD,"same_historical_zip":sha(a.external)==OLD,"workbook_sha256":digest(member(a.external))}
    try:src["official_cloud"]=cloud_workbooks()
    except (OSError,ValueError,KeyError,TypeError,ET.ParseError) as ex:src["official_cloud"]={"status":"HOLD_SOURCE","reason":str(ex)[:180]}
    matches=[v["prefix"] for v in src["official_cloud"].get("versions",[]) if v.get("workbook",{}).get("sha256")==src["workbook_sha256"]]
    src["matching_current_cloud_versions"]=matches
    out={"status":"CURRENT_SOURCE_COMPARISON_NOT_HISTORICAL_EQUIVALENCE","source":src,"historical_byte_identity_status":"MATCH" if src["same_historical_zip"] else "HOLD_CHANGED_ZIP","limits":["Compare exact current official workbook member and archived selected rows; prior full workbook byte identity not established. Current-source parity cannot establish unchanged historical values for unarchived genes.","All 67 targets retained with original model coverage and current selection-stage values; no new panel refit, independent validation, enrichment or causal claims. Donor/pool share cohort."]}
    external=parse_official(a.external)
    ext=external[TM_ALL]
    with zipfile.ZipFile(a.pool) as z:
        donor=read_gene_table(z,"p5_all_gene_results.tsv.gz","logFC","FDR",16050)
        pool=read_gene_table(z,"p5d_poolunit_all_gene_results.tsv.gz","pool_logFC","pool_FDR",21081)
    shared,rows=compare(donor,pool,ext)
    strong={z["gene_id"]:z for z in rows if z["strong_concordant"]}
    discord={z["gene_id"]:z for z in rows if z["sign_discordant"]}
    st=equivalent(archived_rows(a.panel,"strong_241.tsv"),strong)
    di=equivalent(archived_rows(a.panel,"discordant_64.tsv"),discord)
    out["current_counts"]={"shared":shared,"external_fdr_available":len(rows),"strong":len(strong),"discordant":len(discord),"high_confidence_table":len(external[TM_HIGH])}
    out["archived_row_parity"]={"strong241":st,"discordant64":di}
    parity=out["current_counts"]|{}
    good=(shared,len(rows),len(strong),len(discord))==(12251,12168,241,64) and st.get("changed_values")==di.get("changed_values")==0 and set(strong)==set(panel)
    out["current_source_status"]="CURRENT_OFFICIAL_WORKBOOK_AND_ARCHIVED_ROWS_RECONCILED" if good and matches else "HOLD_CURRENT_SOURCE_OR_SELECTED_PARITY"
    rec=[]
    for g in ids:
        u,v,w=donor.get(g),pool.get(g),ext.get(g)
        internal=u is not None and v is not None and u[1]<.05 and v[1]<.05 and u[0]*v[0]>0
        eligible=internal and min(abs(u[0]),abs(v[0]))>=1
        status="NOT_INTERNAL_STRONG"
        if eligible:
            status="CURRENT_48H_NOT_TESTED" if w is None else "CURRENT_48H_FDR_MISSING" if w["fdr"] is None else "CURRENT_48H_FDR_NOT_SIGNIFICANT" if w["fdr"]>=.05 else "CURRENT_48H_DIRECTION_NOT_CONCORDANT" if u[0]*w["log2fc"]<=0 else "CURRENT_ALL_SELECTION_RULES_PASS"
        rec.append({"gene_id":g,"donor":dict(symbol=u[2],log2fc=u[0],fdr=u[1]) if u else None,"pool":dict(symbol=v[2],log2fc=v[0],fdr=v[1]) if v else None,"current_48h":w,"in_current_high_confidence_table":g in external[TM_HIGH],"status":status})
    out["target_counts"]={"current_48h_tested":sum(v["current_48h"] is not None for v in rec),"current_48h_fdr05":sum(v["current_48h"] is not None and v["current_48h"]["fdr"] is not None and v["current_48h"]["fdr"]<.05 for v in rec),"internal_strong_4":sum(v["status"]!="NOT_INTERNAL_STRONG" for v in rec),"statuses":dict(collections.Counter(v["status"] for v in rec))}
    out["internal_strong_details"]=[v for v in rec if v["status"]!="NOT_INTERNAL_STRONG"]
    out["records"]=rec
    return out
def main():
    p=argparse.ArgumentParser()
    for k in ("hold","targets","panel","pool","external","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    try:x=run(a)
    except (OSError,ValueError,KeyError,TypeError,zipfile.BadZipFile) as ex:x={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(ex)[:250]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps({k:v for k,v in x.items() if k!="records"},indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
