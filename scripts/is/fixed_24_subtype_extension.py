#!/usr/bin/env python3
"""Extend the immutable 24 nonpalindromic EAS regional markers to subtype GWAS.
Retrospective same-release comparisons; fixed 192-test family including missing.
No new marker selection, causal gene inference or independent replication.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re,sys,time
from pathlib import Path
from urllib.parse import urljoin
import yaml
from pyliftover import LiftOver
import importlib
legacy_helpers=importlib.import_module('scripts.is.audit_gws_eas_24_leads_original')
read_legacy_180,compare=legacy_helpers.read_legacy_180,legacy_helpers.compare
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"indexed_helpers/scripts/is"))
import fixed_subtype_indexed as base

CHAIN_SHA="5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d"
TRAITS={"AIS":"Ischemic stroke","CES":"Cardioembolic stroke","LAS":"Large artery stroke","SVS":"Small vessel stroke"}
STUDIES=base.STUDIES
FAMILY=192
def require(ok,msg):
    if not ok:raise ValueError(msg)
def sha(p):
    with p.open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()
def write(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def gate_metadata(md,acc,study,name):
    a,t=study.split("_")
    require(md.get("gwas_id")==acc and md.get("data_file_name")==name,"GWAS accession/filename mismatch")
    require(md.get("genome_assembly")=="GRCh38" and md.get("is_harmonised") is True,"Harmonised GRCh38 unavailable")
    require(md.get("trait_description")==[TRAITS[t]],"Unexpected phenotype")
    require({v for row in md.get("samples",[]) for v in row.get("sample_ancestry_category",[])}=={"European" if a=="EUR" else "East Asian"},"Unexpected ancestry")
def panel(legacy,chain):
    require(sha(chain)==CHAIN_SHA,"Previously verified chain differs")
    old=read_legacy_180(legacy)
    lo=LiftOver(str(chain))
    out={}
    for a in old["selected"]:
        rid=a["region_id"];ch=a["chromosome"];pos=a["position"]
        m=lo.convert_coordinate("chr"+ch,pos-1) or []
        if len(m)!=1 or m[0][0]!="chr"+ch or m[0][2]!="+":
            out[rid]={"status":"HOLD_LIFTOVER","historical":a};continue
        out[rid]={"status":"MAPPED_PENDING_EAS_AIS","historical":a,"chr":ch,"pos":int(m[0][1])+1,
                  "ref":a["historical_ref_allele"],"alt":a["historical_alt_allele"]}
    require(len(out)==24,"Fixed panel changed")
    return {"status":"FIXED_24_LIFTOVER_PENDING_SOURCE_PARITY","historical_blob":"b34b4ae84ea4e5b0c72e3982d530a048cc9be5ca","chain_sha256":CHAIN_SHA,"rows":out,"withheld":old["withheld"],"missing":old["missing"]}
def open_source(study):
    acc=STUDIES[study];directory=base.BASE+acc+"/harmonised/"
    html,_=base.get(directory,2000000)
    names=sorted(set(re.findall(r'href=["\']([^"\']+\.h\.tsv\.gz)["\']',html.decode("utf8"))))
    require(len(names)==1 and re.fullmatch(re.escape(acc)+r"[^/]*\.h\.tsv\.gz",names[0]),"Exact harmonised file unresolved")
    url=urljoin(directory,names[0]);require(url.startswith(directory),"Directory escape")
    time.sleep(4)
    meta,_=base.get(url+"-meta.yaml",200000)
    md=yaml.safe_load(meta);gate_metadata(md,acc,study,names[0])
    time.sleep(4)
    idx,_=base.get(url+".tbi",8000000);tbi=base.parse_index(idx)
    time.sleep(4)
    first,cr=base.get(url,65536,(0,65535))
    total=int(cr.rsplit("/",1)[1])
    hdr=base.parse_bgzf_range(first,0,0)["data_lines"].readline().lstrip("#").rstrip("\r\n").split("\t")
    require(len(set(hdr))==len(hdr),"Duplicate source header")
    f=base.source_fields(hdr,tbi);require(f["rsid"] is not None,"Exact rsID field required")
    return (url,tbi,total,hdr,f),{"source":url,"metadata_sha256":base.sha(meta),"index_sha256":base.sha(idx),"header_range_sha256":base.sha(first),"publisher_full_file_md5":md.get("data_file_md5sum"),"trait_description":md["trait_description"],"samples":md["samples"],"source_total_bytes":total,"fields":f}
def one(src,r,rid,anchor):
    url,tbi,total,hdr,f=src;ch,pos=r["chr"],r["pos"]
    span=base.indexed_span(tbi,ch,pos,pos,total)
    require(span is not None,"No indexed chunk")
    a,b,x,y=span
    time.sleep(5) # Single sequential source requests; no automatic retries.
    raw,cr=base.get(url,b-a+1,(a,b))
    stream=base.parse_bgzf_range(raw,x,y)["data_lines"]
    matches=[]
    for row in csv.reader(stream,delimiter="\t"):
        if len(row)!=len(hdr):continue
        z=dict(zip(hdr,row))
        if z[f["chr"]].removeprefix("chr")==ch and z[f["pos"]]==str(pos) and {z[f["ea"]].upper(),z[f["oa"]].upper()}=={r["ref"],r["alt"]}:matches.append(z)
    require(len(matches)==1,"No unique exact allele-pair row")
    z=matches[0];rs=z[f["rsid"]].strip().lower()
    require(re.fullmatch(r"rs[0-9]+",rs) is not None,"No exact rsID")
    if not anchor:require(rs==r["rsid"],"rsID drift from EAS AIS anchor")
    value=base.orient(z,f,r,rs)
    if anchor:
        old=r["historical"]
        check=compare(old,[{"effect_allele":r["alt"],"other_allele":r["ref"],"beta":value["beta_alt"],"se":value["se"],"p":value["p"],"eaf":value["alt_eaf"]}])
        require(check["joint_effect_p_eaf_reproduced"],"Original archived EAS beta/p/EAF parity failed")
        value["historical_parity"]=check
    return {**value,"rsid":rs,"historical_ALT":r["alt"],"historical_REF":r["ref"],"GRCh38_chr":ch,"GRCh38_pos":pos,"source_range":cr,"source_range_sha256":base.sha(raw)}
def extract(study,manifest):
    require(len(manifest["rows"])==24,"Frozen panel size")
    anchor=study=="EAS_AIS"
    source_status="MAPPED_PENDING_EAS_AIS" if anchor else "ANCHOR_VERIFIED"
    out={"study":study,"accession":STUDIES[study],"status":"HOLD_SOURCE","rows":{}}
    try:
        src,meta=open_source(study);out.update(meta)
    except (OSError,ValueError,KeyError,TypeError,IndexError) as e:
        out["reason"]=str(e)[:240]
        out["rows"]={rid:{"status":"HOLD_SOURCE"} for rid in manifest["rows"]}
        return out
    for rid,r in sorted(manifest["rows"].items()):
        if r["status"]!=source_status:
            out["rows"][rid]={"status":"HOLD_ANCHOR","reason":r["status"]};continue
        try:out["rows"][rid]=one(src,r,rid,anchor)
        except (OSError,ValueError,KeyError,TypeError,IndexError) as e:
            out["rows"][rid]={"status":"HOLD_VARIANT","reason":str(e)[:180]}
    out["status"]="INDEXED_FIXED_PANEL_READBACK"
    out["matched"]=sum(v["status"]=="MATCHED_HARMONISED_ROW" for v in out["rows"].values())
    return out
def anchor_manifest(manifest,source):
    require(source["study"]=="EAS_AIS" and set(source["rows"])==set(manifest["rows"]),"Anchor identity")
    for rid,r in manifest["rows"].items():
        s=source["rows"][rid]
        if s["status"]=="MATCHED_HARMONISED_ROW":
            require(s["historical_parity"]["joint_effect_p_eaf_reproduced"] and r["status"]=="MAPPED_PENDING_EAS_AIS","Unapproved anchor")
            r.update(status="ANCHOR_VERIFIED",rsid=s["rsid"])
        else:r.update(status="HOLD_ANCHOR",source_status=s["status"])
    manifest["status"]="ANCHOR_EAS_AIS_REPRODUCTION"
    manifest["verified_count"]=sum(r["status"]=="ANCHOR_VERIFIED" for r in manifest["rows"].values())
    return manifest
def combine(manifest,root):
    require(len(manifest["rows"])==24,"Fixed panel size")
    d={}
    for k,acc in STUDIES.items():
        p=root/k/(k+".json")
        if p.exists():
            v=json.loads(p.read_text())
            require(v["study"]==k and v["accession"]==acc and set(v["rows"])==set(manifest["rows"]),"Source artifact identity")
            for rid,row in v["rows"].items():
                if row["status"]=="MATCHED_HARMONISED_ROW":
                    m=manifest["rows"][rid]
                    require(m["status"]=="ANCHOR_VERIFIED" and row["rsid"]==m["rsid"] and row["GRCh38_pos"]==m["pos"] and row["historical_ALT"]==m["alt"],"Matched row anchor drift")
            d[k]=v
        else:d[k]={"study":k,"accession":acc,"status":"HOLD_NO_ARTIFACT","rows":{rid:{"status":"HOLD_NO_ARTIFACT"} for rid in manifest["rows"]}}
    pairs=[("EAS_"+t,"EUR_"+t) for t in ("AIS","CES","LAS","SVS")]+[(a+"_SVS",a+"_"+t) for a in ("EUR","EAS") for t in ("CES","LAS")]
    rows=[]
    for rid,m in sorted(manifest["rows"].items()):
        for a,b in pairs:
            u,v=d[a]["rows"][rid],d[b]["rows"][rid]
            row={"region_id":rid,"rsid":m.get("rsid"),"first":a,"second":b,"status":"HOLD_UNMATCHED"}
            if u["status"]==v["status"]=="MATCHED_HARMONISED_ROW":
                c=base.contrast(u,v)
                c.pop("bonferroni_40_zero");c.pop("bonferroni_40_max")
                c.update(bonferroni_192_zero=min(1,FAMILY*c["p_zero_covariance"]),bonferroni_192_max=min(1,FAMILY*c["p_max_variance"]))
                row.update(c)
            rows.append(row)
    require(len(rows)==FAMILY,"Fixed test family")
    valid=[v for v in rows if v["status"]=="CONDITIONAL_REPORTED_BETA_CONTRAST"]
    return {"status":"EXPLORATORY_FIXED_24_SUBTYPE_COMPARISONS","family":FAMILY,"anchor_verified":manifest["verified_count"],"matched_cells":sum(sum(r["status"]=="MATCHED_HARMONISED_ROW" for r in v["rows"].values()) for v in d.values()),"matched_contrasts":len(valid),"significant_zero_covariance":[v for v in valid if v["bonferroni_192_zero"]<.05],"significant_max_variance":[v for v in valid if v["bonferroni_192_max"]<.05],"study_counts":{k:sum(r["status"]=="MATCHED_HARMONISED_ROW" for r in v["rows"].values()) for k,v in d.items()},"manifest":manifest,"studies":d,"contrasts":rows,"limits":["24 historically EAS-selected regional representatives, not LD-pruned independent loci. Only rsID/coordinate/allele-pair and EAS archive parity after forward liftover; not independent FASTA validation or whole-file checksum of harmonised sources.","Same GIGASTROKE release; prior 24 AIS and five-marker studies overlap. No claim of confirmatory or independent tests. Bonferroni 192 is within this extension, not project-wide.","Conditional numerical contrasts require comparable reported effect scales; shared cohort covariance unknown. Max variance bound does not fix selection, LD, power, effect-scale or ascertainment differences. Missing is not zero.","No gene causality, MR, coloc, molecular annotation or new discovery claim."]}
def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["anchor","extract","combine","selftest"]);p.add_argument("--legacy",type=Path);p.add_argument("--chain",type=Path);p.add_argument("--manifest",type=Path);p.add_argument("--study",choices=STUDIES);p.add_argument("--root",type=Path);p.add_argument("--out",type=Path);a=p.parse_args()
    if a.mode=="selftest":
        require(FAMILY==24*8 and len(STUDIES)==8,"Test family")
        good={"gwas_id":STUDIES["EAS_AIS"],"data_file_name":"test","genome_assembly":"GRCh38","is_harmonised":True,"trait_description":["Ischemic stroke"],"samples":[{"sample_ancestry_category":["East Asian"]}]}
        gate_metadata(good,STUDIES["EAS_AIS"],"EAS_AIS","test")
        try:gate_metadata(good,STUDIES["EAS_AIS"],"EUR_AIS","test")
        except ValueError:pass
        else:raise AssertionError("Incorrect ancestry accepted")
        print("Focused family/metadata assertions passed; no biological result.");return
    if a.mode=="anchor":
        m=panel(a.legacy,a.chain);s=extract("EAS_AIS",m)
        m=anchor_manifest(m,s);write(a.out/"manifest.json",m);write(a.out/"EAS_AIS.json",s)
        print(json.dumps({"verified_count":m["verified_count"],"source_status":s["status"],"rows":s["rows"]},indent=2,allow_nan=False))
    elif a.mode=="extract":
        x=extract(a.study,json.loads(a.manifest.read_text()));write(a.out,x)
        print(json.dumps({"study":x["study"],"status":x["status"],"matched":x.get("matched",0),"rows":x["rows"]},indent=2,allow_nan=False))
    else:
        x=combine(json.loads(a.manifest.read_text()),a.root);write(a.out,x)
        print(json.dumps({k:v for k,v in x.items() if k not in ("manifest","studies","contrasts")},indent=2,allow_nan=False))
if __name__=="__main__":main()
