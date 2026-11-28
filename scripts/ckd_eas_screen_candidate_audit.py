#!/usr/bin/env python3
"""Verify derived assay-screen results and inspect published reagent evidence, without rerunning GWAS downloads."""
import argparse
import collections
import hashlib
import json
import math
import zipfile
from pathlib import Path
from statistics import NormalDist

PINS = {"egfr":"ad0d9181d3d1b4d598b7c9e5a3d5be11825db8681a877069a4410c96834328c5",
        "bun":"a56ee97acc8ac60231a5ae68c39bcc8357c06651805721fed3302f4454a1cdf2"}
MAPPING_SHA = "b1180c9f633b378b9b9793b3093c702aa393172bcabe28e83fcc3755eaff4bd0"
def digest(p):
    with p.open("rb") as f:
        return hashlib.file_digest(f,"sha256").hexdigest()

def close(a,b):
    return math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-14)

def read_archive(root,trait):
    p=root/(trait+".zip")
    if digest(p)!=PINS[trait]:
        raise ValueError("Artifact digest mismatch: "+trait)
    with zipfile.ZipFile(p) as z:
        if z.testzip() is not None:
            raise ValueError("Archive CRC")
        def read(name):
            m=[n for n in z.namelist() if n.rsplit("/",1)[-1]==name]
            if len(m)!=1:
                raise ValueError("Ambiguous artifact member: "+name)
            return json.loads(z.read(m[0]))
        s,rows=read("summary.json"),read("assay_results.json")
    if s["trait"]!=trait or s["source_rows"]!=4512904:
        raise ValueError("Unexpected source")
    d={(r["platform"],r["assay_id"]):r for r in rows}
    if len(d)!=len(rows):
        raise ValueError("Duplicate assay")
    replay=0
    for platform in ("Olink","SomaScan"):
        a=[r for r in rows if r["platform"]==platform]
        family=len(a)
        if family!=s["platforms"][platform]["family"]:
            raise ValueError("Family discrepancy")
        valid=[r for r in a if r["status"]=="EXPLORATORY_WALD"]
        if len(valid)!=s["platforms"][platform]["matched"]:
            raise ValueError("Match count mismatch")
        for r in valid:
            bx,sx,by,sy=(r[k] for k in ("beta","se","outcome_beta_aligned","outcome_se"))
            v=by/bx
            se=math.hypot(sy, v*sx)/abs(bx)
            # Use erfc rather than subtracting a near-one normal CDF.
            pvalue=math.erfc(abs(v/se)/math.sqrt(2))
            if not all(close(x,y) for x,y in ((v,r["mr_beta"]),(se,r["mr_se_delta_independence"]),
                      (pvalue,r["mr_p"]),(min(1,pvalue*family),r["bonferroni_p"]))):
                raise ValueError("Wald or Bonferroni replay mismatch")
            replay+=1
        ordered=sorted(a,key=lambda r:r.get("mr_p",1.0))
        m=1.0
        for i,r in reversed(list(enumerate(ordered,1))):
            m=min(m,r.get("mr_p",1.0)*family/i)
            if not close(m,r["bh_q_conservative_full_family"]):
                raise ValueError("BH replay mismatch")
    return s,d,replay

def mapping(root, keys):
    from openpyxl import load_workbook
    p=root/"mapping.xlsx"
    if digest(p)!=MAPPING_SHA:
        raise ValueError("Mapping digest mismatch")
    w=load_workbook(p,read_only=True,data_only=True)
    it=w["Supplementary_Data_1"].iter_rows(values_only=True)
    h=[str(x) for x in next(it)]
    required={"uniprot_id","olink_id","somascan_id","rho_olink_soma_non_ANML"}
    if not required<=set(h):
        raise ValueError("Mapping schema")
    matches=collections.defaultdict(list)
    for row in it:
        r=dict(zip(h,row))
        for platform,field in (("Olink","olink_id"),("SomaScan","somascan_id")):
            key=(platform,str(r[field]))
            if key in keys:
                matches[key].append({k:r[k] for k in sorted(required)})
    w.close()
    return matches

def published_inventory(root,keys):
    from openpyxl import load_workbook
    path=root/"published.xlsx"
    w=load_workbook(path,read_only=True,data_only=True)
    sheets=[]
    targets={a for _,a in keys}
    for s in w.worksheets:
        it=s.iter_rows(values_only=True)
        header=[str(x or "") for x in next(it)]
        # Never infer evidence from substrings; identify exact candidate reagent labels in columns.
        matched=[]
        count=0
        for row in it:
            count+=1
            found=sorted({str(x).strip() for x in row if x is not None and str(x).strip() in targets})
            if found:
                # Report only candidate rows, capped per sheet, exact source header strings.
                matched.append({"matched_labels":found,"fields":{k:v for k,v in zip(header,row) if k and v is not None}})
        sheets.append({"sheet":s.title,"header":header,"rows_scanned":count,"matched_rows":len(matched),
                       "candidate_rows":matched[:100],"truncated":len(matched)>100})
    w.close()
    return {"sha256":digest(path),"sheets":sheets,"status":"EXACT_LABEL_INVENTORY_NOT_FORMAL_COLOC"}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    e,er,ne=read_archive(a.root,"egfr")
    b,br,nb=read_archive(a.root,"bun")
    if e["anchor_selection_sha256"]!=b["anchor_selection_sha256"] or set(er)!=set(br):
        raise ValueError("Outcomes have different frozen anchor sets")
    for k,r in er.items():
        t=br[k]
        for field in ("sentinel","beta","se","F"):
            if r.get(field)!=t.get(field):
                raise ValueError("Exposure differs across endpoints")
    keys={k for k in er if er[k].get("bonferroni_p",1)<.05 or br[k].get("bonferroni_p",1)<.05}
    mp=mapping(a.root,keys)
    pub=published_inventory(a.root,keys)
    out=[]
    for platform,assay in sorted(keys):
        key=(platform,assay)
        x,y=er[key],br[key]
        row={"platform":platform,"assay_id":assay,"reagent_pairs":mp[key],
             "sentinel":x.get("sentinel"),"rsid":x.get("rsid"),"F":x.get("F")}
        for trait,r in (("egfr",x),("bun",y)):
            row[trait]={k:r.get(k) for k in ("status","mr_beta","mr_ci95","mr_p","bonferroni_p","bh_q_conservative_full_family")}
        row["both_bonferroni"]=x.get("bonferroni_p",1)<.05 and y.get("bonferroni_p",1)<.05
        row["opposite_direction"]=(x.get("mr_beta",0)*y.get("mr_beta",0)<0)
        out.append(row)
    a.out.mkdir(parents=True,exist_ok=True)
    result={"status":"DERIVED_NUMERIC_REPLAY_AND_PUBLISHED_MAPPING","anchor_sha256":e["anchor_selection_sha256"],
            "wald_rows_replayed":{"egfr":ne,"bun":nb},"candidate_union":len(out),
            "both_significant_opposite_direction":[[r["platform"],r["assay_id"]] for r in out if r["both_bonferroni"] and r["opposite_direction"]],
            "candidate_rows":out,"published_inventory_headers":[{k:v for k,v in s.items() if k!="candidate_rows"} for s in pub["sheets"]],
            "publication_workbook_sha256":pub["sha256"],
            "limitation":"Derived numeric replay is not independent source extraction; exact label rows require schema interpretation."}
    (a.out/"candidate_audit.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    (a.out/"published_label_rows.json").write_text(json.dumps(pub,indent=2,allow_nan=False,default=str)+"\n")
    print(json.dumps(result,indent=2,allow_nan=False,default=str))
    # Log only narrow candidate matched source rows, not whole publisher workbook.
    for s in pub["sheets"]:
        print("PUBLISHED_SOURCE",s["sheet"],json.dumps(s["candidate_rows"][:20],allow_nan=False,default=str))
if __name__=="__main__":
    main()
