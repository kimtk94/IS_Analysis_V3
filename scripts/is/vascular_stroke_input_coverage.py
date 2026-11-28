#!/usr/bin/env python3
"""Resolve published stroke tag/eQTL inclusion, separately from reported coloc."""
import argparse, collections, hashlib, json, re
from pathlib import Path
from openpyxl import load_workbook
from vascular_published_candidate_audit import TARGETS, tokens
PIN="8a7a0cdd8d4d80bca570b80eeace2e3aa65d1a051239a9a08cbcbbf7dc823e14"
def require(ok,msg):
    if not ok:raise ValueError(msg)
def sha(p):
    with p.open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()
def norm(x):return re.sub(r"\s+"," ",str(x or "").strip()).casefold()
def cells(row):return [str(v) if v is not None else "" for v in row]
def read_sheet(w,n):
    s=w[f"Supplementary Data {n}"]
    require(s.max_column<=300 and s.max_row<=100000,"Sheet bound")
    return [tuple(r) for r in s.iter_rows(values_only=True)]
def exact(v):
    if v is None:return ""
    t=str(v).strip().upper()
    return t if t in TARGETS else ""
def eqtl(rows):
    require("stroke" in norm(rows[0][0]) and "eqtl" in norm(rows[0][0]),"Wrong stroke eQTL sheet title")
    choices=[]
    for i,row in enumerate(rows[:12]):
        h=[norm(v) for v in row]
        if all(x in h for x in ("snp id","proxy snp","gene name","in this study","eqtl")):choices.append((i,h))
    require(len(choices)==1,"Unknown eQTL header")
    i,h=choices[0]
    starts=[j for j,x in enumerate(h) if x=="line number"]
    require(len(starts)==3,"Expected three pooled/sex-specific groups")
    out=[]
    for k,j in enumerate(starts):
        end=starts[k+1] if k+1<len(starts) else len(h)
        sub=h[j:end]
        fields=("snp id","proxy snp","in this study","eqtl","gene name","p-value","beta")
        require(all(sub.count(f)==1 for f in fields),"Unexpected eQTL group schema")
        header=rows[i-1][j:end]
        group=next((str(v).strip() for v in header if v),"")
        require("vsmc" in group.casefold(),"Missing VSMC group label")
        counts=collections.Counter(); matches=[]; tags=set()
        for ri,row in enumerate(rows[i+1:],i+2):
            rec={f:row[j+sub.index(f)] for f in fields}
            tag=str(rec["snp id"] or "").strip().upper()
            if not re.fullmatch(r"RS\d+",tag):continue
            tags.add(tag);counts["rows"]+=1
            if norm(rec["in this study"])=="yes":counts["in_study_rows"]+=1
            if norm(rec["eqtl"])=="yes":counts["eqtl_yes_rows"]+=1
            fields_matched={f:exact(rec[f]) for f in ("snp id","proxy snp","gene name") if exact(rec[f])}
            if fields_matched:
                matches.append({"row":ri,"exact_fields":fields_matched,"source_values":{f:str(v or "") for f,v in rec.items()}})
        out.append({"group":group,"unique_tags":len(tags),"counts":dict(counts),"fixed_tag_ids":sorted(tags&TARGETS),"matches":matches})
    return {"title":str(rows[0][0]),"groups":out}
def taglist(rows):
    title=str(rows[0][0])
    # Table titles and actual header names must be inspected rather than assuming rsID column.
    require("stroke" in title.casefold(),"Wrong original stroke tag sheet")
    candidates=[]
    for i,row in enumerate(rows[:15]):
        for j,v in enumerate(row):
            if norm(v) in {"snp id","rsid","rs id","rs number","snp","variant id","tagging snp","tag snp"}:
                ids=[str(x[j] or "").strip().upper() for x in rows[i+1:]]
                ids=[v for v in ids if re.fullmatch(r"RS\d+",v)]
                if len(set(ids))>=100:candidates.append((i,j,ids))
    if len(candidates)!=1:
        return {"status":"HOLD_TAG_COLUMN","title":title,"headers":[cells(r)[:40] for r in rows[:8]],"candidate_columns":[{"row":i+1,"column":j+1,"unique_rsids":len(set(ids))} for i,j,ids in candidates]}
    i,j,ids=candidates[0]
    return {"status":"EXACT_TAG_COLUMN" if len(set(ids))==139 else "HOLD_TAG_DENOMINATOR","title":title,"header":cells(rows[i]),"rsid_column":j+1,"rsid_rows":len(ids),"unique_rsids":len(set(ids)),"fixed_tag_ids":sorted(set(ids)&TARGETS),"fixed_exact_rows":[{"row":k,"values":cells(row)} for k,row in enumerate(rows[i+1:],i+2) if exact(row[j])]}
def tissue(rows):
    require("colocal" in norm(rows[0][0]),"Unexpected tissue sheet")
    h=[norm(x) for x in rows[3]];require(h[:6]==["tissue","gene name","cad","bp/hypertension","stroke","aaa"],"Tissue trait header changed")
    out=[]
    for i,r in enumerate(rows[4:],5):
        if exact(r[1]):
            out.append({"row":i,"tissue":str(r[0]),"gene":str(r[1]),"CAD":str(r[2] or ""),"BP":str(r[3] or ""),"Stroke":str(r[4] or ""),"AAA":str(r[5] or ""),"stroke_statistics":dict(zip(cells(rows[3][22:30]),cells(r[22:30])))})
    return out
def run(a,b):
    require(sha(a)==PIN,"Changed first publisher workbook")
    w=load_workbook(a,read_only=True,data_only=True,keep_links=False)
    v=load_workbook(b,read_only=True,data_only=True,keep_links=False)
    out={"workbook_1_sha256":sha(a),"workbook_2_sha256":sha(b),"fixed_targets":sorted(TARGETS),"status":"PUBLISHED_INPUT_COVERAGE_ONLY","limitations":["Original tag list and reported eQTL rows are not complete regional statistics.","No new association, LD, colocalization or ancestry comparison. No absence-of-effect claim."]}
    try:
        for name,fn in [("stroke_eqtl",lambda:eqtl(read_sheet(w,2))),("stroke_taglist",lambda:taglist(read_sheet(v,31))),("vascular_tissue_traits",lambda:tissue(read_sheet(w,13)))]:
            try:out[name]=fn()
            except (ValueError,KeyError,IndexError) as e:out[name]={"status":"HOLD_SCHEMA","reason":str(e)[:240]}
    finally:w.close();v.close()
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--first",type=Path);p.add_argument("--second",type=Path);p.add_argument("--out",type=Path);p.add_argument("--selftest",action="store_true");a=p.parse_args()
    if a.selftest:
        assert exact(" rs9515201 ")=="RS9515201" and exact("rs95152010")==""
        assert exact("COL4A2;COL4A1")==""
        assert norm("In  this study")=="in this study"
        print("Focused schema/token assertions passed; no scientific claim.");return
    x=run(a.first,a.second);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps(x,indent=2,allow_nan=False))
if __name__=="__main__":main()
