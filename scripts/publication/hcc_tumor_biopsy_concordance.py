#!/usr/bin/env python3
"""Tumor-only public RNA/protein biopsy link and fixed patient-level concordance pilot."""
import argparse,collections,hashlib,io,json,math,re,statistics,sys,time
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import urlsplit
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"prior/scripts/publication"))
sys.path.insert(0,str(ROOT/"orthogonal/scripts/publication"))
import early_track_orthogonal_sources as orth
import early_track_orthogonal_repair as schema
base=orth.base;require=base.require
ARCH=ROOT/"docs/publication/evidence/hcc_actual_schema_readback_20261125.json"
SUP="https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12943-026-02604-x/MediaObjects/12943_2026_2604_MOESM3_ESM.xlsx"
def prior():
    q=json.loads(ARCH.read_text());require(q["run_id"]==45612527929 and q["job_id"]==139199468686,"prior schema provenance");return q["aggregate"]
def hgnc():
    result={}
    for g in orth.PANEL:
        url="https://rest.genenames.org/fetch/symbol/"+g
        z={"status":"HOLD_GENE_MAPPING"}
        try:
            with urlopen(Request(url,headers={"Accept":"application/json","User-Agent":"HCC-fixed-panel/1.0"}),timeout=30) as f:
                require(f.status==200 and urlsplit(f.geturl()).hostname=="rest.genenames.org","HGNC source")
                raw=f.read(300001)
            require(0<len(raw)<=300000,"HGNC byte cap")
            q=json.loads(raw)["response"];d=q["docs"]
            require(q["numFound"]==1 and len(d)==1 and d[0]["symbol"]==g and d[0]["status"]=="Approved","not exact approved symbol")
            ens=d[0].get("ensembl_gene_id");require(isinstance(ens,str) and re.fullmatch(r"ENSG\d+",ens),"no unique Ensembl ID")
            z.update(status="EXACT_HGNC_ENSEMBL",ensembl=ens,hgnc=d[0]["hgnc_id"],source=url,sha256=hashlib.sha256(raw).hexdigest())
        except (OSError,ValueError,KeyError,TypeError,IndexError) as ex:z["reason"]=str(ex)[:150]
        result[g]=z;time.sleep(.5)
    return result
def manifest():
    import openpyxl
    raw,src=base.get(SUP,15_000_000)
    require(src["sha256"]==prior()["manifest"]["source"]["sha256"],"manifest SHA drift")
    w=openpyxl.load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
    try:
        ws=w.worksheets[0];require((ws.max_row,ws.max_column)==(78,32),"published metadata shape")
        rr=ws.iter_rows(min_row=2,values_only=True);h=next(rr)
        require(h[:4]==("Patient Number","Biopsy ID","Tumor Nodule Number (radiological assignment)","Tumor Label"),"exact tumor-biopsy header")
        out={}
        for v in rr:
            p,k=schema.key(v[0]),schema.key(v[1])
            require(p is not None and k is not None and k not in out and v[2] is not None and v[3] is not None,"unresolved or duplicate tumor biopsy")
            out[k]=p
        require(len(out)==76 and len(set(out.values()))==22,"reported 76 tumor biopsies / 22 patients")
        return out,src
    finally:w.close()
def cols(ws):
    h=next(ws.iter_rows(min_row=1,max_row=1,values_only=True));require(h[0] in (None,""),"unlabeled first-column contract changed")
    x=[schema.key(v) for v in h[1:]]
    require(all(v is not None for v in x) and len(set(x))==len(x),"nonunique assay keys")
    return {v:i+1 for i,v in enumerate(x)}
def vectors(ws,keys,c,lookup,is_rna):
    out=collections.defaultdict(list)
    for row in ws.iter_rows(min_row=2,values_only=True):
        if not isinstance(row[0],str):continue
        g=row[0].split(".")[0] if is_rna else row[0]
        if g not in lookup:continue
        v=[]
        for k in keys:
            a=row[c[k]]
            v.append(float(a) if isinstance(a,(int,float)) and not isinstance(a,bool) and math.isfinite(a) else None)
        out[lookup[g]].append(v)
    return out
def rho(x,y):
    from scipy.stats import spearmanr
    if len(x)<5 or len(set(x))<2 or len(set(y))<2:return None
    r=float(spearmanr(x,y).statistic)
    return r if math.isfinite(r) else None
def analyze():
    import openpyxl
    m,ms=manifest();mp=hgnc();q=prior();books={};sources={}
    for name in orth.ZEN:
        raw,src=orth.zen_download(name)
        require(src["sha256"]==q["workbooks"][name]["source"]["sha256"],"prior workbook SHA drift")
        sources[name]=src;books[name]=openpyxl.load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
    rn,pn=list(orth.ZEN)
    try:
        # Exact intersection with published tumor biopsy manifest, not sheet-order tissue assignment.
        candidates=[]
        for i in (2,3):
            for j in (2,3):
                a,b=cols(books[rn].worksheets[i-1]),cols(books[pn].worksheets[j-1])
                common=set(a)&set(b)&set(m)
                if common:candidates.append((i,j,a,b,common))
        out={"status":"HCC_MANIFEST_TUMOR_ONLY_LINK","sources":sources,"manifest_source":ms,"mapping":mp,
             "candidate_intersections":[{"rna_sheet":i,"protein_sheet":j,"manifest_tumor_keys":len(k)} for i,j,a,b,k in candidates]}
        exact=[v for v in candidates if len(v[4])==22 and len({m[k] for k in v[4]})==9]
        if len(exact)!=1:
            out["hold_reason"]="No unique 22-biopsy / 9-patient cross-assay tumor intersection";return out
        i,j,a,b,common=exact[0];keys=sorted(common)
        require(not any(common&v[4] for v in candidates if (v[0],v[1])!=(i,j)),"ambiguous cross-sheet manifest match")
        rlookup={v["ensembl"]:g for g,v in mp.items() if v["status"]=="EXACT_HGNC_ENSEMBL"};require(len(set(rlookup.values()))==len(rlookup),"ambiguous gene mapping")
        rv=vectors(books[rn].worksheets[i-1],keys,a,rlookup,True)
        pv=vectors(books[pn].worksheets[j-1],keys,b,{g:g for g in orth.PANEL},False)
        results={}
        for g in orth.PANEL:
            z={"RNA_rows":len(rv[g]),"protein_rows":len(pv[g]),"status":"HOLD_MISSING_OR_NONUNIQUE_FEATURE"}
            if len(rv[g])==len(pv[g])==1:
                groups=collections.defaultdict(list)
                for k,x,y in zip(keys,rv[g][0],pv[g][0]):
                    if x is not None and y is not None:groups[m[k]].append((x,y))
                p=[(statistics.median(v[0] for v in rr),statistics.median(v[1] for v in rr)) for rr in groups.values()]
                x=[v[0] for v in p];y=[v[1] for v in p];r=rho(x,y)
                z.update(status="HOLD_NONCONSTANT_OR_INSUFFICIENT_PATIENTS",complete_biopsies=sum(len(v) for v in groups.values()),patients=len(p))
                if r is not None:
                    loo=[rho(x[:k]+x[k+1:],y[:k]+y[k+1:]) for k in range(len(p))]
                    loo=[v for v in loo if v is not None]
                    z.update(status="DESCRIPTIVE_TUMOR_PATIENT_MEDIAN_CONCORDANCE",spearman_rho=r,loo_rho_range=[min(loo),max(loo)] if loo else None)
            results[g]=z
        out.update(matched_tumor_biopsies=22,patients=9,results=results)
        out["limits"]=["Multifocal HCC tumor only; 22 exact biopsies of nine patients, within-patient medians of published assay values. Spearman across patients, no raw scale subtraction, p-values, clinical outcome or independent validation. Normal/NTL columns withheld because manifest provenance unavailable. Existing Ng 2022 cohort reused. No patient or specimen keys emitted."]
        return out
    finally:
        for w in books.values():w.close()
def selftest():
    assert len(orth.PANEL)==8
    assert rho([1,2,3,4,5],[5,4,3,2,1]) < -.99
    assert rho([1,2,3,4],[1,2,3,4]) is None
    assert prior()["workbooks"][list(orth.ZEN)[0]]["sheets"][1]["first_column_shapes"]=={"ensembl_gene":16483}
    print("Fixed gene panel, archived source schema, and patient-level rank checks passed.")
def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["selftest","analyze"]);p.add_argument("--out",type=Path);a=p.parse_args()
    if a.mode=="selftest":selftest();return
    try:z=analyze()
    except (OSError,ValueError,KeyError,TypeError,IndexError,EOFError) as e:z={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:200]}
    z["scientific_status"]="HOLD_NOVELTY_AND_INDEPENDENT_VALIDATION"
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n");print(json.dumps(z,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
