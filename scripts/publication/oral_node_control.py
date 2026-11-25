#!/usr/bin/env python3
"""Exploratory same-patient node controls, exact author labels, no individual records exported."""
import argparse,collections,csv,gzip,io,json,math,re,time
from pathlib import Path
import early_track_orthogonal_repair as old
base=old.base;require=base.require;ROOT=old.ROOT
ARCH=ROOT/"docs/publication/evidence/oral_exact_label_readback_20261125.json"
LIBS={"4":[4],"6":[6,7],"7":[8,9]}
SITES=("LN_MET","LN_FREE")
GENES=("CXCR4","SPP1")
MIN_CELLS=50
def evidence():
    a=json.loads(ARCH.read_text())
    require(a["run_id"]==45612527929 and a["job_id"]==139199468685,"archive provenance")
    rows=a["aggregate"]["records"];require(len(rows)==9,"nine libraries required")
    require(sum(v["columns"][0]["exact_barcode_matches"] for v in rows)==41284,"frozen annotation denominator")
    return {v["gsm"]:v for v in rows}
def classify(rr,bars,lib):
    require(rr and all(len(v)==3 for v in rr),"three-column source contract")
    h=rr[0];require(h[:2]==["cell","celltype_classification"],"exact source header")
    # Only author-defined sample/HTO fields and values; no core barcode guesses.
    require(old.old.norm(h[2]) in {"sample","sampleid","hashid","hashidman","hashglobal","hto","htoid","htoclassification","htoclassificationglobal","patient","sampleclassification","sampleclassificationman","sampleorigin","hash","hashing","hashclassification"},"unrecognized public sample-code header: "+str(h[2])[:60])
    index={v:i for i,v in enumerate(bars)};require(len(index)==len(bars),"duplicate raw barcodes")
    rows=rr[1:];require(rows and len({v[0] for v in rows})==len(rows) and all(v[0] in index for v in rows),"not exact complete barcode join")
    masks={s:[] for s in SITES};cnt=collections.Counter();unknown=0
    for barcode,label,hashcode in rows:
        site=old.SITE_CODES[lib].get(hashcode)
        if site is None:
            if re.fullmatch(r"(?:(?:[1-7][ab]?)\.)?(?:Negative|Doublet)",hashcode):site="EXCLUDED_HASH"
            else:unknown+=1;site="UNKNOWN"
        cnt[site]+=1
        if label=="Macrophages" and site in SITES:masks[site].append(index[barcode])
    require(unknown==0,"unrecognized author hash values; no site guessing")
    return masks,dict(cnt),len(rows),sum(v[1]=="Macrophages" for v in rows)
def metadata(lib,prior):
    sid=old.old.GSMS[lib-1];p=prior[sid];folder,names,listing=old.old.filelist(sid);data={};src={}
    for k,pat,cap in (("barcode",r"barcodes\.tsv\.gz$",15_000_000),("label",r"celltype_metadata.*\.csv\.gz$",5_000_000)):
        data[k],src[k]=base.get(folder+old.old.one(names,pat),cap)
        require(src[k]["sha256"]==p["sources"][k]["sha256"],"prior source SHA drift")
    bars=[v[0] for v in csv.reader(io.StringIO(old.old.text_gz(data["barcode"])),delimiter="\t") if v]
    rr=list(csv.reader(io.StringIO(old.old.text_gz(data["label"]))))
    masks,sites,n,macro=classify(rr,bars,lib)
    require(n==p["columns"][0]["exact_barcode_matches"] and macro==p["columns"][1]["recognized_cell_labels"].get("Macrophages",0),"prior exact annotation count drift")
    return {"lib":lib,"sid":sid,"folder":folder,"files":names,"bars":bars,"masks":masks,
            "aggregate":{"gsm":sid,"matched_labels":n,"macrophages":macro,"site_label_counts":sites,"site_macro_counts":{s:len(masks[s]) for s in SITES},"sources":src}}
def quantify(m):
    import numpy as np
    from scipy.io import mmread
    src={};d={}
    for k,p,cap in (("features",r"features\.tsv\.gz$",5_000_000),("matrix",r"matrix\.mtx\.gz$",180_000_000)):
        d[k],src[k]=base.get(m["folder"]+old.old.one(m["files"],p),cap)
    f=[v for v in csv.reader(io.StringIO(old.old.text_gz(d["features"])),delimiter="\t") if v]
    require(f and all(len(v)==3 for v in f) and len({v[0] for v in f})==len(f),"unique typed feature rows required")
    ge=[i for i,v in enumerate(f) if v[2]=="Gene Expression"]
    ix={g:[i for i in ge if f[i][1]==g] for g in GENES}
    require(len(ge)>=10000 and all(len(v)==1 for v in ix.values()),"unique RNA genes or denominator missing")
    mat=mmread(gzip.GzipFile(fileobj=io.BytesIO(d["matrix"])))
    require(mat.shape==(len(f),len(m["bars"])) and np.issubdtype(mat.data.dtype,np.integer) and (mat.data>0).all(),"full integer sparse matrix contract")
    nnz=mat.nnz;mat=mat.tocsc();require(mat.nnz==nnz,"duplicate coordinates")
    out={}
    for site in SITES:
        cols=m["masks"][site]
        if not cols:out[site]={"cells":0,"umis":0,"gene_counts":{g:0 for g in GENES}};continue
        z=mat[:,cols]
        totals=np.asarray(z[ge,:].sum(axis=0)).ravel()
        require((totals>0).all(),"labeled macrophage without RNA counts")
        den=int(totals.sum())
        out[site]={"cells":len(cols),"umis":den,"gene_counts":{g:int(z[ix[g][0],:].sum()) for g in GENES}}
    return out,{"source":src,"matrix_shape":[len(f),len(m["bars"])],"nonzero":int(nnz),"RNA_features":len(ge)}
def patient(p):
    prior=evidence();result={"status":"HOLD_SOURCE_OR_ELIGIBILITY","study_patient_code":p,"libraries":[]}
    ms=[]
    for lib in LIBS[p]:
        try:m=metadata(lib,prior);ms.append(m);result["libraries"].append(m["aggregate"])
        except (OSError,EOFError,ValueError,KeyError,IndexError,TypeError,UnicodeError) as e:result["libraries"].append({"library":lib,"reason":str(e)[:180]})
    if len(ms)!=len(LIBS[p]):return result
    counts={s:sum(len(m["masks"][s]) for m in ms) for s in SITES};result["site_macro_counts"]=counts
    if min(counts.values())<MIN_CELLS:
        result["status"]="HOLD_FEWER_THAN_50_MACROPHAGES_PER_SITE";return result
    raw=[]
    for m in ms:
        try:q,src=quantify(m);raw.append(q);result["libraries"][ms.index(m)]["matrix_verification"]=src
        except (OSError,EOFError,ValueError,KeyError,IndexError,TypeError,UnicodeError) as e:result["matrix_hold"]=str(e)[:180];return result
    sites={}
    for site in SITES:
        n=sum(x[site]["cells"] for x in raw);den=sum(x[site]["umis"] for x in raw)
        require(n==counts[site] and den>0,"aggregate denominator")
        sites[site]={"cells":n,"CPM":{g:sum(x[site]["gene_counts"][g] for x in raw)*1e6/den for g in GENES}}
    # Per-patient data retained ONLY as study-coded aggregate sign/delta; no barcode or cell data.
    result.update(status="EXPLORATORY_PAIRED_NODE_CONTEXT",results={g:math.log2(sites["LN_MET"]["CPM"][g]+1)-math.log2(sites["LN_FREE"]["CPM"][g]+1) for g in GENES})
    return result
def combine(root):
    d={}
    for p in LIBS:
        f=root/(p+".json");z=json.loads(f.read_text()) if f.exists() else {}
        require(not z or z.get("study_patient_code")==p,"artifact identity")
        d[p]=z
    e=[v for v in d.values() if v.get("status")=="EXPLORATORY_PAIRED_NODE_CONTEXT"]
    import statistics
    out={}
    for g in GENES:
        v=[z["results"][g] for z in e]
        out[g]={"eligible_pairs":len(v),"positive":sum(x>0 for x in v),"negative":sum(x<0 for x in v),"ties":sum(x==0 for x in v),"median_delta":statistics.median(v) if v else None}
    return {"status":"ORAL_NODE_CONTROL_DESCRIPTIVE_ONLY","candidate_patients":3,"eligible_pairs":len(e),"patient_gate_status":{p:z.get("status","HOLD_MISSING") for p,z in d.items()},"results":out,
            "limits":["At most three author-coded patients, patient 6 has 59 macrophages total before site assignment and cannot pass both 50-cell gates. No p-values with this exploratory sample size, no threshold relaxation, no new mechanistic claim. Same author labels, exact HTO codes; combine within patient/site, not cells as independent observations."]}
def selftest():
    evidence()
    b=["AAAAAAAAAAAAAAAA-1","CCCCCCCCCCCCCCCC-1"]
    rr=[["cell","celltype_classification","hash.ID"],[b[0],"Macrophages","4.1"],[b[1],"Macrophages","4.3"]]
    m,c,n,ma=classify(rr,b,4);assert m=={"LN_MET":[0],"LN_FREE":[1]} and n==ma==2
    rr[1][2]="7.2"
    try:classify(rr,b,4)
    except ValueError:pass
    else:raise AssertionError("foreign sample code accepted")
    assert sum(len(v) for v in LIBS.values())==5
    print("Prior exact barcode archive, exact author code, wrong-code rejection, and fixed patient set passed.")
def main():
    ap=argparse.ArgumentParser();ap.add_argument("mode",choices=[*LIBS,"combine","selftest"]);ap.add_argument("--out",type=Path);ap.add_argument("--root",type=Path,default=Path("inputs"));a=ap.parse_args()
    if a.mode=="selftest":selftest();return
    try:z=combine(a.root) if a.mode=="combine" else patient(a.mode)
    except (OSError,EOFError,ValueError,KeyError,IndexError,TypeError,UnicodeError) as e:z={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:180],"study_patient_code":a.mode}
    z["scientific_status"]="HOLD_NOVELTY_AND_INDEPENDENT_VALIDATION"
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n");print(json.dumps(z,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
