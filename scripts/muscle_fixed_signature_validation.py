#!/usr/bin/env python3
"""Frozen acute RNA signature portability. Only aggregate public study summaries exported."""
import argparse, collections, csv, gzip, hashlib, io, json, math, re, statistics, urllib.request
from pathlib import Path
import numpy as np
from scipy.stats import rankdata, t as student_t, binomtest

BASE="https://ftp.ncbi.nlm.nih.gov/geo"
META_42507="4894975ae096b7a9da66d1123344df14b84fa9060c6c066604482acbab4f0ad9"
META_277819="12ddfca1171d9fde5fd5f07219a9a6dfed8b16771b9ba677474451527dee6eb3"
COUNTS_277819="7b928762b2b5eeafc19729168df5ae5b041ea974f5634e0b430eef1f59ea4007"
def require(ok,msg):
    if not ok:raise ValueError(msg)
def download(url,path,cap,sha=None):
    req=urllib.request.Request(url,headers={"User-Agent":"Public-muscle-signature/1.0"})
    with urllib.request.urlopen(req,timeout=45) as r,path.open("wb") as f:
        h=hashlib.sha256();n=0
        while b:=r.read(1024*1024):
            n+=len(b);require(n<=cap,"Source byte limit")
            h.update(b);f.write(b)
    require(n>0 and (sha is None or h.hexdigest()==sha),"Source identity mismatch")
    return {"url":url,"bytes":n,"sha256":h.hexdigest()}
def series(acc):
    return BASE+"/series/"+acc[:-3]+"nnn/"+acc
def read_panel(path,key):
    with path.open() as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    require(40<=len(rows)<=200,"Panel count")
    p={x[key]:int(x["sign"]) for x in rows}
    require(len(p)==len(rows) and set(p.values())=={-1,1},"Frozen panel identity")
    require(all(20<=sum(v==s for v in p.values())<=100 for s in (-1,1)),"Sign family")
    require(all(float(x["adj_p_value"])<.05 and int(x["sign"])*float(x["logFC"])>0 for x in rows),"Source panel gate")
    return p
def matrix(path,acc,n,values):
    meta=collections.defaultdict(list);series_fields={};rows={};header=None;table=False
    with gzip.open(path,"rt",encoding="utf-8-sig",newline="") as f:
        for line in f:
            line=line.rstrip("\r\n")
            if line=="!series_matrix_table_begin":table=True;continue
            if line=="!series_matrix_table_end":table=False;break
            if line.startswith("!Sample_"):
                a=next(csv.reader([line],delimiter="\t"));meta[a[0]].append(a[1:])
            elif line.startswith("!Series_geo_accession"):
                a=next(csv.reader([line],delimiter="\t"));series_fields[a[0]]=a[1:]
            elif table:
                a=next(csv.reader([line],delimiter="\t"))
                if header is None:header=a[1:];continue
                if not values:raise ValueError("Unexpected expression in metadata-only source")
                require(len(a)==n+1 and a[0] not in rows and len(rows)<100000,"Expression width/feature identity")
                try:v=np.array([float(x) for x in a[1:]],dtype=float)
                except ValueError:continue
                if np.isfinite(v).all():rows[a[0]]=v
    require(series_fields.get("!Series_geo_accession")==[acc],"GEO accession")
    ids=meta["!Sample_geo_accession"][0]
    require(len(ids)==len(set(ids))==n and header==ids,"Exact sample order")
    require(all(len(a)==n for m in meta.values() for a in m),"Metadata width")
    chars=collections.defaultdict(dict)
    for m in meta["!Sample_characteristics_ch1"]:
        for i,v in enumerate(m):
            if ": " in v:
                k,val=v.split(": ",1);k=k.strip().casefold()
                require(i not in chars[k],"Repeated characteristic")
                chars[k][i]=val.strip()
    return rows,meta,chars
def platform(path):
    # Only original documented source ID and unique Entrez ID. No gene-symbol guess.
    opener=gzip.open if path.suffix==".gz" else open
    with opener(path,"rt",encoding="utf-8-sig",errors="strict",newline="") as f:
        found=False;header=None;rows=collections.defaultdict(set);n=0
        for line in f:
            line=line.rstrip("\r\n")
            if line.startswith(("!platform_table_end","!annotation_table_end")):break
            if not found:
                if line in ("!platform_table_begin","!annotation_table_begin"):continue
                a=next(csv.reader([line],delimiter="\t"))
                if "ID" in a and ("GENE" in a or "Gene ID" in a):
                    header=a;found=True
                continue
            a=next(csv.reader([line],delimiter="\t"))
            require(len(a)==len(header),"Platform width")
            d=dict(zip(header,a));n+=1
            require(n<=100000,"Platform row cap")
            probe=d["ID"];gene=d.get("GENE",d.get("Gene ID","")).strip()
            if not re.fullmatch(r"A_[0-9]+_P[0-9]+",probe):continue
            if "CONTROL_TYPE" in d and d["CONTROL_TYPE"].strip().upper()!="FALSE":continue
            if re.fullmatch(r"[0-9]+",gene):rows[probe].add(gene)
        require(found and n>10000,"Exact platform annotation schema or row count")
    return {k:next(iter(v)) for k,v in rows.items() if len(v)==1},{"table_rows":n,"header":header,"mapped_probes":sum(len(v)==1 for v in rows.values())}
def score(gene_rows,panel):
    g=sorted(gene_rows)
    require(len(g)>=5000,"Mapped expression gene universe")
    mat=np.stack([gene_rows[x] for x in g])
    require(np.isfinite(mat).all(),"Nonfinite expression")
    ranks=rankdata(mat,axis=0,method="average")
    ranks=(ranks-.5)/len(g)
    index={x:i for i,x in enumerate(g)}
    by={s:[index[x] for x,v in panel.items() if v==s and x in index] for s in (-1,1)}
    full={s:sum(v==s for v in panel.values()) for s in (-1,1)}
    require(all(len(by[s])>=max(20,math.ceil(.7*full[s])) for s in by),"Per-direction coverage <70% or <20")
    values=ranks[by[1],:].mean(axis=0)-ranks[by[-1],:].mean(axis=0)
    return values,{"gene_universe":len(g),"mapped_up":len(by[1]),"mapped_down":len(by[-1]),"total_up":full[1],"total_down":full[-1]}
def describe(v):
    return {"n":len(v),"mean":float(np.mean(v)),"median":float(np.median(v)),"sd":float(np.std(v,ddof=1)) if len(v)>1 else None}
def trend(values,labels,blocks):
    # Fixed ordinal slope, within age x sex strata. 20,000 reproducible label permutations.
    groups=[np.array([i for i,b in enumerate(blocks) if b==q]) for q in sorted(set(blocks))]
    x=np.asarray(labels,dtype=float);y=np.asarray(values,dtype=float)
    def center(a):
        b=a.copy()
        for idx in groups:b[idx]-=b[idx].mean()
        return b
    xc=center(x);yc=center(y);ss=float(xc@xc)
    require(ss>0 and len(x)>len(groups)+1,"Insufficient within-stratum response variation")
    beta=float(xc@yc/ss)
    df=len(x)-len(groups)-1
    se=float(np.sqrt(np.sum((yc-beta*xc)**2)/df/ss))
    crit=float(student_t.ppf(.975,df))
    rng=np.random.default_rng(20261125);extreme=0;B=20000
    for _ in range(B):
        a=x.copy()
        for idx in groups:a[idx]=rng.permutation(a[idx])
        # Within-stratum permuting preserves group mean and variance.
        b=float(center(a)@yc/ss)
        extreme+=abs(b)>=abs(beta)-1e-14
    return {"slope_per_ordinal_category":beta,"model_based_ci95":[beta-crit*se,beta+crit*se],
            "permutation_p_two_sided":(extreme+1)/(B+1),"permutations":B,
            "strata":len(groups),"n":len(x),"within_stratum_informative":sum(len(set(x[idx]))>1 for idx in groups)}
def baseline(root,panel_path):
    p=read_panel(panel_path,"entrez")
    m=root/"baseline.gz";src=download(series("GSE42507")+"/matrix/GSE42507_series_matrix.txt.gz",m,20000000,META_42507)
    data,meta,chars=matrix(m,"GSE42507",44,True)
    require(len(data)>=40000,"Unexpected complete baseline features")
    sources=[
        (BASE+"/platforms/GPL6nnn/GPL6480/annot/GPL6480.annot.gz",root/"platform.gz"),
        ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6480&targ=self&form=text&view=full",root/"platform.txt")]
    failures=[];mapping=None
    for url,path in sources:
        try:
            ps=download(url,path,50000000)
            mapping,ann=platform(path);break
        except (OSError,ValueError,EOFError,UnicodeError,csv.Error) as e:failures.append({"url":url,"reason":str(e)[:180]})
    if mapping is None:return {"status":"HOLD_PLATFORM","source":src,"failures":failures}
    by=collections.defaultdict(list)
    for probe,v in data.items():
        if probe in mapping:by[mapping[probe]].append(v)
    g={k:np.median(np.stack(v),axis=0) for k,v in by.items()}
    vals,coverage=score(g,p)
    # Expression score construction above never uses group labels.
    for k in ("cluster","age","sex"):require(len(chars.get(k,{}))==44,"Missing clinical label "+k)
    cmap={"non":0,"non-responder":0,"nonresponder":0,"non-responders":0,"mod":1,"modest":1,"modest-responder":1,"mod-responder":1,"modest responder":1,"modest-responders":1,"xtr":2,"extreme":2,"extreme-responder":2,"xtr-responder":2,"extreme responder":2,"extreme-responders":2}
    lab=[chars["cluster"][i].casefold().strip() for i in range(44)]
    print("Source cluster category counts",dict(collections.Counter(lab)),flush=True)
    require(all(x in cmap for x in lab),"Unrecognized source cluster categories")
    labels=np.array([cmap[x] for x in lab])
    require(list(np.bincount(labels,minlength=3))==[15,14,15],"Original group counts")
    ages=[chars["age"][i].casefold().strip() for i in range(44)]
    sexes=[chars["sex"][i].casefold().strip() for i in range(44)]
    require(len(set(ages))==2 and set(sexes)=={"male","female"},"Age/sex category contract")
    blocks=list(zip(ages,sexes))
    return {"status":"BASELINE_FIXED_SIGNATURE_ORDINAL_ASSOCIATION_NOT_PREDICTION",
      "matrix_source":src,"platform_source":ps,"platform_schema":ann,"previous_platform_failures":failures,
      "coverage":coverage,"source_group_counts":dict(collections.Counter(lab)),
      "age_category_count":len(set(ages)),"sex_counts":dict(collections.Counter(sexes)),
      "groups":{name:describe(vals[labels==i]) for i,name in enumerate(("non","modest","extreme"))},
      "primary_age_sex_stratified_trend":trend(vals,labels,blocks),
      "unadjusted_sensitivity":trend(vals,labels,["all"]*44),
      "limits":["Selected response clusters, not a representative prospective sample or continuous CSA regression.",
        "Score is per-sample signed mean percentile rank (up minus down), not acute response or calibrated risk.",
        "Permutation assumes within-stratum exchangeability, not randomized response groups; no causal inference.",
        "Ordinal coding assumes equally spaced categories for trend, not equally spaced hypertrophy. No model tuning or prediction."]}
def canonical(v):
    m=re.fullmatch(r"([A-Z]+)_?([0-9]+)",v)
    require(m is not None,"Unrecognized sample title")
    return m.groups()
def chronic(root,panel_path):
    p=read_panel(panel_path,"ensembl")
    m=root/"chronic_meta.gz";src=download(series("GSE277819")+"/matrix/GSE277819_series_matrix.txt.gz",m,1000000,META_277819)
    _,meta,chars=matrix(m,"GSE277819",81,False)
    q=root/"chronic_counts.gz";cs=download(series("GSE277819")+"/suppl/GSE277819_raw_counts.txt.gz",q,10000000,COUNTS_277819)
    titles=meta["!Sample_title"][0];title_index={canonical(v):i for i,v in enumerate(titles)}
    require(len(title_index)==81,"Unique metadata titles")
    g={}
    with gzip.open(q,"rt",encoding="utf-8-sig",newline="") as f:
        rd=csv.reader(f,delimiter="\t");h=next(rd)
        require(len(h)==81 and {canonical(v) for v in h}==set(title_index),"Exact frozen source count headers")
        order=[title_index[canonical(v)] for v in h];inv=np.argsort(order)
        for a in rd:
            require(len(a)==82 and re.fullmatch(r"ENSG[0-9]+(?:\.[0-9]+)?",a[0]),"Count identity/width")
            k=a[0].split(".")[0];require(k not in g,"Stable gene duplicate")
            v=np.asarray([float(x) for x in a[1:]],dtype=float)
            require(np.isfinite(v).all() and np.all(v>=0) and np.all(v==np.floor(v)),"Integer nonnegative counts")
            g[k]=v[inv]
    require(len(g)==65217,"Exact count rows")
    vals,coverage=score(g,p)
    groups=collections.defaultdict(dict)
    for i in range(81):
        d=chars["participant id"][i];arm=chars["condition"][i].strip().casefold()
        require(arm in {"baseline","1set","4set"} and arm not in groups[d],"Exact participant arm")
        groups[d][arm]=vals[i]
    require(len(groups)==27 and all(set(a)=={"baseline","1set","4set"} for a in groups.values()),"27 complete participants")
    changes={k:np.array([a[k]-a["baseline"] for a in groups.values()]) for k in ("1set","4set")}
    delta=changes["4set"]-changes["1set"];pos=int((delta>0).sum());neg=int((delta<0).sum())
    return {"status":"CHRONIC_FIXED_SIGNATURE_VOLUME_CONTRAST_NOT_HYPERTROPHY",
        "metadata_source":src,"count_source":cs,"coverage":coverage,"n_participants":27,
        "baseline_to_post_changes":{k:describe(v) for k,v in changes.items()},
        "paired_4set_minus_1set":{**describe(delta),"positive":pos,"negative":neg,"ties":27-pos-neg,
           "exact_two_sided_sign_p":float(binomtest(pos,pos+neg,.5).pvalue) if pos+neg else 1.},
        "limits":["Outcome-selected cohort. Both limbs from same participant; paired molecular contrast only.",
          "Within-sample rank score depends on measured feature universe; cross-platform absolute scores not comparable.",
          "Post-training expression is not a baseline predictor. No continuous qCSA linkage or prediction."]}
def main():
    ap=argparse.ArgumentParser();ap.add_argument("mode",choices=["baseline","chronic","selftest"]);ap.add_argument("--root",type=Path);ap.add_argument("--panel",type=Path);ap.add_argument("--out",type=Path);a=ap.parse_args()
    if a.mode=="selftest":
        assert canonical("A_102")==canonical("A102")
        p=trend(np.array([0.,1.,2.,3.]),np.array([0,1,0,1]),["a","a","b","b"])
        assert abs(p["slope_per_ordinal_category"]-1)<1e-10
        print("Two focused synthetic assertions passed");return
    a.root.mkdir(parents=True,exist_ok=True)
    try:
        r=baseline(a.root,a.panel) if a.mode=="baseline" else chronic(a.root,a.panel)
    except (OSError,ValueError,KeyError,TypeError,EOFError,UnicodeError,csv.Error) as e:
        r={"status":"HOLD_SOURCE_OR_ANALYSIS_GATE","reason":str(e)[:300]}
    r["panel_sha256"]=hashlib.sha256(a.panel.read_bytes()).hexdigest()
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2,allow_nan=False)+"\n")
    print(json.dumps(r,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
