#!/usr/bin/env python3
"""Exact-intersection GSE313736 validation. Aggregate output only; no donor data export."""
import argparse,collections,concurrent.futures,csv,gzip,hashlib,json,math,re,statistics,subprocess
from pathlib import Path
from ocular_gse313736_pool_pilot import load_arm,META_SHA,digest,rows,require,norm
from ocular_gse313736_all_library_preflight import fetch,sample_records,SERIES
from ocular_frozen_241_external import read_panel,PIN

POOLS="ABCDEFGHI"
ARMS=("no-dex","dex")
EXPECTED={"A":8,"B":7,"C":7,"D":7,"E":7,"F":6,"G":7,"H":6,"I":4}
PRIMARY="retention80_min100"
MODES={PRIMARY:(.8,100),"retention80_min50":(.8,50),"retention80_min200":(.8,200),"all_exact_min100":(0.,100)}
def write(out,obj):
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n")
    print(json.dumps(obj,indent=2,allow_nan=False),flush=True)
def metadata(path):
    require(path.stat().st_size==6592040 and digest(path)==META_SHA,"Pinned metadata identity changed")
    labels={norm("HTMC_pool-"+("dex" if a=="dex" else "")+p+"_GRCh38"):(p,a) for p in POOLS for a in ARMS}
    out=collections.defaultdict(list);seen=set();locations=collections.defaultdict(dict)
    with gzip.open(path,"rt",encoding="utf-8-sig",newline="") as f:
        r=csv.DictReader(f,delimiter="\t")
        require(set(("library_id","DEX Treatment","barcode","droplet_type","sampleid"))<=set(r.fieldnames or []),"Metadata schema")
        for x in r:
            key=labels.get(norm(x["library_id"]))
            require(key is not None and x["DEX Treatment"]==key[1] and x["droplet_type"]=="SNG","Metadata arm/library/QC")
            sample=x["sampleid"];bc=x["barcode"]
            require(re.fullmatch(r"HTMC[-_]?(\d+)",sample,re.I) and re.fullmatch(r"[ACGT]{16}-1",bc),"Metadata exact identity")
            k=(key,bc);require(k not in seen,"Duplicate metadata library-barcode");seen.add(k)
            require(key[1] not in locations[sample] or locations[sample][key[1]]==key[0],"Repeated donor in multiple pools")
            locations[sample][key[1]]=key[0]
            out[key].append((sample,bc))
    require(len(seen)==169398 and len(locations)==59 and len(out)==18,"Fixed cohort identity")
    require(all(set(x)==set(ARMS) and x["dex"]==x["no-dex"] for x in locations.values()),"Donor pairing")
    for p in POOLS:
        a,b=[{x[0] for x in out[p,arm]} for arm in ARMS]
        require(a==b and len(a)==EXPECTED[p],"Frozen paired pool coverage")
    return out
def exact_intersection(bc,records):
    require(len(bc)==len(set(bc)) and all(re.fullmatch(r"[ACGT]{16}-1",x) for x in bc),"Unexpected or duplicate matrix barcode")
    source=set(bc)
    require(len({b for d,b in records})==len(records),"Duplicate metadata barcode")
    keep=[(d,b) for d,b in records if b in source]
    initial=collections.Counter(d for d,b in records);matched=collections.Counter(d for d,b in keep)
    require(len(keep)>0 and len(keep)/len(records)>=.90,"Library exact intersection <90 percent")
    return keep,initial,matched
def feature_ids(path):
    f=rows(path)
    require(all(len(x)==3 and x[2]=="Gene Expression" and re.fullmatch(r"ENSG\d+(?:\.\d+)?",x[0]) for x in f),"Unexpected GEX feature schema")
    ids=[x[0].split(".")[0] for x in f]
    require(len(ids)==len(set(ids)),"Duplicate stable feature ID")
    return ids
def library(root,arm,files,records):
    b=root/(arm+"_barcodes.tsv.gz");f=root/(arm+"_features.tsv.gz")
    bs=fetch(files["barcodes"],b,3000000);fs=fetch(files["features"],f,3000000)
    bcs=rows(b)
    require(all(len(x)==1 for x in bcs),"Matrix barcode table schema")
    ids=feature_ids(f)
    keep,initial,matched=exact_intersection([x[0] for x in bcs],records)
    feature_hash=hashlib.sha256(("\n".join(ids)+"\n").encode()).hexdigest()
    summary={"source_metadata_cells":len(records),"exact_matched_cells":len(keep),
         "missing_metadata_cells":len(records)-len(keep),"retention":len(keep)/len(records),
         "matrix_barcodes":len(bcs),"features":len(ids),"feature_order_sha256":feature_hash,
         "barcodes_sha256":bs,"features_sha256":fs,
         "donors_with_match":len(matched),"minimum_donor_retention":min(matched[d]/n for d,n in initial.items()),
         "donor_arms_below_80pct":sum(matched[d]<.8*n for d,n in initial.items()),
         "minimum_matched_cells":min(matched[d] for d in initial)}
    return keep,initial,matched,summary,ids
def eligibility(a,b):
    require(set(a[1])==set(b[1]),"Original sample sets differ")
    result={}
    for mode,(ret,cut) in MODES.items():
        result[mode]={d for d in a[1] if all(x[2][d]>=cut and x[2][d]>=ret*x[1][d] for x in (a,b))}
    return result
def preflight(root,out):
    meta=root/"metadata.gz";soft=root/"series.soft.gz"
    fetch(SERIES+"/suppl/GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",meta,10000000)
    soft_sha=fetch(SERIES+"/soft/GSE313736_family.soft.gz",soft,25000000)
    specs=sample_records(soft);md=metadata(meta)
    def check(p):
        sub=root/p;sub.mkdir(exist_ok=True)
        arms={}
        for arm in ARMS:
            spec=specs[p,arm]
            require(all(int(re.fullmatch(r"HTMC[-_]?(\d+)",d,re.I).group(1)) in spec["donors"] for d,b in md[p,arm]),"Sample not in GEO pool")
            arms[arm]=library(sub,arm,spec["files"],md[p,arm])
        a,b=[arms[arm] for arm in ARMS]
        require(a[4]==b[4],"Pool feature universe differs")
        eligible=eligibility(a,b)
        return {"pool":p,"source_metadata_pairs":EXPECTED[p],
            "eligible_pairs":{k:len(v) for k,v in eligible.items()},
            "arms":{arm:{**arms[arm][3],"gsm":specs[p,arm]["gsm"],"files":specs[p,arm]["files"],
                    "original_source_donor_count":len(specs[p,arm]["donors"])} for arm in ARMS},
            "original_GEO_donor_sets_equal":specs[p,"dex"]["donors"]==specs[p,"no-dex"]["donors"]}
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:pool_results=list(ex.map(check,POOLS))
    hashes={a["feature_order_sha256"] for p in pool_results for a in p["arms"].values()}
    totals={arm:sum(p["arms"][arm]["exact_matched_cells"] for p in pool_results) for arm in ARMS}
    passed=len(hashes)==1 and all(p["eligible_pairs"][PRIMARY]>=3 for p in pool_results)
    r={"status":"EXACT_INTERSECTION_PREFLIGHT_PASS" if passed else "HOLD_PREFLIGHT",
         "metadata_sha256":META_SHA,"source_soft_sha256":soft_sha,"source_metadata_cells":169398,
         "matched_cells":totals,"source_metadata_pairs":59,"pool_count":9,
         "eligible_pairs":{k:sum(p["eligible_pairs"][k] for p in pool_results) for k in MODES},
         "feature_order_sha256":next(iter(hashes)) if len(hashes)==1 else None,"pools":pool_results,
         "limits":["Intersection of exact source library-barcode keys, not fuzzy matching or missing-value imputation.",
                   "Thresholds set after metadata-only diagnostics but before expression scores. Missingness may be informative.",
                   "Each of nine pools must have at least three paired samples meeting primary thresholds. No outcome-based exclusions.",
                   "No clinical or donor IDs or cell barcodes emitted. No expression inference in preflight."]}
    write(out,r);require(passed,"Preflight gate failed")
def matrix_download(url,path,expected):
    require(url.startswith("https://ftp.ncbi.nlm.nih.gov/geo/samples/"),"Matrix source host")
    subprocess.run(["curl","-fLsS","--connect-timeout","20","--max-time","480","--max-filesize","1400000000",url,"-o",str(path)],check=True,timeout=500)
    with gzip.open(path,"rt",encoding="ascii") as f:
        h=f.readline().strip().lower()
        require(h=="%%matrixmarket matrix coordinate integer general","Unexpected matrix type")
        for i in range(100):
            line=f.readline()
            if not line.startswith("%"):break
        else:raise ValueError("Matrix comments exceed bound")
        shape=list(map(int,line.split()))
        require(len(shape)==3 and shape[:2]==[expected["features"],expected["matrix_barcodes"]] and 0<shape[2]<=230_000_000,"Unexpected matrix dimensions or allocation bound")
    print("Verified matrix dimensions",path.name,shape,flush=True)
def score(root,gate_path,pool,out):
    import numpy as np
    gate=json.loads(gate_path.read_text())
    require(gate["status"]=="EXACT_INTERSECTION_PREFLIGHT_PASS" and gate["metadata_sha256"]==META_SHA and gate["pool_count"]==9,"Preflight not passed")
    q=[p for p in gate["pools"] if p["pool"]==pool];require(len(q)==1,"Missing pool preflight");q=q[0]
    mdpath=root/"metadata.gz";fetch(SERIES+"/suppl/GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",mdpath,10000000)
    md=metadata(mdpath);panel=read_panel(root/"panel.zip")
    libs={}
    for arm in ARMS:
        g=q["arms"][arm];x=library(root,arm,g["files"],md[pool,arm]);libs[arm]=x
        require(all(x[3][k]==g[k] for k in x[3]),"Exact join or small-file hashes changed")
    a,b=[libs[arm] for arm in ARMS];eligible=eligibility(a,b)
    require({k:len(v) for k,v in eligible.items()}==q["eligible_pairs"],"Eligibility drift")
    require(len(eligible[PRIMARY])>=3,"Primary pool pair gate")
    # All samples with nonzero exact cells in both arms; masks are frozen pre-expression.
    donors=sorted(set(a[2])&set(b[2]))
    ids=sorted(set(a[4])&set(panel))
    require(len(ids)>=math.ceil(.9*241) and sum(panel[g]["ecm_member"] for g in ids)==16,"Frozen panel mapping")
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        tasks=[ex.submit(matrix_download,q["arms"][arm]["files"]["matrix"],root/(arm+"_matrix.mtx.gz"),q["arms"][arm]) for arm in ARMS]
        for t in tasks:t.result()
    vals={}
    for arm in ARMS:
        r=[(d,bc) for d,bc in libs[arm][0] if d in donors]
        vals[arm]=load_arm(root,arm,r,donors,panel)
        require(vals[arm][0]==ids and vals[arm][1]==a[4],"Matrix feature mapping changed")
        require(all(int(n)==libs[arm][2][d] for d,n in zip(donors,vals[arm][2])),"Matrix retained counts changed")
    d=vals["dex"][3]-vals["no-dex"][3]
    signed=d*np.array([1 if panel[g]["acute_log2fc"]>0 else -1 for g in ids])[:,None]
    ix={"all_241":list(range(len(ids))),"ECM_16":[i for i,g in enumerate(ids) if panel[g]["ecm_member"]]}
    results={}
    for mode,keep in eligible.items():
        mask=[i for i,d in enumerate(donors) if d in keep];results[mode]={}
        for name,indices in ix.items():
            if not mask:results[mode][name]={"status":"HOLD_EMPTY"};continue
            m=signed[np.ix_(indices,mask)]
            v=m.mean(axis=0)
            results[mode][name]={"status":"SCORED","pairs":len(mask),"positive":int((v>0).sum()),"negative":int((v<0).sum()),"ties":int((v==0).sum()),
              "pool_equal_donor_mean":float(v.mean()),"donor_median":float(np.median(v)),
              "mapped_genes":len(indices),"positive_mean_genes":int((m.mean(axis=1)>0).sum())}
    write(out,{"status":"EXACT_INTERSECTION_POOL_SCORED","pool":pool,"panel_sha256":PIN,"metadata_sha256":META_SHA,
       "mapped_panel_id_sha256":hashlib.sha256(("\n".join(ids)+"\n").encode()).hexdigest(),
       "feature_order_sha256":gate["feature_order_sha256"],"source_metadata_pairs":EXPECTED[pool],
       "eligible_pairs":q["eligible_pairs"],"preflight_arms":q["arms"],
       "matrix_qc":{arm:vals[arm][4] for arm in ARMS},"results":results,
       "limits":["Exact source intersection, full GEX denominator, fixed panel signs, all SNG cell states.",
                 "No per-donor expression/identifiers or clinical data exported. No individual pool significance test.",
                 "Counts are conditional on public metadata and count-matrix overlap, not clinical IOP response."]})
def signp(values):
    p=sum(x>0 for x in values);n=sum(x<0 for x in values);k=p+n
    return min(1.,2*sum(math.comb(k,i) for i in range(min(p,n)+1))/2**k) if k else 1.
def holm(p):
    order=sorted(p,key=lambda k:(p[k],k));adj={};prev=0.
    for i,k in enumerate(order):
        prev=max(prev,min(1.,(len(order)-i)*p[k]));adj[k]=prev
    return adj
def combine(root,out):
    paths=list(root.rglob("pool-*.json"));data=[json.loads(p.read_text()) for p in paths]
    require(len(data)==9 and {x["pool"] for x in data}==set(POOLS),"Require all nine unique scored pools")
    require(all(x["status"]=="EXACT_INTERSECTION_POOL_SCORED" and x["panel_sha256"]==PIN and x["metadata_sha256"]==META_SHA for x in data),"Pool provenance")
    require(len({x["mapped_panel_id_sha256"] for x in data})==len({x["feature_order_sha256"] for x in data})==1,"Common panel/feature universe")
    data.sort(key=lambda x:x["pool"]);result={}
    for mode in MODES:
        r={}
        for name in ("all_241","ECM_16"):
            rows_=[x["results"][mode][name] for x in data]
            require(all(x["status"]=="SCORED" and x["pairs"]>=3 for x in rows_),"Missing eligible pool for "+mode)
            v=[x["pool_equal_donor_mean"] for x in rows_];loo=[statistics.mean(v[:i]+v[i+1:]) for i in range(9)]
            r[name]={"pools":9,"donor_pairs":sum(x["pairs"] for x in rows_),"pool_positive":sum(x>0 for x in v),
                "pool_negative":sum(x<0 for x in v),"pool_ties":sum(x==0 for x in v),
                "donor_positive":sum(x["positive"] for x in rows_),"donor_negative":sum(x["negative"] for x in rows_),"donor_ties":sum(x["ties"] for x in rows_),
                "equal_pool_mean":statistics.mean(v),"pool_median":statistics.median(v),
                "minimum_leave_one_pool_out_mean":min(loo),"maximum_leave_one_pool_out_mean":max(loo),
                "all_leave_one_pool_out_positive":min(loo)>0,"exact_two_sided_pool_sign_p":signp(v),
                "pool_scores":[{"pool":x["pool"],"pairs":y["pairs"],"mean":z,"positive_donors":y["positive"],"mapped_genes":y["mapped_genes"]} for x,y,z in zip(data,rows_,v)]}
        p=holm({k:v["exact_two_sided_pool_sign_p"] for k,v in r.items()})
        for k in r:r[k]["holm_two_panels"]=p[k]
        result[mode]=r
    write(out,{"status":"NINE_POOL_FIXED_PANEL_EXTERNAL_VALIDATION","primary_mode":PRIMARY,
        "source_metadata_pairs":59,"source_metadata_cells":169398,"panel_sha256":PIN,
        "matched_cells":{arm:sum(x["preflight_arms"][arm]["exact_matched_cells"] for x in data) for arm in ARMS},
        "results":result,"limits":["Primary exact two-sided sign test across nine pool means, conditional on independent pool pairs and sign-test null; not 59 independent treatment replicates.",
          "Positive direction defined by frozen prior panel. Holm across two panels; sensitivity modes exploratory, not additional primary hypotheses.",
          "Exact matched subset; missing metadata cells may be informative. Pool E original donor set differed; analysis retains only released matched sample IDs.",
          "All released SNG cell states together; relative expression normalization and cell-composition effects remain possible.",
          "Cohort provenance/overlap and source treatment duration still need final manuscript audit. No clinical IOP outcome, no claim of new causal mechanism."]})
def selftest():
    A="A"*16+"-1";C="C"*16+"-1";G="G"*16+"-1"
    a,i,m=exact_intersection([A,C],[("x",A),("x",C)])
    assert len(a)==2 and i["x"]==m["x"]==2
    try:exact_intersection([A,C],[("x",A),("x",G)])
    except ValueError:pass
    else:raise AssertionError("Missing barcode improperly accepted")
    assert signp([1]*9)==.00390625
    assert signp([1]*8+[-1])==.0390625
    assert holm({"a":.00390625,"b":.00390625})=={"a":.0078125,"b":.0078125}
    print("Five focused exact-join, sign-test and Holm assertions passed")
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["selftest","preflight","score","combine"])
    p.add_argument("--root",type=Path,default=Path("."));p.add_argument("--gate",type=Path)
    p.add_argument("--pool",choices=list(POOLS));p.add_argument("--out",type=Path,default=Path("outputs/result.json"));a=p.parse_args()
    if a.mode=="selftest":selftest()
    elif a.mode=="preflight":preflight(a.root,a.out)
    elif a.mode=="score":score(a.root,a.gate,a.pool,a.out)
    else:combine(a.root,a.out)
