#!/usr/bin/env python3
"""Fixed-panel GSE313736 pilot, source alignment before expression scoring.

Only aggregate pool summaries are exported. No clinical fields, genotype, cell
barcodes, donor IDs, or donor-gene matrices are emitted.
"""
import argparse,collections,csv,gzip,hashlib,json,math,re,statistics
from pathlib import Path
from ocular_frozen_241_external import read_panel

META_SHA="41e4a6bf5792fe58e023b92f0129119dcf2b49877f62dd580b9dab5de48f5192"
# Exact GEO GSM records, vehicle and 500 nM dex pool C.
POOLS={"C":{"no-dex":"GSM9373278","dex":"GSM9373287","official_donors":{24,25,26,28,31,32,34,76}}}
CUTS=(50,100,200)
def require(ok,msg):
    if not ok:raise ValueError(msg)
def digest(p):
    with p.open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()
def norm(s):return re.sub(r"[^a-z0-9]","",s.casefold())
def rows(path):
    with gzip.open(path,"rt",encoding="utf-8-sig",newline="") as f:
        return list(csv.reader(f,delimiter="\t"))
def pool_library(arm,pool,gsm):
    return {norm("HTMC_pool-"+("dex" if arm=="dex" else "")+pool+"_GRCh38")}
def barcode_index(source,selected):
    # Exact barcode, or one explicit reversible Cell Ranger -1 suffix convention.
    require(len(source)==len(set(source)),"Duplicate matrix barcodes")
    require(len(selected)==len(set(selected)),"Duplicate metadata barcodes")
    strategies={"exact":lambda x:x,"remove_terminal_dash1":lambda x:re.sub(r"-1$","",x)}
    matched=[]
    for name,fn in strategies.items():
        ids=[fn(x) for x in source];wanted=[fn(x) for x in selected]
        if len(ids)!=len(set(ids)) or len(wanted)!=len(set(wanted)):continue
        index={v:i for i,v in enumerate(ids)}
        if all(v in index for v in wanted):
            matched.append((name,[index[v] for v in wanted]))
    require(matched,"No complete exact/reversible barcode alignment")
    require(all(ix==matched[0][1] for _,ix in matched),"Ambiguous barcode alignment")
    return matched[0]
def meta(root,pool):
    path=root/"metadata.gz"
    require(path.stat().st_size==6592040 and digest(path)==META_SHA,"Public metadata identity changed")
    spec=POOLS[pool];by={arm:[] for arm in ("no-dex","dex")};seen=set();libnames=collections.defaultdict(set)
    global_pool=collections.defaultdict(lambda:collections.defaultdict(set));n=0
    with gzip.open(path,"rt",encoding="utf-8-sig",newline="") as f:
        rd=csv.DictReader(f,delimiter="\t")
        require(set(("library_id","DEX Treatment","barcode","droplet_type","sampleid"))<=set(rd.fieldnames or []),"Metadata schema")
        for r in rd:
            n+=1
            arm=r["DEX Treatment"];lib=r["library_id"];barcode=r["barcode"];sample=r["sampleid"]
            m=re.fullmatch(r"HTMC[-_]?(\d+)",sample,re.I)
            require(m is not None and r["droplet_type"]=="SNG" and arm in by,"Metadata identity/QC")
            k=(lib,barcode);require(k not in seen,"Duplicate library-barcode");seen.add(k)
            global_pool[sample][arm].add(lib);libnames[arm].add(lib)
            if norm(lib) in pool_library(arm,pool,spec[arm]):
                require(int(m.group(1)) in spec["official_donors"],"Metadata donor outside publisher pool")
                by[arm].append((sample,barcode))
    require(n==169398 and len(global_pool)==59,"Frozen metadata coverage changed")
    require(all(set(x)==set(by) and all(len(v)==1 for v in x.values()) for x in global_pool.values()),"Sample not confined to one library per arm")
    # Library names are public technical labels, not donor identifiers.
    print("Source library labels",json.dumps({k:sorted(v) for k,v in libnames.items()}),flush=True)
    require(all(by.values()),"Expected exact pool C library identifiers not found")
    ids={arm:{d for d,b in v} for arm,v in by.items()}
    require(ids["dex"]==ids["no-dex"],"Pilot donor sets differ")
    return by,sorted(ids["dex"])
def load_arm(root,arm,records,donors,panel):
    import numpy as np
    from scipy.io import mmread
    from scipy.sparse import issparse
    from threadpoolctl import threadpool_limits
    stem=root/arm
    paths={k:Path(str(stem)+"_"+k+ext) for k,ext in (("barcodes",".tsv.gz"),("features",".tsv.gz"),("matrix",".mtx.gz"))}
    b=rows(paths["barcodes"]);features=rows(paths["features"])
    require(all(len(x)==1 and re.fullmatch("[ACGT]+-[0-9]+",x[0]) for x in b),"Source barcode schema")
    require(all(len(x)==3 for x in features),"Source features schema")
    require(all(x[2]=="Gene Expression" for x in features),"Non-GEX feature requires explicit denominator")
    require(all(re.fullmatch(r"ENSG\d+(?:\.\d+)?",x[0]) for x in features),"Unexpected feature ID")
    fids=[x[0].split(".")[0] for x in features]
    require(len(set(fids))==len(fids),"Duplicate stable Ensembl ID")
    panel_ids=sorted(set(fids)&set(panel))
    ecm=[g for g in panel_ids if panel[g]["ecm_member"]]
    require(len(panel_ids)>=math.ceil(.9*241) and len(ecm)==16,"Frozen panel mapping gate failed")
    bc=[r[0] for r in b]
    mode,selected=barcode_index(bc,[barcode for _,barcode in records])
    donor_index={d:i for i,d in enumerate(donors)}
    group=np.full(len(bc),-1,dtype=np.int32)
    for (d,_),i in zip(records,selected):group[i]=donor_index[d]
    ncells=np.bincount(group[group>=0],minlength=len(donors))
    require(int(ncells.sum())==len(records) and np.min(ncells)>0,"Selected donor/cell mismatch")
    gene_index={g:i for i,g in enumerate(panel_ids)}
    rowmap=np.array([gene_index.get(g,-1) for g in fids],dtype=np.int32)
    # Read original sparse matrix; bound memory by processing entries in chunks.
    with gzip.open(paths["matrix"],"rb") as f:
        with threadpool_limits(limits=2):
            m=mmread(f)
        require(f.read(1)==b"","Trailing matrix payload")
    require(issparse(m) and m.shape==(len(features),len(bc)),"Matrix feature/barcode shape mismatch")
    coo=m.tocoo(copy=False)
    require(coo.nnz>0 and np.issubdtype(coo.dtype,np.integer),"Expected integer UMI count matrix")
    require(coo.nnz<=300_000_000,"Matrix nonzero bound")
    totals=np.zeros(len(donors),dtype=np.float64)
    gene_counts=np.zeros(len(panel_ids)*len(donors),dtype=np.float64)
    for start in range(0,coo.nnz,2_000_000):
        sl=slice(start,min(start+2_000_000,coo.nnz))
        val=coo.data[sl];r=rowmap[coo.row[sl]];g=group[coo.col[sl]]
        require(np.min(val)>=0,"Negative UMI")
        keep=g>=0
        totals+=np.bincount(g[keep],weights=val[keep],minlength=len(donors))
        keep&=r>=0
        gene_counts+=np.bincount(r[keep]*len(donors)+g[keep],weights=val[keep],minlength=len(gene_counts))
    counts=gene_counts.reshape(len(panel_ids),len(donors))
    require(np.min(totals)>0 and np.all(counts.sum(axis=0)<=totals) and np.max(totals)<2**53,"Pseudobulk totals invalid")
    expression=np.log2(1+counts*1e6/totals[None,:])
    summary={"matrix_shape":[int(x) for x in coo.shape],"matrix_nonzeros":int(coo.nnz),
        "source_bytes":{k:p.stat().st_size for k,p in paths.items()},
        "source_sha256":{k:digest(p) for k,p in paths.items()},
        "barcode_alignment":mode,"metadata_cells_matched":len(records),"metadata_cells_unmatched":0,
        "other_matrix_cells_excluded":len(bc)-len(records),"sample_ids":len(donors),
        "mapped_panel_genes":len(panel_ids),"mapped_ECM_genes":len(ecm),
        "min_cells_per_sample":int(ncells.min()),"max_cells_per_sample":int(ncells.max()),
        "all_selected_UMI":int(totals.sum())}
    return panel_ids,fids,ncells,expression,summary
def pilot(root,out,pool):
    import numpy as np
    panel=read_panel(root/"panel.zip")
    by,donors=meta(root,pool)
    # Establish all design/identity constraints before interpreting expression.
    a=load_arm(root,"no-dex",by["no-dex"],donors,panel)
    b=load_arm(root,"dex",by["dex"],donors,panel)
    require(a[0]==b[0] and a[1]==b[1],"Feature universe differs between arms")
    ids=a[0];eligible=np.minimum(a[2],b[2]);d=b[3]-a[3]
    signs=np.array([1 if panel[g]["acute_log2fc"]>0 else -1 for g in ids])
    aligned=d*signs[:,None];panels={"all_241":list(range(len(ids))),"ECM_16":[i for i,g in enumerate(ids) if panel[g]["ecm_member"]]}
    result={}
    for cut in CUTS:
        mask=eligible>=cut
        sub={}
        for name,ix in panels.items():
            if not np.any(mask):sub[name]={"status":"HOLD_NO_ELIGIBLE_PAIR"};continue
            v=aligned[ix,:][:,mask].mean(axis=0)
            sub[name]={"pairs":int(mask.sum()),"positive":int((v>0).sum()),"negative":int((v<0).sum()),"ties":int((v==0).sum()),
                       "mean_signed_log2_cpm_difference":float(v.mean()),"median_signed_log2_cpm_difference":float(np.median(v)),
                       "genes_with_positive_mean_signed_change":int((aligned[ix,:][:,mask].mean(axis=1)>0).sum()),
                       "mapped_genes":len(ix)}
        result[str(cut)]=sub
    r={"status":"ONE_POOL_PILOT_DESCRIPTIVE_NOT_REPLICATION","pool":pool,"pair_count":len(donors),
       "arms":{"vehicle":a[4],"dex":b[4]},"scores_by_min_cells":result,
       "limits":["Pilot only, not all nine paired pools. No donor-level confirmatory P or independent biological replication claim.",
                 "Retained all released SNG rows with exact metadata/matrix alignment; no reclassification or response-based exclusion.",
                 "Source pool is exposure unit; donors within pool share library and treatment. Across-pool inference required.",
                 "UMI counts summed per donor and arm, full GEX library denominator, original signs; log2(CPM+1), not clinical IOP or donor genotype.",
                 "No donor identifiers, clinical fields, raw matrix or per-donor expression exported. Local source digests, not publisher checksums."]}
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(r,indent=2,allow_nan=False)+"\n");print(json.dumps(r,indent=2,allow_nan=False))
def selftest():
    assert barcode_index(["AA-1","TT-1"],["TT-1"])[1]==[1]
    assert barcode_index(["AA-1","TT-1"],["TT"])[1]==[1]
    try:barcode_index(["AA-1","AA"],["AA","AA-1","TT"])
    except ValueError:pass
    else:raise AssertionError("Invalid barcode alignment accepted")
    assert norm("pool-dexC")=="pooldexc"
    print("Four focused barcode/label checks passed")
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["selftest","pilot"]);p.add_argument("--pool",choices=sorted(POOLS),default="C");p.add_argument("--root",type=Path,default=Path("."));p.add_argument("--out",type=Path,default=Path("outputs/pilot.json"));a=p.parse_args()
    if a.mode=="selftest":selftest()
    else:pilot(a.root,a.out,a.pool)
