#!/usr/bin/env python3
"""Source-bound quantitative readout. Fixed panels; aggregate public outputs only."""
import argparse,collections,concurrent.futures,csv,gzip,hashlib,io,json,math,re,statistics,sys,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
for d in ("schema_helpers","paired_helpers","prior_helpers"):
    sys.path.insert(0,str(ROOT/d/"scripts/publication"))
import early_track_schema_repair as repair
old=repair.old;base=old.base;require=base.require
ARCH=ROOT/"docs/publication/evidence/early_track_schema_readback_20261125.json"
HCC_GENES=("SOAT1","HMGCR","SQLE","LDLR","FASN","ACACA","ACLY","SREBF2")
SC_ENDPOINTS=(("epithelial","pEMT"),("myeloid","SPP1"),("myeloid","CXCR4"),("T_cell","CXCR4"),("B_cell","CXCR4"))
def prior(mode):
    a=json.loads(ARCH.read_text());require(a["run_id"]==45606353863,"source run drift")
    hits=[x["aggregate"] for x in a["records"] if x["aggregate"]["mode"]==mode];require(len(hits)==1,"prior source archive")
    return hits[0]
def norm(v):return re.sub(r"[^a-z0-9]","",v.lower())
def header_profile(h,keys):
    # Report only counts and recognized non-sample annotation words.
    known={"","gene","geneid","genename","genesymbol","symbol","ensembl","ensemblgeneid","chr","chromosome","start","end","strand","length","unnamed0","x","biotype","description","genetype"}
    return {"width":len(h),"exact_matches":{k:sum(v in set(q) for v in h) for k,q in keys.items()},"recognized_annotation_fields":[norm(v) for v in h if norm(v) in known]}
def exact_columns(h,keys):
    # Exact publisher identifiers only; non-sample columns may be explicit annotation.
    opts=[]
    for k,q in keys.items():
        if len(set(q))!=len(q):continue
        ix=[i for i,v in enumerate(h) if v in set(q)]
        if len(ix)==len(q) and len({h[i] for i in ix})==len(q):
            opts.append((k,ix,[q.index(h[i]) for i in ix]))
    require(opts,"no complete exact GEO sample column set")
    require(all(v[1:]==opts[0][1:] for v in opts),"multiple incompatible sample column mappings")
    return opts[0]
def bulk():
    import numpy as np
    rows,ids,titles,meta=base.geo_metadata("GSE275870");desc,keys=repair.bulk_sample_keys(rows,ids,titles)
    raw,src=base.get(base.geo_folder("GSE275870")+"suppl/GSE275870_counts.csv.gz",15_000_000)
    # The publisher calls the file CSV; detect actual tab/comma/semicolon only by an exact full join.
    opts=[];profiles={}
    fields={"GSM":ids,"description":desc,"title":titles}
    for delim in (",","\t",";"):
        rr=csv.reader(io.TextIOWrapper(gzip.GzipFile(fileobj=io.BytesIO(raw)),encoding="utf-8-sig"),delimiter=delim)
        h=next(rr)
        profiles[repr(delim)]=header_profile(h,fields)
        try:ident,ix,order=exact_columns(h,fields);opts.append((delim,h,ident,ix,order))
        except ValueError:pass
    out={"source":src,"metadata":meta,"header_profiles":profiles}
    if len(opts)!=1:return {"status":"HOLD_EXACT_SAMPLE_COLUMNS",**out}
    delim,h,ident,ix,order=opts[0]
    other=[i for i in range(len(h)) if i not in ix]
    gene_headers={"gene","geneid","ensembl","ensemblgeneid","genesymbol","symbol","genename"}
    genes=[i for i in other if norm(h[i]) in gene_headers]
    # An unnamed first identifier column is valid only if all rows establish unique gene keys.
    if not genes and 0 in other and norm(h[0]) in {"","x","unnamed0"}:genes=[0]
    # Prefer explicitly named Ensembl/gene-ID if both symbol and ID occur.
    stable=[i for i in genes if norm(h[i]) in {"geneid","ensembl","ensemblgeneid"}]
    if len(stable)==1:genes=stable
    require(len(genes)==1,"ambiguous count feature identifier")
    require(len(set(h))==len(h),"duplicate source columns")
    gi=genes[0];out.update(column_mapping=ident,annotation_columns=len(other),feature_header=norm(h[gi]),delimiter=repr(delim))
    rr=csv.reader(io.TextIOWrapper(gzip.GzipFile(fileobj=io.BytesIO(raw)),encoding="utf-8-sig"),delimiter=delim);next(rr)
    totals=np.zeros(108);seen=set();vectors={};ensg=0
    for row in rr:
        require(len(row)==len(h),"count row width drift")
        key=old.strip_version(row[gi].strip());require(key and key not in seen,"duplicate/empty normalized feature")
        seen.add(key);ensg+=bool(re.fullmatch(r"ENSG\d+",key))
        v=np.array([row[i] for i in ix],dtype=float)
        require(np.isfinite(v).all() and (v>=0).all() and (v==np.floor(v)).all(),"not nonnegative integer counts")
        totals+=v;vectors[key]=v
    require(len(seen)>=10000 and (totals>0).all(),"incomplete counts")
    if ensg>=.9*len(seen):mapping,ann=old.gtf_mapping();mode="Ensembl_GRCh38_release_100"
    else:mapping={g:g for g in old.GENES if g in seen};ann=None;mode="exact_unique_symbol"
    buckets=collections.defaultdict(list)
    for j,k in enumerate(order):buckets[keys[k]].append(j)
    pats=sorted({p for p,g in buckets if (p,"PT") in buckets and (p,"MET") in buckets})
    result={}
    for name,gs in old.PANELS.items():
        usable=[g for g in gs if g in mapping and mapping[g] in vectors]
        if len(usable)!=len(gs):
            result[name]={"status":"HOLD_INCOMPLETE_FIXED_PANEL","mapped":len(usable),"required":len(gs),"missing":sorted(set(gs)-set(usable))};continue
        score=np.mean([np.log2(vectors[mapping[g]]/totals*1e6+1) for g in gs],axis=0)
        delta=[float(np.mean(score[buckets[p,"MET"]])-np.mean(score[buckets[p,"PT"]])) for p in pats]
        z=old.sign_summary(delta)
        if "p_exact_sign" in z:z["bonferroni_6"]=min(1,6*z["p_exact_sign"])
        result[name]={**z,"mapped":len(usable),"required":len(gs)}
    return {"status":"ORAL_PAIRED_BULK_FIXED_PANEL",**out,"gene_rows":len(seen),"samples":108,"patients":72,"paired_patients":len(pats),"mapping":mode,"annotation":ann,"results":result,
            "limits":["Same six frozen endpoints as #353; whole tissue and already published paired cohort; sample-average log2(CPM+1), equal patient weight, exact sign test, fixed family 6. No individual expression values or novel causal claim."]}
def sc_feature_universe():
    url=base.geo_folder("GSE198315")+"suppl/GSE198315_features.tsv.gz"
    raw,src=base.get(url,1_000_000)
    rows=[l.split("\t") for l in gzip.decompress(raw).decode("utf8").splitlines() if l]
    require(len(rows)>=10000 and all(1<=len(v)<=3 for v in rows),"feature source invalid")
    ids=set(v[0] for v in rows);symbols={v[1] for v in rows if len(v)>=2}
    require(len(ids)==len(rows),"duplicate combined feature IDs")
    return ids|symbols,src
def strict_orientation(h,a,universe):
    # Prior readback proved equal header/row width for all 20 specimens.
    h=[repair.strip_field(v) for v in h];a=[repair.strip_field(v) for v in a]
    require(len(h)==len(a) and len(a)>=51,"equal complete header/row width required")
    n=len(a)-1;cols=h[1:]
    require(len(set(cols))==n and all(cols),"column keys duplicate/empty")
    # Leading header label may vary, but must not be a barcode. Exact first data feature must
    # occur in the separately published feature inventory; all other columns strict barcodes.
    require(not repair.barcode(h[0]),"unexpected first header barcode")
    if a[0] in universe and all(repair.barcode(v) for v in cols):
        orientation="genes_by_cells"
    elif repair.barcode(a[0]) and old.SC_GENES<=set(cols) and set(cols)<=universe:
        orientation="cells_by_genes"
    else:raise ValueError("unresolved independent feature/barcode orientation")
    repair.numeric(a[1:],n)
    return orientation,cols,a
def sc_headers():
    keys,files,meta=old.sc_mapping();univ,fs=sc_feature_universe()
    obs={}
    for site in ("TC","mLN"):
        p=prior("sc_"+site);obs.update({r["source"]["url"]:r for r in p["records"] if "source" in r})
    records=[]
    for p,g in sorted(keys):
        if g not in {"TC","mLN"}:continue
        _,sid,_=keys[p,g];z={"site":g,"status":"HOLD_HEADER"}
        try:
            url=files[sid];prev=obs[url];raw,src=base.get(url,2_097_152,(0,2_097_151))
            require(src["sha256"]==prev["source"]["sha256"],"prior prefix SHA drift")
            h,a=repair.head_rows(raw);z.update(header_width=len(h),first_row_width=len(a),source=src)
            z["first_header_shape"]={"empty":not bool(h[0]),"letters_only":bool(re.fullmatch(r"[A-Za-z_ .-]+",h[0])),"feature_inventory":h[0] in univ}
            z["barcode_columns"]=sum(repair.barcode(v) for v in h[1:])
            orient,cols,_=strict_orientation(h,a,univ)
            z.update(status="DENSE_HEADER_VERIFIED",orientation=orient,columns=len(cols))
        except (OSError,EOFError,ValueError,TypeError,KeyError,UnicodeError,IndexError) as e:z["reason"]=str(e)[:180]
        records.append(z);time.sleep(.5)
    return {"status":"ORAL_HEADER_FULL_CONTRACT","feature_source":fs,"metadata":meta,"records":records,"verified":sum(x["status"]=="DENSE_HEADER_VERIFIED" for x in records),"target":20}
def sc_read_pair(p,keys,files,univ,prefixes,full_prior):
    result={}
    for site in ("TC","mLN"):
        _,sid,_=keys[p,site];url=files[sid]
        z={"site":site,"status":"HOLD_PILOT"}
        try:
            raw,src=base.get(url,300_000_000);z["source"]=src
            require(hashlib.sha256(raw[:2097152]).hexdigest()==prefixes[url],"prior prefix SHA drift")
            if url in full_prior:require(src["sha256"]==full_prior[url],"previous full matrix SHA drift")
            # Independent source feature identity also for full-body replay, not just prefix.
            out=repair.sc_counts(raw)
            z.update(status="FULL_SOURCE_AND_MARKER_GATE",**out)
        except (OSError,EOFError,ValueError,TypeError,KeyError,UnicodeError,StopIteration,IndexError) as e:z["reason"]=str(e)[:160]
        result[site]=z
    return result
def sc_inputs():
    keys,files,meta=old.sc_mapping();univ,fs=sc_feature_universe()
    repair.identify=lambda h,a:strict_orientation(h,a,univ)
    prefixes={r["source"]["url"]:r["source"]["sha256"] for site in ("TC","mLN") for r in prior("sc_"+site)["records"] if "source" in r}
    full={r["source"]["url"]:r["source"]["sha256"] for r in prior("sc_pilot")["records"] if "source" in r}
    require(len(prefixes)==20 and len(full)==2,"prior source hashes missing")
    return keys,files,meta,univ,fs,prefixes,full
def sc_pilot():
    keys,files,meta,univ,fs,prefixes,full=sc_inputs()
    r=sc_read_pair("P01",keys,files,univ,prefixes,full)
    success=all(v["status"]=="FULL_SOURCE_AND_MARKER_GATE" for v in r.values())
    # Report counts and eligibility, never cell-level/individual expression.
    for z in r.values():
        if "thresholds" in z:z["thresholds"]={th:{"assigned":t["assigned"],"compartments":{k:{"cells":v["cells"],"status":v["status"]} for k,v in t["compartments"].items()}} for th,t in z["thresholds"].items()}
    return {"status":"ORAL_FIRST_PAIR_FULL_REPLAY","full_pair_pass":success,"records":list(r.values()),"metadata":meta,"feature_source":fs}
def sc_full():
    keys,files,meta,univ,fs,prefixes,full=sc_inputs()
    pats=sorted({p for p,g in keys if g=="TC"});require(len(pats)==10,"ten pairs required")
    # Two threads bound CPU/memory; only aggregate results emitted. No raw files in artifacts.
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        rr=list(ex.map(lambda p:sc_read_pair(p,keys,files,univ,prefixes,full),pats))
    complete=[v for v in rr if all(z["status"]=="FULL_SOURCE_AND_MARKER_GATE" for z in v.values())]
    ret={}
    for th in ("2","3"):
        ret[th]={}
        for compartment,gene in SC_ENDPOINTS:
            delta=[]
            for pair in complete:
                a,b=(pair[s]["thresholds"][th]["compartments"][compartment] for s in ("TC","mLN"))
                if any(z["status"]!="DESCRIPTIVE_MARKER_GATE" for z in (a,b)):continue
                if gene=="pEMT":u,v=(z["pemt_mean_log2_cpm_plus1"] for z in (a,b))
                else:u,v=(math.log2(z["cpm"][gene]+1) for z in (a,b))
                delta.append(v-u)
            z=old.sign_summary(delta)
            if th=="2" and "p_exact_sign" in z:z["bonferroni_5"]=min(1,5*z["p_exact_sign"])
            if th=="3":z.pop("p_exact_sign",None)
            ret[th][compartment+"_"+gene]=z
    failed=collections.Counter(z.get("reason","HOLD") for p in rr for z in p.values() if z["status"]!="FULL_SOURCE_AND_MARKER_GATE")
    return {"status":"ORAL_EXPLORATORY_PAIRED_COMPARTMENT_CONTEXT","pairs_target":10,"complete_pairs":len(complete),"failed_source_reasons":dict(failed),"results":ret,
            "qc_cells_by_site":{s:sum(p[s]["qc_cells"] for p in complete) for s in ("TC","mLN")},"metadata":meta,"feature_source":fs,
            "limits":["Primary two-marker gates; three-marker descriptive sensitivity, fixed five primary endpoints, >=50 cells per compartment per site and >=5 complete pairs. Author-coded pairs, not CNV malignant annotation; known context programs, original study already addressed immune/pEMT. Selection and composition confounding remain; no donor-level exports, no novel causal or independently validated association."]}
def hcc_rna():
    rows,ids,titles,aliases,meta=repair.hcc_keys()
    p=prior("hcc_rna")["source"];raw,src=base.get(p["url"],p["bytes"])
    require(src["sha256"]==p["sha256"],"prior HCC source drift")
    rr=csv.reader(io.TextIOWrapper(gzip.GzipFile(fileobj=io.BytesIO(raw)),encoding="utf-8-sig"),delimiter="\t")
    h=next(rr);label,ix,order=exact_columns(h,{"title":titles,"alias":aliases,"GSM":ids})
    require(len(ix)==70 and label=="title","HCC sample mapping drift")
    other=[i for i in range(len(h)) if i not in ix]
    anno=[i for i in other if norm(h[i]) in {"genesymbol","symbol","genename","gene","hgncsymbol","externalgenename"}]
    info={"source":src,"metadata":meta,"non_sample_header":[norm(h[i]) for i in other],"matched_samples":70,"pairs":35}
    if len(anno)!=1:return {"status":"HOLD_GENE_SYMBOL_SCHEMA",**info}
    pool=collections.defaultdict(list);seen=0
    idx={titles[order[j]]:i for j,i in enumerate(ix)};patients=sorted({t[:-1] for t in titles})
    for row in rr:
        seen+=1;require(len(row)==len(h),"FPKM row width drift")
        symbol=row[anno[0]].strip()
        if symbol not in HCC_GENES:continue
        a=[float(row[idx[k+"P"]]) for k in patients];b=[float(row[idx[k+"T"]]) for k in patients]
        require(all(math.isfinite(v) and v>=0 for v in a+b),"invalid FPKM")
        pool[symbol].append([math.log2(y+1)-math.log2(x+1) for x,y in zip(a,b)])
    require(seen==63652,"previous HCC row count drift")
    results={}
    for gene in HCC_GENES:
        v=pool[gene]
        if len(v)!=1:results[gene]={"status":"HOLD_MISSING_OR_NONUNIQUE_GENE","rows":len(v)};continue
        z=old.sign_summary(v[0]);z["bonferroni_8"]=min(1,8*z["p_exact_sign"])
        for k in ("median_delta","mean_delta"):z[k+"_log2_fpkm_plus1"]=z.pop(k+"_log2_cpm_plus1")
        results[gene]=z
    return {"status":"HCC_FIXED_LIPID_CONTEXT_RNA","feature_rows":seen,"results":results,**info,
            "limits":["Fixed eight cholesterol/lipid context genes including original paper SOAT1, not selected by this readout. Paired log2(FPKM+1) tumor-minus-nontumor, exact sign tests, family 8. FPKM not count-based DE; no protein linkage, no clinical outcome or novelty claim."]}
def hcc_protein():
    import openpyxl
    rows,ids,titles,aliases,meta=repair.hcc_keys()
    p=prior("hcc_protein")["publisher_supplement"]["source"];raw,src=base.get(p["url"],p["bytes"])
    require(src["sha256"]==p["sha256"],"supplement source SHA drift")
    # Exact latest official ProteomeXchange revision 4 path.
    pride={}
    try:
        b,pr=base.get("https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/05/PXD006512/",4_000_000)
        names=base.links(b);pride={"source":pr,"raw_links":sum(n.lower().endswith(".raw") for n in names),
            "processed_candidates":[n for n in names if re.fullmatch(r"[A-Za-z0-9_.-]+\.(?:zip|txt|tsv|xlsx)",n,re.I) and re.search(r"maxquant|proteinGroups|quant|result|sdrf",n,re.I)]}
    except (OSError,ValueError,KeyError) as e:pride={"status":"HOLD_DIRECTORY","reason":str(e)[:180]}
    vocab={"gene","geneid","genename","genesymbol","symbol","ensembl","ensemblgeneid","protein","proteinid","proteinids","uniprot","uniprotid","uniprotaccession","accession","description","fpkm","tpm","lfq","intensity","ratio","proteinabundance","ibaq","geneexpression","proteinexpression"}
    inventory=[]
    with zipfile.ZipFile(io.BytesIO(raw)) as zz:
        members=zz.infolist();require(len(members)==18 and sum(v.file_size for v in members)<=600_000_000,"archive inventory drift")
        for f in members:
            name=Path(f.filename).name
            z={"file":name,"bytes":f.file_size}
            if f.filename.lower().endswith(".xlsx") and f.file_size<=80_000_000:
                wb=openpyxl.load_workbook(io.BytesIO(zz.read(f)),read_only=True,data_only=True)
                sheets=[]
                for si,ws in enumerate(wb,1):
                    heads=[]
                    for row in ws.iter_rows(min_row=1,max_row=min(8,ws.max_row or 1),values_only=True):
                        h=[str(v) if v is not None else "" for v in row]
                        t=header_profile(h,{"title":titles,"alias":aliases})
                        anno=[norm(v) for v in h if norm(v) in vocab]
                        if anno or max(t["exact_matches"].values())>0:
                            heads.append({"annotation_terms":anno,"matches":t["exact_matches"]})
                    sheets.append({"index":si,"dimensions":[ws.max_row,ws.max_column],"header_profiles":heads[:3]})
                wb.close();z["sheets"]=sheets
            inventory.append(z)
    return {"status":"HCC_PROTEOME_TABLE_AND_LATEST_PRIDE_GATE","source":src,"pride":pride,"archive_inventory":inventory,
            "limits":["Publisher table labels and sample key counts only, not confirmed RNA/protein analyte types or biological aliquot mapping. PRIDE project revision 4: 2019/05; no raw-file reprocessing or restricted genetic data."]}
def selftest():
    assert exact_columns(["gene","length","b","a"],{"GSM":["a","b"]})==("GSM",[2,3],[1,0])
    bars=["".join(x)+"-1" for x in __import__("itertools").islice(__import__("itertools").product("ACGT",repeat=16),50)]
    assert strict_orientation(["GeneID"]+bars,["SOAT1"]+["0"]*50,{"SOAT1"})[0]=="genes_by_cells"
    assert len(HCC_GENES)==8 and len(SC_ENDPOINTS)==5 and len(old.PANELS)==6
    assert prior("hcc_rna")["exact_sample_columns"]==70
    assert len(prior("sc_TC")["records"])==len(prior("sc_mLN")["records"])==10
    print("Five focused exact-column, orientation, fixed-family and archived source assertions passed.")
def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["bulk","sc_headers","sc_pilot","sc_full","hcc_rna","hcc_protein","selftest"]);p.add_argument("--out",type=Path);a=p.parse_args()
    if a.mode=="selftest":selftest();return
    try:z=globals()[a.mode]()
    except (OSError,EOFError,ValueError,TypeError,KeyError,IndexError,UnicodeError,StopIteration,zipfile.BadZipFile) as e:z={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:220]}
    z.update(mode=a.mode,scientific_status="HOLD_NOVELTY_AND_INDEPENDENT_VALIDATION")
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n");print(json.dumps(z,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
