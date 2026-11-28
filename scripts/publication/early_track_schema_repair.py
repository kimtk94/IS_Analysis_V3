#!/usr/bin/env python3
"""Repair source-derived schema, preserving #353 fixed hypotheses and scientific HOLD gates."""
import argparse,collections,csv,gzip,hashlib,io,json,math,re,statistics,sys,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"paired_helpers/scripts/publication"))
sys.path.insert(0,str(ROOT/"prior_helpers/scripts/publication"))
import early_track_paired_cohorts as old
base=old.base
require=old.require
# Only metadata-confirmed keys; never derive pairing from matrix order or approximate matching.
def characteristics(rows,key,n):
    pattern=re.compile(re.escape(key)+r":\s*(.+)",re.I)
    matches=[]
    for rr in rows.get("!Sample_characteristics_ch1",[]):
        if len(rr)==n and all(pattern.fullmatch(v.strip()) for v in rr):
            matches.append([pattern.fullmatch(v.strip())[1].strip() for v in rr])
    require(len(matches)==1,"missing or ambiguous characteristic "+key)
    return matches[0]
def normalize(v):return re.sub(r"\s+"," ",v.strip().lower())
def bulk_sample_keys(rows,ids,titles):
    n=len(ids);desc=old.row_field(rows,"!Sample_description",n);src=old.row_field(rows,"!Sample_source_name_ch1",n)
    tissue=characteristics(rows,"tissue",n)
    require(n==108 and len(set(desc))==108,"108 unique descriptions")
    keys=[]
    for title,d,s,t in zip(titles,desc,src,tissue):
        a=re.fullmatch(r"patient (\d{3}) (primary tumor sample|metastasis) (\d+)",title.strip(),re.I)
        b=re.fullmatch(r"P(\d{3})(Tm|Me)(\d+)",d.strip(),re.I)
        require(a is not None and b is not None,"unrecognized title/description convention")
        require(a[1]==b[1] and int(a[3])==int(b[3]),"patient or biopsy number mismatch")
        group="PT" if a[2].lower().startswith("primary") else "MET"
        require((b[2].lower()=="tm")==(group=="PT"),"title/description site mismatch")
        s,t=normalize(s),normalize(t)
        require(s==t,"source name and tissue characteristic disagree")
        # Official GSM8486852 "tumor of the tongue"; GSM8486850 "metastasis of the lymph node".
        require((bool(re.search(r"\b(?:tumou?r|carcinoma)\b",t)) and not re.search(r"metasta|lymph node|non.tumou?r|normal",t)) if group=="PT" else bool(re.search(r"metasta",t) and re.search(r"lymph",t)),"unrecognized site tissue: "+t[:80])
        keys.append((a[1],group))
    require(len({p for p,g in keys})==72 and collections.Counter(g for p,g in keys)=={"PT":87,"MET":21},"publisher totals not reproduced")
    return desc,keys
def hcc_keys():
    rows,ids,titles,meta=base.geo_metadata("GSE124535")
    n=len(ids);require(n==70 and len(set(titles))==70,"70 unique HCC titles")
    source=old.row_field(rows,"!Sample_source_name_ch1",n)
    tissue=characteristics(rows,"tissue",n)
    require(source==titles,"GEO source names do not equal author sample titles")
    require(all(re.fullmatch(r"L\d{3}[PT]",t) for t in titles),"invalid title")
    require(all(normalize(v)==("non-tumor" if t[-1]=="P" else "hcc") for t,v in zip(titles,tissue)),"tissue/suffix mismatch")
    require(len({t[:-1] for t in titles})==35 and all(t[:-1]+("T" if t[-1]=="P" else "P") in titles for t in titles),"35 exact pairs missing")
    aliases=[rr for rr in rows.get("!Sample_description",[]) if len(rr)==n and all(re.fullmatch(r"\d+[PT]",v.strip()) for v in rr)]
    require(len(aliases)==1,"no unique alternate RNA identifier field")
    alias=[v.strip() for v in aliases[0]]
    require(len(set(alias))==n and all(t[-1]==a[-1] for t,a in zip(titles,alias)),"alias duplicate or site conflict")
    require(len({a[:-1] for a in alias})==35,"alias patient cardinality")
    by=collections.defaultdict(set)
    for t,a in zip(titles,alias):by[t[:-1]].add(a[:-1])
    require(all(len(v)==1 for v in by.values()) and len({next(iter(v)) for v in by.values()})==35,"alias patient mapping conflict")
    return rows,ids,titles,alias,meta
def hcc_rna():
    rows,ids,titles,aliases,meta=hcc_keys()
    raw,src=base.get(base.geo_folder("GSE124535")+"suppl/GSE124535_HCC.RNA-seq.35.samples.fpkm.txt.gz",25_000_000)
    rr=csv.reader(io.TextIOWrapper(gzip.GzipFile(fileobj=io.BytesIO(raw)),encoding="utf-8-sig"),delimiter="\t")
    header=next(rr)
    candidates=[]
    for label,keys in (("title",titles),("GEO_description_alias",aliases),("GSM",ids)):
        hits=[i for i,v in enumerate(header) if v in set(keys)]
        if len(hits)==70 and len({header[i] for i in hits})==70:
            candidates.append((label,hits,[keys.index(header[i]) for i in hits]))
    out={"metadata":meta,"source":src,"title_pairs":35,"alias_pairs":35,"header_width":len(header),"candidate_column_matches":{k:sum(v in set(q) for v in header) for k,q in (("title",titles),("description_alias",aliases),("GSM",ids))}}
    if len(candidates)!=1:
        return {"status":"HOLD_EXACT_RNA_HEADER",**out,"reason":"no unique exact 70-column identity; no alias normalization or inferred order"}
    label,hits,order=candidates[0];valid=seen=other=0
    for row in rr:
        seen+=1;require(len(row)==len(header),"FPKM row width mismatch")
        try:
            v=[float(row[i]) for i in hits]
            if all(math.isfinite(q) and q>=0 for q in v):valid+=1
            else:other+=1
        except ValueError:other+=1
    require(seen>=10000 and valid>=10000,"insufficient expression matrix")
    return {"status":"HCC_RNA_EXACT_SAMPLE_MATRIX","column_mapping":label,**out,"exact_sample_columns":70,"feature_rows":seen,"finite_nonnegative_rows":valid,"other_rows":other,
            "limits":["FPKM only; no raw-count DE. Exact title/description alias validated by GEO; still need independent protein aliquot linkage. No patient rows exported."]}
def hcc_protein():
    import openpyxl
    rows,ids,titles,aliases,meta=hcc_keys()
    out={"metadata":meta,"title_pairs":35,"alias_pairs":35,"pride":{}}
    # Official PRIDE directories confirmed in ProteomeXchange; raw names are NOT biological matching.
    for acc,month in (("PXD006512","02"),("PXD008373","03")):
        try:
            raw,src=base.get("https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/"+month+"/"+acc+"/",4_000_000)
            names=base.links(raw);archives=[n for n in names if re.fullmatch(r"[A-Za-z0-9_.-]+\.(?:zip|txt|tsv|xlsx)",n,re.I) and re.search(r"maxquant|proteinGroups|quant|result|sdrf",n,re.I)]
            out["pride"][acc]={"source":src,"raw_file_links":sum(n.lower().endswith(".raw") for n in names),"processed_candidates":archives}
        except (OSError,ValueError,KeyError) as e:out["pride"][acc]={"status":"HOLD_DIRECTORY","reason":str(e)[:150]}
    url="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-019-0987-8/MediaObjects/41586_2019_987_MOESM4_ESM.zip"
    out["publisher_supplement"]={"status":"HOLD_SOURCE","url":url}
    try:
        raw,src=base.get(url,100_000_000)
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            infos=z.infolist();require(len(infos)<=200 and sum(v.file_size for v in infos)<=600_000_000,"archive expansion cap")
            sheets=[];nwb=0
            for f in infos:
                if not f.filename.lower().endswith(".xlsx") or f.file_size>80_000_000:continue
                nwb+=1;wb=openpyxl.load_workbook(io.BytesIO(z.read(f)),read_only=True,data_only=True)
                for si,ws in enumerate(wb,1):
                    best={k:0 for k in ("title","description_alias")};pairs=best.copy()
                    for row in ws.iter_rows(min_row=1,max_row=min(8,ws.max_row or 1),values_only=True):
                        for k,q in (("title",titles),("description_alias",aliases)):
                            vals=[v for v in row if isinstance(v,str) and v in set(q)]
                            if len(vals)>best[k] and len(vals)==len(set(vals)):
                                best[k]=len(vals);pairs[k]=sum(p+"P" in vals and p+"T" in vals for p in {t[:-1] for t in vals})
                    if max(best.values()) or (ws.max_column or 0)>=70:
                        sheets.append({"table_file":Path(f.filename).name,"sheet_index":si,"dimensions":[ws.max_row,ws.max_column],"exact_sample_keys":best,"complete_pairs":pairs})
                wb.close()
        out["publisher_supplement"]={"status":"PUBLIC_ZIP_SCHEMA_READ","source":src,"members":len(infos),"xlsx_workbooks":nwb,"candidate_sheets":sheets}
    except (OSError,ValueError,KeyError,zipfile.BadZipFile,TypeError) as e:out["publisher_supplement"]["reason"]=str(e)[:160]
    return {"status":"HCC_PROTEOME_LINKAGE_QUALIFICATION",**out,"limits":["RNA title/alias relationship from GEO only. Candidate protein column labels are not confirmed biological aliquot or quantitative units. No raw or individual data outputs."]}
def barcode(v):
    # 10x barcode with optional publisher sample prefix/suffix; no arbitrary strings.
    return bool(re.fullmatch(r"(?:[A-Za-z0-9_.-]+[_:])?[ACGT]{16}(?:-\d+)?(?:[_:.][A-Za-z0-9.-]+)?",v))
def strip_field(v):return v.strip().strip('"')
def identify(header,first):
    h=[strip_field(v) for v in header];a=[strip_field(v) for v in first]
    require(len(a)>=51,"short data row")
    n=len(a)-1
    if len(h)==n:features=h
    elif len(h)==n+1 and h[0].lower() in {"","gene","gene_id","gene_name","genes","symbol","barcode","cell","cell_id","cell barcode","x"}:features=h[1:]
    else:raise ValueError("header/row width mismatch")
    require(len(set(features))==n and all(features),"nonunique or empty column keys")
    if all(barcode(v) for v in features) and not barcode(a[0]):orientation="genes_by_cells"
    elif barcode(a[0]) and old.SC_GENES<=set(features):orientation="cells_by_genes"
    else:raise ValueError("unresolved orientation; no barcode/gene identity assumption")
    return orientation,features,a
def head_rows(raw):
    # gzip.GzipFile handles concatenated members (zlib single-member prefix did not).
    with gzip.GzipFile(fileobj=io.BytesIO(raw)) as f:
        a=f.readline(16_000_001);b=f.readline(16_000_001)
        require(a.endswith(b"\n") and b.endswith(b"\n") and len(a)<=16_000_000 and len(b)<=16_000_000,"incomplete first two lines")
    return [next(csv.reader([v.decode("utf-8-sig").rstrip("\r\n")],delimiter="\t")) for v in (a,b)]
def sc_meta(site):
    keys,files,meta=old.sc_mapping()
    out=[]
    for p,g in sorted(keys):
        if g!=site:continue
        _,sid,_=keys[p,g];z={"site":g,"status":"HOLD_HEADER"}
        try:
            require(sid in files,"missing official sample file")
            raw,src=base.get(files[sid],2_097_152,(0,2_097_151))
            z.update(source=src)
            h,a=head_rows(raw)
            z.update(header_width=len(h),first_row_width=len(a))
            mode,features,_=identify(h,a)
            z.update(status="DENSE_HEADER_VERIFIED",orientation=mode,columns=len(features))
        except (OSError,EOFError,ValueError,KeyError,UnicodeError) as e:z["reason"]=str(e)[:150]
        out.append(z);time.sleep(1)
    return {"status":"ORAL_PAIRED_SOURCE_HEADER_REPAIR","site":site,"metadata":meta,"target_matrices":10,"verified":sum(z["status"]=="DENSE_HEADER_VERIFIED" for z in out),"records":out,
            "limits":["Header orientation only; no patient inference. Exact barcode and fixed marker identifiers required. No cell identifiers exported."]}
def numeric(fields,n):
    import numpy as np
    require(len(fields)==n,"expression row width changed")
    # numpy cannot silently accept malformed values: token count checked and numeric conversion strict.
    v=np.array(fields,dtype=float)
    require(np.isfinite(v).all() and (v>=0).all(),"nonfinite or negative expression")
    return v
def score_cells(markers,tot,det,mt):
    import numpy as np
    require(set(markers)==old.SC_GENES,"fixed marker set incomplete")
    good=(det>=200)&(tot>=500)&(mt<=.2*tot)
    require(int(good.sum())>0,"no QC cells")
    ans={}
    for th in (2,3):
        lineage=np.vstack([sum((markers[g]>0).astype(int) for g in gs)>=th for gs in old.SC_MARKERS.values()])
        n=0;by={}
        for j,k in enumerate(old.SC_MARKERS):
            mask=good&(lineage.sum(axis=0)==1)&lineage[j]
            if k in {"epithelial","fibroblast","endothelial"}:mask&=markers["PTPRC"]==0
            nc=int(mask.sum());n+=nc
            z={"cells":nc,"status":"HOLD_FEWER_THAN_50_CELLS"}
            if nc>=50:
                den=float(tot[mask].sum());genes=old.PANELS["pEMT"]+["SPP1","CXCR4"] if k=="epithelial" else ["SPP1","CXCR4"]
                cpm={g:float(markers[g][mask].sum()*1e6/den) for g in genes}
                z.update(status="DESCRIPTIVE_MARKER_GATE",cpm=cpm)
                if k=="epithelial":z["pemt_mean_log2_cpm_plus1"]=statistics.mean(math.log2(cpm[g]+1) for g in old.PANELS["pEMT"])
            by[k]=z
        ans[str(th)]={"assigned":n,"compartments":by}
    return {"qc_cells":int(good.sum()),"thresholds":ans}
def sc_counts(raw):
    import itertools,numpy as np
    with gzip.GzipFile(fileobj=io.BytesIO(raw)) as f:
        rr=csv.reader(io.TextIOWrapper(f,encoding="utf-8-sig"),delimiter="\t")
        header=next(rr);first=next(rr);mode,cols,a=identify(header,first);n=len(cols)
        seen=set()
        if mode=="genes_by_cells":
            totals=np.zeros(n);det=np.zeros(n,dtype=np.int32);mt=np.zeros(n);markers={}
            for row in itertools.chain([first],rr):
                k=strip_field(row[0]);require(k and k not in seen,"duplicate/empty feature");seen.add(k)
                v=numeric(row[1:],n);totals+=v;det+=v>0
                if k.upper().startswith("MT-"):mt+=v
                if k in old.SC_GENES:markers[k]=v
            require(len(seen)>=10000 and any(k.upper().startswith("MT-") for k in seen),"feature universe incomplete")
            nf,nc=len(seen),n
        else:
            require(n>=10000,"feature universe incomplete")
            indices={g:cols.index(g) for g in old.SC_GENES};mi=np.array([g.upper().startswith("MT-") for g in cols]);require(mi.any(),"no mitochondrial genes")
            totals=[];det=[];mt=[];markers={g:[] for g in indices}
            for row in itertools.chain([first],rr):
                k=strip_field(row[0]);require(barcode(k) and k not in seen,"duplicate or invalid barcode");seen.add(k)
                v=numeric(row[1:],n);totals.append(v.sum());det.append(np.count_nonzero(v));mt.append(v[mi].sum())
                for g,i in indices.items():markers[g].append(v[i])
            totals,det,mt=map(np.asarray,(totals,det,mt));markers={k:np.asarray(v) for k,v in markers.items()}
            nf,nc=n,len(seen)
        require(nc>=50,"cell count too low")
        return {"orientation":mode,"features":nf,"cells":nc,**score_cells(markers,totals,det,mt)}
def sc_pilot():
    keys,files,meta=old.sc_mapping();out=[]
    for g in ("TC","mLN"):
        _,sid,_=keys["P01",g];z={"site":g,"status":"HOLD_PILOT"}
        try:
            require(sid in files,"missing exact official source")
            raw,src=base.get(files[sid],300_000_000)
            z.update(source=src)
            z.update(status="FULL_PAIRED_SPECIMEN_PILOT",**sc_counts(raw))
        except (OSError,EOFError,ValueError,KeyError,UnicodeError,StopIteration) as e:z.update(status="HOLD_PILOT",reason=str(e)[:160])
        out.append(z)
    return {"status":"ORAL_FIRST_PAIR_SCHEMA_REPAIR","metadata":meta,"records":out,"limits":["Preselected P01 only; no new inference. Cap increased because original GEO TC file exceeds 220 MB; no combined 5-7 GB download. No cell-level exports; no malignant label."]}
def selftest():
    assert characteristics({"!Sample_characteristics_ch1":[["tissue: HCC","tissue: non-tumor"]]},"tissue",2)==["HCC","non-tumor"]
    import itertools
    bars=["".join(v)+"-1" for v in itertools.islice(itertools.product("ACGT",repeat=16),50)]
    raw=gzip.compress(("\t".join([""]+bars)+"\n").encode())+gzip.compress(("\t".join(["GENE"]+["0"]*50)+"\n").encode())
    h,a=head_rows(raw);assert identify(h,a)[0]=="genes_by_cells"
    genes=sorted(old.SC_GENES)+["GENE%d"%i for i in range(50)]
    assert identify([""]+genes,["ACGTACGTACGTACGT-1"]+["0"]*len(genes))[0]=="cells_by_genes"
    assert old.sign_summary([1]*5)["p_exact_sign"]==.0625 and len(old.PANELS)==6
    arch=json.loads((ROOT/"docs/publication/evidence/early_track_paired_readback_20261125.json").read_text())
    assert arch["run_id"]==45601731870 and len(arch["records"])==5
    print("Five source schema, gzip member, orientation, fixed-panel and archive checks passed.")
def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["bulk","sc_TC","sc_mLN","sc_pilot","hcc_rna","hcc_protein","selftest"]);p.add_argument("--out",type=Path);a=p.parse_args()
    if a.mode=="selftest":selftest();return
    try:
        if a.mode=="bulk":
            old.bulk_sample_keys=bulk_sample_keys;x=old.bulk()
        elif a.mode in {"sc_TC","sc_mLN"}:x=sc_meta(a.mode[3:])
        else:x=globals()[a.mode]()
    except (OSError,EOFError,ValueError,KeyError,IndexError,TypeError,UnicodeError,StopIteration) as e:x={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:220]}
    x.update(mode=a.mode,scientific_status="HOLD_INDEPENDENT_VALIDATION_AND_NOVELTY")
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps(x,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
