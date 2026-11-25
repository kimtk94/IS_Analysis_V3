#!/usr/bin/env python3
"""New independent-cohort gates; fixed paired bulk test, scRNA pilot, HCC exact source keys."""
import argparse,collections,csv,gzip,hashlib,io,json,math,re,statistics,sys,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"prior_helpers/scripts/publication"))
import early_track_acceleration as base
base.HOSTS.update({"ftp.ensembl.org","ftp.pride.ebi.ac.uk","media.springernature.com"})
PANELS={"SPP1":["SPP1"],"CXCR4":["CXCR4"],"pEMT":["SERPINE1","TGFBI","MMP10","LAMC2","P4HA2","PDPN","ITGA5","LAMA3","CDH13","TNC","MMP2","EMP3","INHBA","LAMB3","VIM"],
        "epithelial":["EPCAM","KRT5","KRT14","KRT19"],"T_cell":["CD3D","CD3E","TRAC"],"B_cell":["MS4A1","CD79A","CD74"]}
GENES=set(sum(PANELS.values(),[]))
SC_MARKERS={"epithelial":["EPCAM","KRT5","KRT14","KRT8","KRT18","KRT19"],"myeloid":["LST1","TYROBP","FCER1G","CD68"],
"T_cell":["CD3D","CD3E","TRAC"],"B_cell":["MS4A1","CD79A","CD79B","CD19"],
"fibroblast":["COL1A1","COL1A2","DCN","LUM"],"endothelial":["PECAM1","VWF","KDR","EMCN"]}
SC_GENES=GENES|set(sum(SC_MARKERS.values(),[]))|{"PTPRC"}
def require(x,msg):base.require(x,msg)
def strip_version(x):return re.sub(r"\.\d+$","",x)
def row_field(rows,key,n):
    v=rows.get(key,[])
    require(len(v)==1 and len(v[0])==n,"missing or ambiguous "+key)
    return v[0]
def exact_header_match(header,fields):
    # Every data column maps uniquely to one GEO sample via an explicit publisher field.
    options=[]
    for key,vals in fields.items():
        if len(vals)==len(header) and len(set(vals))==len(vals) and set(vals)==set(header):
            options.append((key,[vals.index(v) for v in header]))
    require(options,"count columns do not exactly match any complete unique GEO identifier field")
    require(all(v[1]==options[0][1] for v in options),"identifier mappings conflict")
    return options[0]
def bulk_sample_keys(rows,ids,titles):
    n=len(ids); desc=row_field(rows,"!Sample_description",n); tissue=row_field(rows,"!Sample_source_name_ch1",n)
    require(n==108 and len(set(desc))==108,"expected 108 GEO samples and unique descriptions")
    keys=[]
    for t,d,s in zip(titles,desc,tissue):
        a=re.fullmatch(r"patient (\d{3}) (primary tumor sample|metastasis) \d+",t.strip(),re.I)
        b=re.fullmatch(r"P(\d{3})(Me|Tu|Pr|T|M|PT|Met|Tumor)\d+",d.strip(),re.I)
        # The explicit patient title and GEO tissue establish identity; description must independently encode same patient.
        require(a is not None,"unrecognized bulk sample title")
        require(re.match(r"^P"+a[1]+r"[^0-9]",d,re.I) is not None,"description patient mismatch")
        group="PT" if a[2].lower().startswith("primary") else "MET"
        s=s.lower()
        require(("primary" in s and "tumor" in s and "metasta" not in s) if group=="PT" else ("metasta" in s and "lymph" in s),"GEO tissue/title conflict")
        if b:require((b[2].lower() in {"tu","pr","t","pt","tumor"})==(group=="PT"),"description group conflict")
        keys.append((a[1],group))
    require(len({p for p,g in keys})==72 and collections.Counter(g for p,g in keys)=={"PT":87,"MET":21},"publisher totals not reproduced")
    return desc,keys
def gtf_mapping():
    # Use exactly the gene annotation release named by original GEO, not today's HGNC aliases.
    url="https://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"
    raw,src=base.get(url,100_000_000)
    symbols=collections.defaultdict(set)
    for line in gzip.GzipFile(fileobj=io.BytesIO(raw)):
        if line.startswith(b"#"):continue
        a=line.decode().rstrip("\n").split("\t")
        if len(a)!=9 or a[2]!="gene":continue
        x=dict(re.findall(r'(\w+) "([^"]+)"',a[8]))
        if x.get("gene_name") in GENES and re.fullmatch(r"ENSG\d+",x.get("gene_id","")):
            symbols[x["gene_name"]].add(x["gene_id"])
    return {g:next(iter(v)) for g,v in symbols.items() if len(v)==1},src
def sign_summary(delta):
    from scipy.stats import binomtest
    n=len(delta); nonzero=[x for x in delta if x!=0];p=sum(x>0 for x in nonzero)
    if n<5:return {"status":"HOLD_FEWER_THAN_5_PAIRED_PATIENTS","pairs":n}
    med=statistics.median(delta);loo=[statistics.median(delta[:i]+delta[i+1:]) for i in range(n)]
    return {"status":"EXPLORATORY_PAIRED_WHOLE_TISSUE","pairs":n,"positive":p,"negative":len(nonzero)-p,"ties":n-len(nonzero),
            "median_delta_log2_cpm_plus1":med,"mean_delta_log2_cpm_plus1":statistics.mean(delta),
            "loo_median_range":[min(loo),max(loo)],"p_exact_sign":float(binomtest(p,len(nonzero),.5).pvalue) if nonzero else 1.0}
def bulk():
    import numpy as np
    rows,ids,titles,meta=base.geo_metadata("GSE275870")
    desc,keys=bulk_sample_keys(rows,ids,titles)
    url=base.geo_folder("GSE275870")+"suppl/GSE275870_counts.csv.gz"
    raw,src=base.get(url,15_000_000)
    rr=csv.reader(io.TextIOWrapper(gzip.GzipFile(fileobj=io.BytesIO(raw)),encoding="utf-8-sig"))
    header=next(rr)
    require(len(header)==109 and len(set(header[1:]))==108,"expected one feature and 108 count columns")
    ident,order=exact_header_match(header[1:],{"GSM":ids,"description":desc,"title":titles})
    totals=np.zeros(108,dtype=np.float64);seen=set();vectors={};ensg=0
    for row in rr:
        require(len(row)==109,"invalid count width")
        key=strip_version(row[0].strip())
        require(bool(key) and key not in seen,"duplicate normalized gene identifier")
        seen.add(key);ensg+=bool(re.fullmatch(r"ENSG\d+",key))
        v=np.asarray(row[1:],dtype=float)
        require(np.isfinite(v).all() and (v>=0).all() and (v==np.floor(v)).all(),"noninteger/negative raw counts")
        totals+=v;vectors[key]=v
    require(len(seen)>=10000 and (totals>0).all(),"count matrix incomplete")
    if ensg>=.9*len(seen):
        mapping,ann=gtf_mapping();mode="Ensembl_GRCh38_release_100"
    else:
        mapping={g:g for g in GENES if g in seen};ann=None;mode="exact_unique_symbol"
    buckets=collections.defaultdict(list)
    for j,k in enumerate(order):buckets[keys[k]].append(j)
    patients=sorted({p for p,g in buckets if (p,"PT") in buckets and (p,"MET") in buckets})
    out={}
    for name,genes in PANELS.items():
        usable=[g for g in genes if g in mapping and mapping[g] in vectors]
        if len(usable)!=len(genes):
            out[name]={"status":"HOLD_INCOMPLETE_FIXED_PANEL","mapped":len(usable),"required":len(genes),"missing":sorted(set(genes)-set(usable))};continue
        score=np.mean([np.log2(vectors[mapping[g]]/totals*1e6+1) for g in genes],axis=0)
        delta=[float(np.mean(score[buckets[p,"MET"]])-np.mean(score[buckets[p,"PT"]])) for p in patients]
        out[name]={**sign_summary(delta),"mapped":len(usable),"required":len(genes)}
    # Fixed six endpoints, including untestable endpoints in Bonferroni denominator.
    for z in out.values():
        if "p_exact_sign" in z:z["bonferroni_6"]=min(1,6*z["p_exact_sign"])
    return {"status":"ORAL_PAIRED_BULK_FIXED_PANEL","source":src,"metadata":meta,"mapping":mode,"annotation":ann,
            "gene_rows":len(seen),"samples":len(ids),"patients":len({p for p,g in keys}),"paired_patients":len(patients),
            "column_mapping":ident,"group_counts":dict(collections.Counter(g for p,g in keys)),"results":out,
            "limits":["Original Pyko et al. DOI 10.1002/cac2.12660 already compares paired metastases. Fixed context panels from earlier studies, not novel discoveries. Mean sample log2(CPM+1) per patient/site, equal patient weight, exact sign test and fixed family 6. Tissue composition unresolved; no malignant-cell or causal inference. No individual-level output."]}
def sc_mapping():
    rows,ids,titles,meta=base.geo_metadata("GSE198315")
    require(len(ids)==36,"expected 36 specimens")
    tissue=row_field(rows,"!Sample_source_name_ch1",36)
    # Actual official title encodes P01..P10 and tissue, cross-check if tissue field is informative.
    keys={}
    for i,t in enumerate(titles):
        m=re.fullmatch(r"(P\d{2})-(TC|TP|NT|mLN)",t)
        require(m is not None and (m[1],m[2]) not in keys,"ambiguous sample title")
        keys[m[1],m[2]]=(i,ids[i],t)
    patients=sorted({p for p,g in keys})
    require(patients==["P%02d"%i for i in range(1,11)] and all((p,"TC") in keys and (p,"mLN") in keys for p in patients),"ten TC/mLN pairs not established")
    # Per-sample official GEO text, from metadata only (no inferred path from titles).
    files={}
    for k,rrs in rows.items():
        if not k.startswith("!Sample_supplementary_file"):continue
        for rr in rrs:
            require(len(rr)==36,"sample file row misalignment")
            for i,(sid,v) in enumerate(zip(ids,rr)):
                u=v.replace("ftp://","https://",1)
                if re.fullmatch(r"https://ftp\.ncbi\.nlm\.nih\.gov/geo/samples/GSM\d+nnn/"+sid+r"/suppl/"+sid+r"_"+re.escape(titles[i])+r"_UMI_count_matrix\.txt\.gz",u):
                    require(sid not in files,"multiple sample matrices")
                    files[sid]=u
    return keys,files,meta
def parse_dense_header(raw):
    txt=base.gzip_prefix(raw,32_000_000);lines=txt.splitlines()
    require(len(lines)>=2,"dense source header truncated")
    header=lines[0].split("\t");first=lines[1].split("\t")
    # Explicit leading blank/name or R table header with one fewer field.
    if len(first)==len(header)+1:bars=header
    elif len(first)==len(header) and header[0].strip('"').lower() in {"","gene","gene_id","gene_name","genes","symbol"}:bars=header[1:]
    else:raise ValueError("dense header width not explicit")
    require(len(bars)>=50 and len(set(bars))==len(bars) and all(v.strip('"') for v in bars),"nonunique cell columns")
    require(len(first)==len(bars)+1,"first expression row width")
    return len(bars)
def sc_meta():
    keys,files,meta=sc_mapping();out=[]
    for p,g in sorted(keys):
        if g not in {"TC","mLN"}:continue
        _,sid,_=keys[p,g]
        z={"site":g,"status":"HOLD_FILE"}
        try:
            require(sid in files,"exact per-sample GEO URL not in metadata")
            raw,src=base.get(files[sid],524288,(0,524287))
            z.update(status="DENSE_HEADER_VERIFIED",cells=parse_dense_header(raw),source=src)
        except (OSError,ValueError,KeyError,UnicodeError) as e:z["reason"]=str(e)[:140]
        out.append(z);time.sleep(.5)
    return {"status":"ORAL_SINGLE_CELL_PAIR_SOURCE_GATE","metadata":meta,"patients":10,"total_samples":len(keys),
            "site_counts":dict(collections.Counter(g for p,g in keys)),"target_matrices":20,"header_verified":sum(z["status"]=="DENSE_HEADER_VERIFIED" for z in out),"records":out,
            "limits":["Exact author-coded patient and tissue title, not clinical records. Original GSE198315 DOI 10.1016/j.isci.2025.112988 already analyzes pEMT and immune states; no new novelty claim. Header only, no expression-body or cell-type validation."]}
def sc_counts(raw):
    import numpy as np
    stream=io.TextIOWrapper(gzip.GzipFile(fileobj=io.BytesIO(raw)),encoding="utf8")
    header=next(stream).rstrip("\r\n").split("\t");first=next(stream).rstrip("\r\n").split("\t")
    n=len(first)-1
    require(n>=50 and (len(header)==n or (len(header)==n+1 and header[0].strip('"').lower() in {"","gene","gene_id","gene_name","genes","symbol"})),"matrix header contract")
    bars=header[-n:];require(len(set(bars))==n,"duplicate cell columns")
    import itertools
    totals=np.zeros(n,dtype=float);det=np.zeros(n,dtype=np.int32);mt=np.zeros(n);seen=set();markers={}
    for line in itertools.chain(["\t".join(first)],stream):
        key,sep,rest=line.rstrip("\r\n").partition("\t");key=key.strip('"')
        require(sep and key and key not in seen,"duplicate/missing gene key")
        seen.add(key)
        v=np.fromstring(rest,sep="\t")
        require(rest.count("\t")==n-1 and len(v)==n and np.isfinite(v).all() and (v>=0).all(),"invalid expression row")
        # SoupX-derived data may be fractional. Do not coerce to integers.
        totals+=v;det+=(v>0)
        if key.upper().startswith("MT-"):mt+=v
        if key in SC_GENES:markers[key]=v
    require(len(seen)>=10000 and set(markers)==SC_GENES and any(g.upper().startswith("MT-") for g in seen),"gene universe / fixed markers incomplete")
    good=(det>=200)&(totals>=500)&(mt<=.2*totals)
    require(int(good.sum())>0,"no QC cells")
    result={}
    for thresh in (2,3):
        lineage=np.vstack([sum((markers[g]>0).astype(int) for g in genes)>=thresh for genes in SC_MARKERS.values()])
        assigned=lineage.sum(axis=0)
        compartments={}
        for j,key in enumerate(SC_MARKERS):
            mask=good&(assigned==1)&lineage[j]
            if key in {"epithelial","fibroblast","endothelial"}:mask&=(markers["PTPRC"]==0)
            z={"cells":int(mask.sum()),"status":"HOLD_FEWER_THAN_50_CELLS"}
            if mask.sum()>=50:
                den=float(totals[mask].sum());genes=PANELS["pEMT"]+["SPP1","CXCR4"] if key=="epithelial" else ["SPP1","CXCR4"]
                cpm={g:float(markers[g][mask].sum()*1e6/den) for g in genes}
                z.update(status="DESCRIPTIVE_MARKER_GATE",cpm=cpm)
                if key=="epithelial":z["pemt_mean_log2_cpm_plus1"]=statistics.mean(math.log2(cpm[g]+1) for g in PANELS["pEMT"])
            compartments[key]=z
        result[str(thresh)]=compartments
    return {"features":len(seen),"cells":n,"qc_cells":int(good.sum()),"thresholds":result}
def sc_pilot():
    keys,files,meta=sc_mapping();out=[]
    # Outcome-blind first author-coded pair only for complete file/analysis feasibility, no inference.
    for g in ("TC","mLN"):
        _,sid,_=keys["P01",g];z={"site":g,"status":"HOLD_PILOT"}
        try:
            require(sid in files,"no exact sample file")
            raw,src=base.get(files[sid],220_000_000)
            z.update(status="FULL_PAIRED_SPECIMEN_PILOT",source=src,**sc_counts(raw))
        except (OSError,ValueError,KeyError,UnicodeError) as e:z["reason"]=str(e)[:150]
        out.append(z)
    return {"status":"ORAL_SINGLE_CELL_FIRST_PAIR_PILOT","metadata":meta,"records":out,"limits":["Exactly first author-coded pair, not representative evidence. No test or cell-level export. Original SoupX processed values, marker gates only, not malignant cell annotations or proof of independence."]}
def hcc_keys():
    rows,ids,titles,meta=base.geo_metadata("GSE124535")
    require(len(ids)==70 and len(set(titles))==70,"expected 70 unique HCC titles")
    keys=[]
    for t in titles:
        require(re.fullmatch(r"L\d{3}[PT]",t) is not None,"HCC title not exact LxxxP/T")
        keys.append((t[:-1],t[-1]))
    require(len({p for p,g in keys})==35 and all((p,"P") in keys and (p,"T") in keys for p,g in keys),"35 exact title pairs not established")
    # Confirm two source tissue classes and consistent correspondence across all titles.
    tissue=row_field(rows,"!Sample_source_name_ch1",70)
    groups={g:set(tissue[i] for i,k in enumerate(keys) if k[1]==g) for g in ("P","T")}
    require(len(groups["P"])==len(groups["T"])==1 and groups["P"]!=groups["T"],"tissue classes not consistent")
    return rows,ids,titles,meta
def hcc_rna():
    rows,ids,titles,meta=hcc_keys()
    raw,src=base.get(base.geo_folder("GSE124535")+"suppl/GSE124535_HCC.RNA-seq.35.samples.fpkm.txt.gz",25_000_000)
    stream=io.TextIOWrapper(gzip.GzipFile(fileobj=io.BytesIO(raw)),encoding="utf-8-sig")
    header=next(stream).rstrip("\r\n").split("\t")
    # Accept only the explicit GEO sample title/ID set, never partial or column-order matching.
    hits=[i for i,v in enumerate(header) if v in set(titles)]
    require(len(hits)==70 and len({header[i] for i in hits})==70,"RNA header does not contain all 70 exact publisher sample keys")
    seen=0;valid=0;nonfinite=0
    for line in stream:
        a=line.rstrip("\r\n").split("\t");seen+=1
        require(len(a)==len(header),"RNA row width mismatch")
        try:
            v=[float(a[i]) for i in hits]
            if all(math.isfinite(q) and q>=0 for q in v):valid+=1
            else:nonfinite+=1
        except ValueError:nonfinite+=1
    return {"status":"HCC_RNA_EXACT_SAMPLE_MATRIX","metadata":meta,"source":src,"exact_sample_columns":70,"title_pairs":35,"feature_rows":seen,"finite_nonnegative_rows":valid,"other_rows":nonfinite,
            "limits":["FPKM as published; no count-based differential expression. No patient data exported; still no protein aliquot linkage."]}
def hcc_protein():
    import openpyxl
    rows,ids,titles,meta=hcc_keys();title_set=set(titles)
    url="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-019-0987-8/MediaObjects/41586_2019_987_MOESM4_ESM.zip"
    out={"metadata":meta,"source_zip_url":url,"pride":{},"publisher_supplement":{"status":"HOLD_SOURCE"}}
    for acc,month in (("PXD006512","02"),("PXD008373","03")):
        try:
            raw,src=base.get("https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/"+month+"/"+acc+"/",4_000_000)
            names=base.links(raw);archives=[n for n in names if re.fullmatch(r"[A-Za-z0-9_.-]+\.(?:zip|txt|tsv|xlsx)",n,re.I) and re.search(r"maxquant|proteinGroups|quant|result|sdrf",n,re.I)]
            out["pride"][acc]={"source":src,"raw_file_links":sum(n.lower().endswith(".raw") for n in names),"processed_candidates":archives}
        except (OSError,ValueError,KeyError) as e:out["pride"][acc]={"status":"HOLD_DIRECTORY","reason":str(e)[:150]}
    try:
        raw,src=base.get(url,100_000_000)
        z=zipfile.ZipFile(io.BytesIO(raw));infos=z.infolist()
        require(len(infos)<=200 and sum(v.file_size for v in infos)<=600_000_000,"supplement expansion cap")
        sheets=[];workbooks=0
        for f in infos:
            if not f.filename.lower().endswith(".xlsx") or f.file_size>80_000_000:continue
            workbooks+=1
            wb=openpyxl.load_workbook(io.BytesIO(z.read(f)),read_only=True,data_only=True)
            for si,ws in enumerate(wb,1):
                best=0;pair_n=0;group_types=[]
                for row in ws.iter_rows(min_row=1,max_row=min(8,ws.max_row or 1),values_only=True):
                    # Count exact RNA sample column keys only; no patient aliases or matching raw-file names.
                    vals=[v for v in row if isinstance(v,str) and v in title_set]
                    if len(vals)>best and len(vals)==len(set(vals)):
                        best=len(vals);pair_n=sum(p+"P" in vals and p+"T" in vals for p in {t[:-1] for t in vals})
                        group_types=sorted({t[-1] for t in vals})
                if best or (ws.max_column or 0)>=70:
                    sheets.append({"table_file":Path(f.filename).name,"sheet_index":si,"rows":ws.max_row,"columns":ws.max_column,
                                   "exact_RNA_title_columns":best,"complete_title_pairs":pair_n,"site_suffix_types":group_types})
            wb.close()
        out["publisher_supplement"]={"status":"PUBLIC_ZIP_SCHEMA_READ","source":src,"members":len(infos),"xlsx_workbooks":workbooks,"candidate_sheets":sheets}
    except (OSError,ValueError,KeyError,zipfile.BadZipFile,TypeError) as e:out["publisher_supplement"]["reason"]=str(e)[:160]
    return {"status":"HCC_PROTEOME_LINKAGE_QUALIFICATION",**out,"limits":["Exact title matches establish candidate columns, not protein quantitative units, biological aliquot identity or cohort independence. PRIDE raw filename overlap does not prove RNA/protein linkage. No restricted human genetic data. No raw files or clinical values exported."]}
def selftest():
    assert exact_header_match(["b","a"],{"gsm":["a","b"]})==("gsm",[1,0])
    try:exact_header_match(["x","b"],{"gsm":["a","b"]})
    except ValueError:pass
    else:raise AssertionError("invalid join accepted")
    assert strip_version("ENSG00001.17")=="ENSG00001" and len(PANELS)==6 and len(PANELS["pEMT"])==15
    assert sign_summary([1,1,1,1,1])["p_exact_sign"]==.0625
    assert sign_summary([1,-1,1])["status"].startswith("HOLD")
    assert len(json.loads((ROOT/"docs/publication/evidence/early_track_specificity_readback_20261125.json").read_text())["records"])==4
    print("Six focused join, panel, sign-test and archive checks passed.")
def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["bulk","sc_meta","sc_pilot","hcc_rna","hcc_protein","selftest"]);p.add_argument("--out",type=Path);a=p.parse_args()
    if a.mode=="selftest":selftest();return
    try:z=globals()[a.mode]()
    except (OSError,ValueError,KeyError,IndexError,TypeError,UnicodeError) as e:z={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:220]}
    z["scientific_status"]="HOLD_INDEPENDENT_VALIDATION_AND_NOVELTY";z["mode"]=a.mode
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n");print(json.dumps(z,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
