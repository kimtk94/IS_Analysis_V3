#!/usr/bin/env python3
"""Use exact public source schemas and author-documented barcode transformation; no individual outputs."""
import argparse,collections,csv,gzip,hashlib,io,json,math,re,statistics,time
from pathlib import Path
import early_track_orthogonal_sources as old
base=old.base;require=base.require;ROOT=old.ROOT
base.HOSTS.add("media.springernature.com")
SUP="https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12943-026-02604-x/MediaObjects/12943_2026_2604_MOESM3_ESM.xlsx"
SAMPLE={"sampleid","sample","samplename","sampleidentifier","samplecode","biopsyid","biopsy","biopsynumber","biopsycode","specimenid","specimen","tissueid"}
PATIENT={"patientid","patient","patientnumber","patientcode","subjectid","subject"}
TISSUE={"tissue","tissuetype","tissueclass","sampletype","biopsytype","histology","tissueorigin"}
GENES=old.GENE_HEADERS|{"genenames","genesymbol","symbol","hgnc","hgncgenesymbol","genenamesprimary","geneprimary","gene_name","genename"}
ANNO=old.ANNO|GENES|{"geneaccession","ensemblid","proteinaccession","proteinaccessions","proteinname","proteinid","proteinids","proteinname","protein group","protein groupid","fastaheaders","entrezgeneid","entrezid","description","transcriptid","organism","peptides","uniquepeptides","uniprotaccession","qvalue"}
def key(v):
    # Preserve exact publisher text and case. Numeric identifiers must be truly numeric integers.
    if isinstance(v,str) and v.strip()==v and v:return ("str",v)
    if isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v) and float(v).is_integer():return ("int",int(v))
    return None
def numeric(v):return isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v)
def canonical_tissue(v):
    k=old.norm(v)
    if k in {"tumor","tumour","hcc","tumorliver","tumourliver","tum","tumortissue","tumourbiopsy"}:return "tumor"
    if k in {"nontumorliver","nontumourliver","nontumor","nontumour","ntl","adjacentnontumor","adjacentnontumour","liverparenchyma","nontumortissue"}:return "non_tumor"
    if k in {"normal","normalliver","healthy","healthyliver"}:return "normal"
    return None
def read_supp():
    import openpyxl
    raw,src=base.get(SUP,15_000_000);wb=openpyxl.load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
    maps=[];profiles=[]
    try:
        for wi,ws in enumerate(wb,1):
            for ri,row in enumerate(ws.iter_rows(min_row=1,max_row=min(ws.max_row or 1,14),values_only=True),1):
                h=[old.norm(x) for x in row]
                a=[i for i,x in enumerate(h) if x in SAMPLE];b=[i for i,x in enumerate(h) if x in PATIENT];c=[i for i,x in enumerate(h) if x in TISSUE]
                if not a:continue
                profiles.append({"sheet":wi,"row":ri,"dimensions":[ws.max_row,ws.max_column],"sample_cols":len(a),"patient_cols":len(b),"tissue_cols":len(c)})
                for ai in a:
                    q=collections.defaultdict(list)
                    for v in ws.iter_rows(min_row=ri+1,values_only=True):
                        if ai>=len(v) or key(v[ai]) is None:continue
                        sid=key(v[ai]);pat=key(v[b[0]]) if len(b)==1 and b[0]<len(v) else None
                        tis=canonical_tissue(v[c[0]]) if len(c)==1 and c[0]<len(v) else None
                        q[sid].append((pat,tis))
                    if len(q)>=20 and all(len(v)==1 for v in q.values()):
                        maps.append((wi,ri,ai,{k:v[0] for k,v in q.items()}))
        return src,profiles,maps
    finally:wb.close()
def headers(wb,known):
    result=[];profile=[]
    for wi,ws in enumerate(wb,1):
        if (ws.max_row or 0)<100:continue
        rows=list(ws.iter_rows(min_row=1,max_row=min(ws.max_row,44),values_only=True))
        for ri,row in enumerate(rows[:12],1):
            following=rows[ri:ri+24]
            if len(following)<8:continue
            h=[old.norm(x) for x in row]
            genes=[i for i,v in enumerate(h) if v in GENES]
            # Candidate assay columns must have exact header keys, no known annotation terms,
            # and predominantly numerical expression in subsequent rows.
            opts=collections.defaultdict(list)
            for i,v in enumerate(row):
                k=key(v)
                if k is None or h[i] in ANNO:continue
                if k[0]=="int" and k not in known:continue
                if not isinstance(v,str) and k not in known:continue
                if sum(i<len(rr) and numeric(rr[i]) for rr in following)<.75*len(following):continue
                opts[k].append(i)
            ndup=sum(len(v)>1 for v in opts.values())
            columns={k:v[0] for k,v in opts.items() if len(v)==1}
            if len(columns)<20:continue
            profile.append({"sheet":wi,"header_row":ri,"dimensions":[ws.max_row,ws.max_column],
                "numeric_backed_keys":len(columns),"duplicate_keys":ndup,"exact_supplement_keys":len(set(columns)&known),
                "explicit_gene_columns":len(genes),"recognized_annotation_terms":sorted({v for v in h if v in ANNO})})
            if ndup==0:result.append((wi,ri,genes,columns))
    return result,profile
def vector_rows(wb,cand,keys):
    wi,ri,genes,cols=cand
    require(len(genes)==1,"unique explicitly labeled symbol column required")
    ws=wb.worksheets[wi-1];gi=genes[0]
    data=collections.defaultdict(list)
    for row in ws.iter_rows(min_row=ri+1,values_only=True):
        if gi>=len(row) or not isinstance(row[gi],str):continue
        g=row[gi].strip()
        if g not in old.PANEL:continue
        data[g].append([float(row[cols[k]]) if cols[k]<len(row) and numeric(row[cols[k]]) else None for k in keys])
    return data
def hcc_schema():
    import openpyxl
    from scipy.stats import spearmanr
    srcs={};books={}
    first=json.loads((ROOT/"docs/publication/evidence/early_track_orthogonal_first_readback_20261125.json").read_text())
    q=next(z["aggregate"] for z in first["records"] if z["aggregate"]["mode"]=="hcc")
    for name in old.ZEN:
        raw,src=old.zen_download(name)
        require(src["sha256"]==q["sources"][name]["sha256"],"previous workbook SHA drift")
        srcs[name]=src;books[name]=openpyxl.load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
    try:
        sup={}
        try:ss,sp,sm=read_supp();sup={"source":ss,"header_profiles":sp,"eligible_sample_maps":len(sm)}
        except (OSError,ValueError,KeyError,TypeError) as e:sm=[];sup={"status":"HOLD_SUPPLEMENT","reason":str(e)[:180]}
        known=set().union(*(set(m[3]) for m in sm)) if sm else set()
        rn,pn=list(old.ZEN);rc,rp=headers(books[rn],known);pc,pp=headers(books[pn],known)
        combos=[]
        for a in rc:
            for b in pc:
                common=set(a[3])&set(b[3])
                if len(common)>=1:combos.append((a,b,common))
        out={"status":"HCC_MULTISHEET_EXACT_KEY_AUDIT","sources":srcs,"sample_manifest":sup,"RNA_profiles":rp,"protein_profiles":pp,
             "candidate_cross_sheet_joins":[{"rna_sheet":a[0],"rna_header":a[1],"protein_sheet":b[0],"protein_header":b[1],"exact_keys":len(c),"manifest_keys":len(c&known)} for a,b,c in combos],
             "confirmed_biopsy_labels":0,"patient_tissue_context":[]}
        # Require exact manifest identifiers with nonmissing patient and tissue; multiple maps
        # may corroborate, but conflicting values are not resolved by guessing.
        meta=collections.defaultdict(set)
        for _,_,_,m in sm:
            for k,v in m.items():
                if v[0] is not None and v[1] is not None:meta[k].add(v)
        meta={k:next(iter(v)) for k,v in meta.items() if len(v)==1}
        eligible=[(a,b,c&set(meta)) for a,b,c in combos if len(c&set(meta))>=20 and len(a[2])==len(b[2])==1]
        # No double counting the same RNA or protein sheet, even if alternate header candidates match.
        eligible=[(a,b,c) for a,b,c in eligible if sum(x[0]==a[0] or y[0]==b[0] for x,y,z in eligible)==1]
        used=set()
        for a,b,c in eligible:
            if c&used:continue
            used.update(c);keys=sorted(c)
            rv=vector_rows(books[rn],a,keys);pv=vector_rows(books[pn],b,keys)
            report={"rna_sheet":a[0],"protein_sheet":b[0],"exact_manifest_labels":len(keys),"genes":{}}
            for g in old.PANEL:
                z={"RNA_rows":len(rv[g]),"protein_rows":len(pv[g]),"status":"HOLD_MISSING_OR_NONUNIQUE_GENE"}
                if len(rv[g])==len(pv[g])==1:
                    groups=collections.defaultdict(list)
                    for k,x,y in zip(keys,rv[g][0],pv[g][0]):
                        if x is not None and y is not None:groups[meta[k]].append((x,y))
                    ctx={}
                    for t in ("tumor","non_tumor"):
                        v=[(statistics.mean(x for x,y in rr),statistics.mean(y for x,y in rr)) for (p,tis),rr in groups.items() if tis==t]
                        if len(v)>=5 and len({x for x,y in v})>1 and len({y for x,y in v})>1:
                            rho=float(spearmanr([x for x,y in v],[y for x,y in v]).statistic)
                            if math.isfinite(rho):ctx[t]={"patients":len(v),"spearman_rho":rho}
                    z.update(status="DESCRIPTIVE_PATIENT_AGGREGATE_CONTEXT",patient_tissue_correlations=ctx)
                report["genes"][g]=z
            out["patient_tissue_context"].append(report)
        out["confirmed_biopsy_labels"]=len(used)
        out["limits"]=["Exact sample manifest/headers only, no forced aliases or across-sheet row-order matching. Patient means within tissue are descriptive, no p-values or causal claims. 2026 cohort explicitly reuses Ng 2022 patients and assays; potential overlap with GSE64041 remains unresolved. No repeated biopsies counted as independent patients."]
        return out
    finally:
        for w in books.values():w.close()
# Author merge order from R/OSCC-01-mergecells.Rmd, source commit bc7333009fce5c5746eb7ea75d546826683613ed.
LIBS=dict(zip(old.GSMS,range(1,10)))
SITE_CODES={
  1:{"1":"LN_MET"},2:{"2.1":"PT","2.2":"LN_MET","2.3":"LN_MET"},3:{"3":"PT"},
  4:{"4.1":"LN_MET","4.2":"LN_MET","4.3":"LN_FREE"},5:{"5":"PT"},
  6:{"6.1":"PT","6.2":"LN_MET","6.3":"LN_FREE"},7:{"6.4":"PT","6.5":"LN_FREE","6.6":"LN_FREE"},
  8:{"7.1":"PT","7.2":"LN_MET","7.3":"LN_MET"},9:{"7.4":"PT","7.5":"LN_FREE","7.6":"LN_MET"}}
def labels(raw,bars,lib):
    text=old.text_gz(raw);bi={v:i for i,v in enumerate(bars)}
    require(len(bi)==len(bars),"duplicate raw barcodes")
    maps={"identity":bi,"author_merge_index_suffix":{v+"_"+str(lib):i for v,i in bi.items()}}
    candidates=[];profiles=[]
    for delim in (",","\t",";"):
        rr=[r for r in csv.reader(io.StringIO(text),delimiter=delim) if r and any(v.strip() for v in r)]
        if not rr:continue
        widths=collections.Counter(map(len,rr));profiles.append({"delimiter":repr(delim),"row_width_counts":dict(widths),"rows":len(rr)})
        if len(widths)!=1 or len(rr[0])<2 or len(rr[0])>20:continue
        for k,mp in maps.items():
            for col in range(len(rr[0])):
                a=rr;first=rr[0][col].strip()
                hdr=None
                if first not in mp and old.norm(first) in {"","x","cell","cellid","barcode","barcodes","cellbarcode","rownames"}:
                    hdr=rr[0];a=rr[1:]
                if not a or not all(r[col] in mp for r in a) or len({r[col] for r in a})!=len(a):continue
                # A cell type column must contain known source-author Macro designation.
                lc=[i for i in range(len(a[0])) if i!=col and any(r[i]=="Macro" for r in a) and all(0<len(r[i])<=80 for r in a)]
                if len(lc)!=1:continue
                ci=lc[0];celltypes={mp[r[col]]:r[ci] for r in a}
                site={}
                if len(SITE_CODES[lib])==1:site={i:next(iter(SITE_CODES[lib].values())) for i in celltypes}
                elif hdr:
                    sc=[i for i,v in enumerate(hdr) if old.norm(v) in {"sample","sampleid","hashid","hashidman","hto","htoid"}]
                    if len(sc)==1 and all(r[sc[0]] in SITE_CODES[lib] for r in a):
                        site={mp[r[col]]:SITE_CODES[lib][r[sc[0]]] for r in a}
                candidates.append((repr(delim),k,col,ci,celltypes,site))
    if len(candidates)!=1:return None,{"delimiter_profiles":profiles,"complete_candidate_joins":len(candidates)}
    d,k,c,i,ct,site=candidates[0]
    return (ct,site),{"delimiter_profiles":profiles,"complete_candidate_joins":1,"mapping":k,"delimiter":d,"cell_column":c+1,"label_column":i+1,
                        "labeled_cells":len(ct),"Macro_cells":sum(x=="Macro" for x in ct.values()),"explicit_site_counts":dict(collections.Counter(site.values()))}
def get_files(sid,full=False):
    folder,names,listing=old.filelist(sid)
    parts=[("barcode",r"barcodes\.tsv\.gz$",15_000_000),("label",r"celltype_metadata.*\.csv\.gz$",5_000_000)]
    if full:parts.extend([("features",r"features\.tsv\.gz$",5_000_000),("matrix",r"matrix\.mtx\.gz$",60_000_000)])
    data={};src={}
    for k,p,cap in parts:data[k],src[k]=base.get(folder+old.one(names,p),cap)
    bars=[v[0] for v in csv.reader(io.StringIO(old.text_gz(data["barcode"])),delimiter="\t") if v]
    result,profile=labels(data["label"],bars,LIBS[sid])
    return data,src,bars,result,profile
def oral_labels():
    out=[]
    for sid in old.GSMS:
        z={"gsm":sid}
        try:
            _,src,bars,m,p=get_files(sid);z.update(status="AUTHOR_LABEL_JOIN" if m else "HOLD_LABEL_SCHEMA",sources=src,raw_barcodes=len(bars),**p)
        except (OSError,EOFError,ValueError,KeyError,IndexError,UnicodeError) as e:z.update(status="HOLD_SOURCE",reason=str(e)[:150])
        out.append(z);time.sleep(.4)
    return {"status":"ORAL_AUTHOR_LABEL_SCHEMA_AUDIT","libraries":9,"complete_joins":sum(v["status"]=="AUTHOR_LABEL_JOIN" for v in out),
            "explicit_site_libraries":sum(bool(v.get("explicit_site_counts")) for v in out),"records":out,
            "limits":["Author merge order, identity or exact library index suffix only; complete bijection required. Site assignment only single-site source or explicitly labeled author HTO sample code. No guessing from barcode/HTO maxima. Patients 4/6/7 have both LN_MET and LN_FREE per author code; maximum paired sign-test sample size three before eligibility, minimum two-sided p=.25."]}
def oral_pilot():
    import numpy as np
    from scipy.io import mmread
    sid=old.GSMS[0];data,src,bars,m,p=get_files(sid,True)
    out={"status":"ORAL_AUTHOR_LABEL_REPAIRED_PILOT","gsm":sid,"sources":src,"raw_barcodes":len(bars),**p}
    if m is None:return out
    f=[v for v in csv.reader(io.StringIO(old.text_gz(data["features"])),delimiter="\t") if v]
    require(all(len(v)>=2 for v in f),"feature schema")
    mat=mmread(io.BytesIO(gzip.decompress(data["matrix"])))
    require(mat.shape==(len(f),len(bars)) and np.isfinite(mat.data).all() and (mat.data>=0).all() and (mat.data==np.floor(mat.data)).all(),"matrix dimensions/count contract")
    nnz=mat.nnz;mat=mat.tocsc();require(mat.nnz==nnz,"duplicate coordinates")
    ct,site=m;macro=[i for i,g in ct.items() if g=="Macro"];genes=[i for i,v in enumerate(f) if len(v)<3 or v[2]=="Gene Expression"]
    ix={g:[i for i in genes if f[i][1]==g] for g in ("SPP1","CXCR4")}
    out.update(matrix_nonzero=int(nnz),features=len(f),mapped_label_cells=len(ct))
    if len(macro)>=50 and all(len(v)==1 for v in ix.values()):
        den=float(mat[genes,:][:,macro].sum());require(den>0,"macro RNA count denominator")
        out["macro_context"]={g:{"CPM":float(mat[v[0],:][:,macro].sum())*1e6/den,"cells_detected":int(mat[v[0],:][:,macro].count_nonzero())} for g,v in ix.items()}
    out["limits"]=["One metastatic lymph-node library, no paired effect. Full sparse matrix and exact labels; no patient/cell level output. No inference about all other libraries."]
    return out
def selftest():
    old.selftest()
    b=["AAAAAAAAAAAAAAAA-1","CCCCCCCCCCCCCCCC-1"]
    raw=gzip.compress(b'barcode,celltype,extra\nAAAAAAAAAAAAAAAA-1_1,Macro,1\nCCCCCCCCCCCCCCCC-1_1,T,2\n')
    m,p=labels(raw,b,1);assert m and p["mapping"]=="author_merge_index_suffix" and p["Macro_cells"]==1
    bad=gzip.compress(b'barcode,celltype\nAAAAAAAAAAAAAAAA-1_2,Macro\nCCCCCCCCCCCCCCCC-1_2,T\n')
    assert labels(bad,b,1)[0] is None
    assert sum(any(x=="LN_FREE" for x in c.values()) for c in SITE_CODES.values())==3
    assert canonical_tissue("Non-tumor liver")=="non_tumor"
    a=json.loads((ROOT/"docs/publication/evidence/early_track_orthogonal_first_readback_20261125.json").read_text())
    assert a["run_id"]==45611128274 and len(a["records"])==3
    print("Focused exact author suffix, reject wrong library, tissue schema and archive checks passed.")
def main():
    ap=argparse.ArgumentParser();ap.add_argument("mode",choices=["hcc_schema","oral_labels","oral_pilot","selftest"]);ap.add_argument("--out",type=Path);a=ap.parse_args()
    if a.mode=="selftest":selftest();return
    try:z=globals()[a.mode]()
    except (OSError,EOFError,ValueError,KeyError,IndexError,TypeError,UnicodeError) as e:z={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:220]}
    z.update(mode=a.mode,scientific_status="HOLD_PATIENT_LEVEL_INDEPENDENCE_AND_NOVELTY")
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n");print(json.dumps(z,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
