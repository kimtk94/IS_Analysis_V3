#!/usr/bin/env python3
"""Independent public data qualification. No patient or cell identifiers in output."""
import argparse,collections,csv,gzip,hashlib,io,json,math,re,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"prior/scripts/publication"))
import early_track_acceleration as base
require=base.require
base.HOSTS.update({"zenodo.org","www.ncbi.nlm.nih.gov","ftp.ncbi.nlm.nih.gov"})
PANEL=("SOAT1","HMGCR","SQLE","LDLR","FASN","ACACA","ACLY","SREBF2")
ZEN={"Transcriptomics_Tums_NTL_Normals.xlsx":"9012a374588a273bcb930b83643bf1b0","Proteomics_Tums_NTL_Normals.xlsx":"aee43ad590214e713b395c4f60d827a1"}
GSMS=("GSM5844648","GSM6553453","GSM6553455","GSM6553456","GSM6553457","GSM6553458","GSM6553459","GSM6553461","GSM6553462")
def norm(v):return re.sub(r"[^a-z0-9]","",str(v).lower())
GENE_HEADERS={"gene","genes","genesymbol","genesymbols","symbol","genename","genenames","hgncsymbol"}
ANNO=GENE_HEADERS|{"geneid","ensemblgeneid","ensembl","uniprot","uniprotid","protein","proteinid","proteinids","proteinnames","description","biotype","chromosome","chr","length","accession","majorityproteinids","id","index"}
def zen_download(name):
    require(name in ZEN,"exact public filename")
    url="https://zenodo.org/records/16992982/files/"+name+"?download=1"
    raw,src=base.get(url,90_000_000)
    require(hashlib.md5(raw).hexdigest()==ZEN[name],"published MD5 mismatch")
    return raw,{**src,"published_md5":ZEN[name]}
def workbook_headers(wb):
    cand=[];profiles=[]
    for wi,ws in enumerate(wb,1):
        hs=[]
        for ri,row in enumerate(ws.iter_rows(min_row=1,max_row=min(ws.max_row or 1,12),values_only=True),1):
            v=[str(x).strip() if x is not None else "" for x in row]
            ix=[i for i,x in enumerate(v) if norm(x) in GENE_HEADERS]
            labels={x:i for i,x in enumerate(v) if x and norm(x) not in ANNO}
            duplicates=len([x for x in v if x and norm(x) not in ANNO])-len(labels)
            if ix:hs.append({"row":ri,"explicit_gene_columns":len(ix),"candidate_keys":len(labels),"duplicate_candidate_keys":duplicates})
            if len(ix)==1 and len(labels)>=20 and duplicates==0:cand.append((wi,ri,ix[0],labels))
        profiles.append({"sheet":wi,"dimensions":[ws.max_row,ws.max_column],"header_candidates":hs})
    return cand,profiles
def gene_vectors(wb,cand,keys):
    wi,ri,gi,labels=cand;ws=wb.worksheets[wi-1];out=collections.defaultdict(list);counts=collections.Counter();n=0
    for row in ws.iter_rows(min_row=ri+1,values_only=True):
        n+=1
        if gi>=len(row) or not isinstance(row[gi],str):continue
        g=row[gi].strip()
        if g not in PANEL:continue
        vals=[]
        for k in keys:
            i=labels[k];x=row[i] if i<len(row) else None
            vals.append(float(x) if isinstance(x,(float,int)) and not isinstance(x,bool) and math.isfinite(x) else None)
        out[g].append(vals);counts[g]+=1
    return out,dict(counts),n
def hcc():
    import openpyxl
    from scipy.stats import spearmanr
    sources={};books={}
    for name in ZEN:
        raw,src=zen_download(name);sources[name]=src
        books[name]=openpyxl.load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
    rn,pn=list(ZEN)
    try:
        rc,rp=workbook_headers(books[rn]);pc,pp=workbook_headers(books[pn]);pairs=[]
        for a in rc:
            for b in pc:
                overlap=set(a[3])&set(b[3])
                if overlap:pairs.append((len(overlap),a,b,overlap))
        out={"status":"HCC_ZENODO_MATCHED_SOURCE_GATE","sources":sources,"RNA_workbook":rp,"protein_workbook":pp,
             "documented_RNA_specimens":231,"documented_protein_specimens":87,
             "candidate_overlap_sizes":sorted({x[0] for x in pairs},reverse=True),"matched_specimen_labels":0}
        # Strict expected complete public specimen coverage and unambiguous data sheets.
        exact=[x for x in pairs if x[0]==87 and len(x[1][3])==231 and len(x[2][3])==87]
        if len(exact)!=1:
            out["hold_reason"]="No unique complete 231 RNA / 87 protein exact-column join. No row-order or guessed identifier matching."
            return out
        n,a,b,keys=exact[0];keys=sorted(keys)
        rv,rnrows,nr=gene_vectors(books[rn],a,keys);pv,pnrows,np=gene_vectors(books[pn],b,keys)
        result={}
        for g in PANEL:
            z={"RNA_rows":rnrows.get(g,0),"protein_rows":pnrows.get(g,0),"status":"HOLD_MISSING_OR_NONUNIQUE_GENE"}
            if len(rv[g])==len(pv[g])==1:
                v=[(x,y) for x,y in zip(rv[g][0],pv[g][0]) if x is not None and y is not None]
                z.update(status="HOLD_INSUFFICIENT_FINITE_NONCONSTANT_VALUES",finite_shared=len(v))
                if len(v)>=20 and len({x for x,y in v})>1 and len({y for x,y in v})>1:
                    rho=float(spearmanr([x for x,y in v],[y for x,y in v]).statistic)
                    if math.isfinite(rho):z.update(status="DESCRIPTIVE_EXACT_SPECIMEN_LABEL_CORRELATION",spearman_rho=rho)
            result[g]=z
        out.update(matched_specimen_labels=n,RNA_rows_scanned=nr,protein_rows_scanned=np,results=result,
                   status="HCC_EXACT_LABEL_MATCH_DESCRIPTIVE")
        return out
    finally:
        for wb in books.values():wb.close()
def filelist(sid):
    folder="https://ftp.ncbi.nlm.nih.gov/geo/samples/"+sid[:-3]+"nnn/"+sid+"/suppl/"
    raw,src=base.get(folder,2_000_000);names=base.links(raw)
    files=[x for x in names if x.startswith(sid+"_") and re.fullmatch(r"[A-Za-z0-9_.-]+\.(?:csv|tsv|mtx)(?:\.gz)?",x)]
    return folder,files,src
def text_gz(raw):
    return gzip.decompress(raw).decode("utf-8-sig")
def one(files,pattern):
    z=[n for n in files if re.search(pattern,n,re.I)]
    require(len(z)==1,"nonunique "+pattern)
    return z[0]
def oral():
    inventory=[]
    for sid in GSMS:
        try:
            folder,files,src=filelist(sid)
            kinds={k:sum(bool(re.search(p,f,re.I)) for f in files) for k,p in {
                "barcodes":r"barcodes\.tsv\.gz$","features":r"features\.tsv\.gz$","matrix":r"matrix\.mtx\.gz$","author_labels":r"celltype_metadata.*\.csv\.gz$"}.items()}
            inventory.append({"gsm":sid,"status":"PUBLIC_FILE_INVENTORY","files":files,"kinds":kinds,"source":src})
        except (OSError,ValueError,KeyError) as e:inventory.append({"gsm":sid,"status":"HOLD_FILELIST","reason":str(e)[:150]})
        time.sleep(.4)
    return {"status":"ORAL_GERMAN_PUBLIC_SOURCE_INVENTORY","GSE":"GSE195655","libraries":len(GSMS),"inventory":inventory,
            "complete_library_file_sets":sum(all(v["kinds"].get(k)==1 for k in ("barcodes","features","matrix","author_labels")) for v in inventory if "kinds" in v),
            "limits":["Nine scRNA libraries, not nine individual patients; multi-region pooled HTO samples require author-verified demultiplexing. Reported 5 primary tumors, 9 metastatic and 5 tumor-free nodes. No raw FASTQ access or guessed specimen assignments."]}
def oral_pilot():
    import numpy as np
    from scipy.io import mmread
    sid=GSMS[0];folder,files,listing=filelist(sid)
    data={};sources={}
    for k,p,cap in (("barcode",r"barcodes\.tsv\.gz$",15_000_000),("features",r"features\.tsv\.gz$",5_000_000),("matrix",r"matrix\.mtx\.gz$",60_000_000),("label",r"celltype_metadata.*\.csv\.gz$",5_000_000)):
        name=one(files,p);data[k],sources[k]=base.get(folder+name,cap)
    bars=[v[0] for v in csv.reader(io.StringIO(text_gz(data["barcode"])),delimiter="\t") if v]
    feats=[v for v in csv.reader(io.StringIO(text_gz(data["features"])),delimiter="\t") if v]
    require(len(set(bars))==len(bars) and all(len(v)>=2 for v in feats),"barcode/feature identity")
    lab=list(csv.reader(io.StringIO(text_gz(data["label"]))))
    require(lab and all(len(v)==2 for v in lab),"two-column author metadata contract")
    if [norm(v) for v in lab[0]]==["cell","celltypeclassification"]:lab=lab[1:]
    require(len({v[0] for v in lab})==len(lab),"duplicate author cell labels")
    bi={v:i for i,v in enumerate(bars)};match={v[0]:v[1] for v in lab if v[0] in bi}
    out={"status":"ORAL_AUTHOR_LABEL_EXACT_JOIN_PILOT","gsm":sid,"listing":listing,"sources":sources,
         "barcode_columns":len(bars),"feature_rows":len(feats),"author_labeled_cells":len(lab),"exact_label_matches":len(match),
         "author_label_counts":dict(collections.Counter(v[1] for v in lab))}
    # Full integer sparse matrix shape and coordinate validation; no forced barcode normalization.
    mat=mmread(io.BytesIO(gzip.decompress(data["matrix"]))).tocsc()
    require(mat.shape==(len(feats),len(bars)) and np.isfinite(mat.data).all() and (mat.data>=0).all() and (mat.data==np.floor(mat.data)).all(),"matrix shape/count contract")
    out["matrix_nonzero"]=int(mat.nnz)
    if len(match)!=len(lab):
        out["hold_reason"]="Exact source barcode and author annotation key mismatch; no prefix or suffix guessing."
        return out
    macro=[bi[k] for k,v in match.items() if v=="Macro"]
    out["author_Macro_cells"]=len(macro)
    gi={g:[i for i,v in enumerate(feats) if v[1]==g and (len(v)<3 or v[2]=="Gene Expression")] for g in ("CXCR4","SPP1")}
    ge=[i for i,v in enumerate(feats) if len(v)<3 or v[2]=="Gene Expression"]
    if len(macro)>=50 and all(len(v)==1 for v in gi.values()):
        denom=float(mat[ge,:][:,macro].sum())
        require(denom>0,"nonzero aggregate RNA counts")
        out["macro_context"]={g:{"CPM":float(mat[v[0],:][:,macro].sum())*1e6/denom,"cells_detected":int(mat[v[0],:][:,macro].count_nonzero())} for g,v in gi.items()}
    out["limits"]=["One preselected library, not an independent patient-level effect. Author Macro label exact only. No inference of primary/metastatic/tumor-free site from library name or barcode."]
    return out
def selftest():
    import openpyxl
    w=openpyxl.Workbook();s=w.active;s.append(["Gene Symbol"]+["sample_"+str(i) for i in range(25)])
    s.append(["SOAT1"]+list(range(25)))
    c,p=workbook_headers(w)
    assert len(c)==1 and len(c[0][3])==25
    v,n,k=gene_vectors(w,c[0],["sample_1","sample_2"])
    assert n=={"SOAT1":1} and v["SOAT1"]==[[1.0,2.0]]
    assert len(PANEL)==8 and len(GSMS)==9 and len(ZEN)==2
    a=json.loads((ROOT/"docs/publication/evidence/early_track_genetic_first_readback_20261125.json").read_text())
    assert a["run_id"]==45610423960 and len(a["records"])==7
    print("Exact sample key, fixed gene, public source and archive checks passed.")
def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["hcc","oral","oral_pilot","selftest"]);p.add_argument("--out",type=Path);a=p.parse_args()
    if a.mode=="selftest":selftest();return
    try:z=globals()[a.mode]()
    except (OSError,EOFError,ValueError,KeyError,TypeError,IndexError,UnicodeError) as e:z={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:220]}
    z.update(mode=a.mode,scientific_status="HOLD_PATIENT_LEVEL_INDEPENDENCE_AND_NOVELTY")
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n");print(json.dumps(z,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
