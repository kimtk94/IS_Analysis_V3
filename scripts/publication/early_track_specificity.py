#!/usr/bin/env python3
"""Third-wave fixed compartment-specific descriptive oral analysis and HCC published-effects audit."""
import argparse,collections,gzip,hashlib,io,json,math,re,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"prior_helpers/scripts/publication"))
import early_track_acceleration as base
ARCH=ROOT/"docs/publication/evidence/early_track_second_readback_20261125.json"
FIRST=ROOT/"prior_helpers/docs/publication/evidence/early_track_first_readback_20261125.json"
GATES={
 "epithelial":["EPCAM","KRT5","KRT14","KRT8","KRT18","KRT19"],
 "myeloid":["LST1","TYROBP","FCER1G","CD68"],
 "T_cell":["CD3D","CD3E","TRAC"],
 "B_cell":["MS4A1","CD79A","CD79B","CD19"],
 "fibroblast":["COL1A1","COL1A2","DCN","LUM"],
 "endothelial":["PECAM1","VWF","KDR","EMCN"]}
# Literature-defined contextual marker score, not a novel signature; separate from lineage gates.
PEMT=["SERPINE1","TGFBI","MMP10","LAMC2","P4HA2","PDPN","ITGA5","LAMA3","CDH13","TNC","MMP2","EMP3","INHBA","LAMB3","VIM"]
CONTEXT=["SPP1","CXCR4"]
def archive():
    x=json.loads(ARCH.read_text())
    base.require(x["pinned_prior_head"]=="9703492179f0a29ecb7a332a24d474c49a7c10f7","prior head mismatch")
    return x["records"]
def prior_group(group):
    hits=[x["aggregate"] for x in archive() if x["aggregate"].get("group")==group and x["run_id"]==45598152083]
    base.require(len(hits)==1 and len(hits[0]["records"])==5,"fixed prior group")
    return hits[0]["records"]
def prior_hcc():
    hits=[x["aggregate"] for x in archive() if x["job_id"]==139149611159]
    base.require(len(hits)==1 and hits[0]["verified"]==15,"HCC prior archive")
    return hits[0]
def marker_index(syms):
    counter=collections.Counter(syms)
    out={}
    for g in set(sum(GATES.values(),[])+PEMT+CONTEXT+["PTPRC"]):
        base.require(counter[g]==1,"missing or nonunique fixed gene "+g)
        out[g]=syms.index(g)
    return out
def gate_cells(m,ix,good,min_markers):
    import numpy as np
    # Detection only, no outcome genes in assignment. Exclusive gates; mixed/unclassified withheld.
    lineage=np.vstack([np.asarray((m[[ix[g] for g in gs],:]>0).sum(axis=0)).ravel()>=min_markers for gs in GATES.values()])
    immune=np.asarray(m[ix["PTPRC"],:].toarray()).ravel()>0
    hits=lineage.sum(axis=0)
    out={}
    for j,key in enumerate(GATES):
        mask=good & (hits==1) & lineage[j]
        if key in ("epithelial","fibroblast","endothelial"):mask &= ~immune
        out[key]=mask
    total=sum(int(v.sum()) for v in out.values())
    base.require(total<=int(good.sum()),"exclusive gates violated")
    return out,{"qc_cells":int(good.sum()),"assigned":total,"mixed":int((good&(hits>1)).sum()),"unassigned":int(good.sum())-total}
def read_specimen(row,files,first):
    import numpy as np
    import scipy.io
    from scipy import sparse
    f,fm=base.get(files["features"],3_000_000);b,bm=base.get(files["barcodes"],3_000_000)
    base.require(fm["sha256"]==first["feature_sha256"] and bm["sha256"]==first["barcode_sha256"],"feature/barcode changed")
    feat=[v.split("\t") for v in gzip.decompress(f).decode().splitlines() if v.strip()]
    bars=gzip.decompress(b).decode().splitlines()
    base.require(len(feat)==row["genes"] and len(bars)==row["cells"] and len(set(bars))==len(bars),"dimensions changed")
    base.require(all(len(v)==3 and v[2]=="Gene Expression" for v in feat) and len({v[0] for v in feat})==len(feat),"feature contract")
    syms=[v[1] for v in feat];ix=marker_index(syms)
    raw,src=base.get(files["matrix"],150_000_000)
    base.require(src["sha256"]==row["source"]["sha256"] and src["bytes"]==row["source"]["bytes"],"full matrix SHA/size drift")
    x=scipy.io.mmread(gzip.GzipFile(fileobj=io.BytesIO(raw)))
    base.require(sparse.issparse(x) and x.shape==(row["genes"],row["cells"]) and x.nnz==row["declared_nonzero"] and np.issubdtype(x.data.dtype,np.integer) and np.all(x.data>0),"matrix invalid")
    m=x.tocsc();base.require(m.nnz==x.nnz,"duplicate coordinates");del x,raw
    totals=np.asarray(m.sum(axis=0)).ravel();detected=np.diff(m.indptr)
    mt=np.array([s.upper().startswith("MT-") for s in syms]);mito=np.asarray(m[mt,:].sum(axis=0)).ravel()
    good=(detected>=200)&(totals>=500)&(mito<=.2*totals)
    base.require(int(good.sum())==row["qc_cells"],"prior QC replay differs")
    results={}
    for minimum in (2,3):
        masks,info=gate_cells(m,ix,good,minimum)
        vals={}
        for key,mask in masks.items():
            n=int(mask.sum())
            d={"cells":n,"fraction_qc":n/int(good.sum()),"status":"HOLD_FEWER_THAN_50_CELLS"}
            if n>=50:
                # Aggregate raw counts; denominator across all genes, not panel.
                denom=int(totals[mask].sum())
                genes=PEMT+CONTEXT if key=="epithelial" else CONTEXT
                counts=np.asarray(m[[ix[g] for g in genes],:][:,mask].sum(axis=1)).ravel()
                cpms={g:int(v)*1e6/denom for g,v in zip(genes,counts)}
                d.update(status="DESCRIPTIVE_MARKER_GATE_PSEUDOBULK",umis=denom,cpm=cpms)
                if key=="epithelial":d["pemt_mean_log2_cpm_plus1"]=sum(math.log2(cpms[g]+1) for g in PEMT)/len(PEMT)
            vals[key]=d
        results[str(minimum)]={"summary":info,"compartments":vals}
    return {"status":"EXACT_MATRIX_REPLAY_AND_MARKER_GATES","specimen":row["specimen"],"group":row["group"],"matrix_sha256":src["sha256"],"thresholds":results}
def oral(group):
    old=prior_group(group);first=json.loads(FIRST.read_text())
    q=[v["aggregate"] for v in first["completed"] if v["aggregate"].get("track")=="oral"]
    base.require(len(q)==1 and len(q[0]["records"])==10,"first oral archive")
    meta,ids,titles,src=base.geo_metadata("GSE322620")
    base.require(src["matrix"]["sha256"]==q[0]["source"]["matrix"]["sha256"],"GEO source drift")
    base.require(titles==["PT%d"%i for i in range(1,6)]+["MET%d"%i for i in range(1,6)],"title mismatch")
    files=base.sample_files(meta,ids);results=[]
    for row in old:
        try:
            base.require(row["status"]=="FULL_MATRIX_QC_AND_SPECIMEN_PSEUDOBULK","prior full matrix missing")
            i=titles.index(row["specimen"])
            z=read_specimen(row,files[ids[i]],q[0]["records"][i])
        except (OSError,ValueError,KeyError,IndexError,TypeError,UnicodeError) as e:
            z={"status":"HOLD_SPECIMEN","specimen":row["specimen"],"reason":str(e)[:180]}
        results.append(z);time.sleep(1)
    return {"status":"ORAL_MARKER_DEFINED_COMPARTMENT_SENSITIVITY","group":group,"gates":GATES,"pemt_context":PEMT,"context":CONTEXT,
            "records":results,"limits":["GEO design wording ambiguous; no donor pairing or independence claim. Exclusive detected-lineage markers, not validated clustering, malignant CNV or cell-type reference label transfer. Mixed/unknown cells excluded; dropout and tissue composition remain. Fixed 2-marker primary and 3-marker sensitivity, 50 cells minimum; a 3-marker gate for T cells is stricter than others. No cell-level exports, patient-level inference, tests, ligand-receptor claims, or newly discovered signature. pEMT context already published (PMC13081553, Puram 2017)."]}
def norm(v):
    return re.sub(r"[^a-z0-9]+"," ",str(v).lower()).strip()
def exact_effect_header(rr):
    # Strict explicit RNA and protein signed log-fold-change columns; no positional guessing.
    out={}
    for i,v in enumerate(rr):
        s=norm(v)
        if not s:continue
        k=None
        if s in {"gene","gene name","gene symbol","symbol","hgnc symbol"}:k="gene"
        elif re.search(r"\b(?:mrna|rna|transcriptome|transcript)\b",s) and re.search(r"\b(?:log2 ?fc|log ?fc|log2 fold change|log fold change)\b",s):k="rna"
        elif re.search(r"\b(?:protein|proteome|proteomic)\b",s) and re.search(r"\b(?:log2 ?fc|log ?fc|log2 fold change|log fold change)\b",s):k="protein"
        if k:
            if k in out:return {}
            out[k]=i
    return out if set(out)=={"gene","rna","protein"} else {}
def finite(v):
    return isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v)
def effect_summary(rows,idx):
    from scipy.stats import spearmanr
    pool=collections.defaultdict(list);invalid=0
    for rr in rows:
        if len(rr)<=max(idx.values()):invalid+=1;continue
        g,a,b=(rr[idx[k]] for k in ("gene","rna","protein"))
        if not isinstance(g,str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,39}",g.strip()) or not finite(a) or not finite(b):invalid+=1;continue
        pool[g.strip()].append((float(a),float(b)))
    unique={k:v[0] for k,v in pool.items() if len(v)==1}
    base.require(len(unique)>=100,"fewer than 100 unique effect rows")
    a,b=zip(*unique.values());q=collections.Counter(("up" if u>0 else "down" if u<0 else "zero")+"_"+("up" if v>0 else "down" if v<0 else "zero") for u,v in zip(a,b))
    rho=float(spearmanr(a,b).statistic)
    base.require(math.isfinite(rho),"nonfinite correlation")
    return {"unique_genes":len(unique),"duplicate_gene_keys_withheld":sum(len(v)>1 for v in pool.values()),"invalid_rows":invalid,
            "quadrants":dict(q),"quadrant_fractions":{k:v/len(unique) for k,v in q.items()},"spearman_rho":rho,
            "note":"Descriptive across published gene effects, not biological replicates or a new association."}
def hcc():
    import openpyxl
    p=prior_hcc();selected=[v for v in p["workbooks"] if v["caption"] in {"Source Data","Supplementary Data 2","Supplementary Data 6"}]
    base.require(len(selected)==3,"pinned candidate workbooks")
    out=[]
    for w in selected:
        raw,src=base.get(w["source"]["url"],w["source"]["bytes"])
        base.require(src["sha256"]==w["source"]["sha256"],"workbook SHA drift")
        wb=openpyxl.load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
        for si,ws in enumerate(wb,1):
            rr=list(ws.iter_rows(min_row=1,max_row=min(8,ws.max_row),values_only=True))
            # Only allow whitelisted aggregate figure/gene-statistics header strings, no participant IDs.
            recognized=[]
            for ri,row in enumerate(rr,1):
                terms=[{"column":i+1,"label":str(v)} for i,v in enumerate(row) if isinstance(v,str) and len(v)<=80 and re.search(r"\b(?:gene symbol|gene name|log2fc|logfc|log2 fold change|mrna|transcriptome|proteome|protein|uniprot|spearman|fdr|adjusted p|correlation)\b",v,re.I) and not re.search(r"patient|sample|tumou?r\s*\d|case|donor",v,re.I)]
                if terms:recognized.append({"row":ri,"terms":terms})
            name=ws.title
            figure1f=bool(re.fullmatch(r"(?:fig(?:ure)?[ ._-]*)?1[ ._-]*f",name.strip(),re.I))
            z={"workbook":w["caption"],"sheet_index":si,"dimensions":[ws.max_row,ws.max_column],"figure_1f_exact_title":figure1f,"aggregate_header_candidates":recognized[:3]}
            if figure1f:
                candidates=[(i,exact_effect_header(row)) for i,row in enumerate(rr,1) if exact_effect_header(row)]
                if len(candidates)==1:
                    i,ix=candidates[0]
                    try:z["effect_replay"]=effect_summary(ws.iter_rows(min_row=i+1,values_only=True),ix)
                    except (ValueError,TypeError,IndexError) as e:z["effect_hold"]=str(e)[:150]
                else:z["effect_hold"]="no unique explicit signed RNA/protein effect header; no guessing"
            out.append(z)
        wb.close()
    return {"status":"HCC_PUBLISHED_EFFECT_SOURCE_REPLAY_GATE","sheets":out,"replayed":sum("effect_replay" in v for v in out),
            "limits":["Ng et al. DOI 10.1038/s41467-022-29960-8 Fig 1f already reports rho 0.33 and four directional quadrants. Only exact labeled figure source with explicit signed columns can be replayed; no substitute sheet, no patient values or clinical association. Header inspection is not paired RNA-protein individual data. Independent cohort and matched sample provenance remain required."]}
def external_hcc():
    acc="GSE124535"
    meta,ids,titles,src=base.geo_metadata(acc)
    folder=base.geo_folder(acc)+"suppl/"
    raw,listing=base.get(folder,2_000_000)
    names=[v for v in base.links(raw) if re.fullmatch(re.escape(acc)+r"_[A-Za-z0-9_.-]+",v)]
    return {"status":"HCC_JIANG_INDEPENDENT_CANDIDATE_SOURCE","geo":src,"samples":len(ids),"characteristic_fields":base.field_profile(meta,ids),
            "supplement_filenames":names,"supplement_directory":listing,"limits":["Jiang 2019 Nature DOI 10.1038/s41586-019-0987-8 lists GSE124535 and PRIDE PXD006512/PXD008373. RNA metadata does not establish protein pairing or independence from other Chinese cohorts; no full RNA/protein files or patient data exported."]}
def selftest():
    import numpy as np
    from scipy import sparse
    assert not set(sum(GATES.values(),[]))&set(PEMT+CONTEXT)
    syms=list(dict.fromkeys(sum(GATES.values(),[])+PEMT+CONTEXT+["PTPRC"]))
    ix=marker_index(syms);a=np.zeros((len(syms),4),dtype=int)
    for g in ["EPCAM","KRT5"]:a[ix[g],0]=1;a[ix[g],1]=1
    for g in ["LST1","TYROBP"]:a[ix[g],1]=1;a[ix[g],2]=1
    a[ix["PTPRC"],2]=1
    masks,info=gate_cells(sparse.csc_matrix(a),ix,np.ones(4,dtype=bool),2)
    assert masks["epithelial"].tolist()==[True,False,False,False] and masks["myeloid"].tolist()==[False,False,True,False] and info["assigned"]==2
    assert exact_effect_header(["Gene Symbol","mRNA log2FC","Protein log2FC"])=={"gene":0,"rna":1,"protein":2}
    assert exact_effect_header(["Gene","RNA","Protein"])=={}
    assert len(prior_group("PT"))==5 and len(prior_group("MET"))==5
    print("Fixed marker gates, exclusive assignment, signed-effect header and pinned archive assertions passed.")
def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["PT","MET","hcc","external_hcc","selftest"]);p.add_argument("--out",type=Path);a=p.parse_args()
    if a.mode=="selftest":selftest();return
    try:x=oral(a.mode) if a.mode in {"PT","MET"} else globals()[a.mode]()
    except (OSError,ValueError,KeyError,TypeError,IndexError,UnicodeError) as e:x={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:200]}
    x["scientific_status"]="HOLD_NO_INDEPENDENT_INFERENCE"
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps(x,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
