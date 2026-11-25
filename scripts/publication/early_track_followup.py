#!/usr/bin/env python3
"""Follow verified first-wave inputs: full public oral matrices and HCC supplement schema.
No cell/participant identifiers or row-level clinical values are exported.
"""
import argparse,collections,gzip,hashlib,io,json,math,re,sys,time
from pathlib import Path
import early_track_acceleration as base

ROOT=Path(__file__).resolve().parents[2]
ARCH=ROOT/"docs/publication/evidence/early_track_first_readback_20261125.json"
PANEL={"known_metastasis_context":["SPP1","CXCR4"],"epithelial":["EPCAM","KRT5","KRT14","KRT19"],
       "myeloid":["LST1","TYROBP","FCER1G","CD68"],"T_cell":["CD3D","CD3E","TRAC"],
       "B_cell":["MS4A1","CD79A","CD74"],"fibroblast":["COL1A1","COL1A2","DCN"]}
def prior(track):
    data=json.loads(ARCH.read_text())
    base.require(data["provenance"]["run_id"]==45597748445,"first run identity")
    hits=[j["aggregate"] for j in data["completed"] if j["aggregate"].get("track")==track]
    base.require(len(hits)==1,"archived track missing")
    return hits[0]

def oral(group):
    import numpy as np
    import scipy.io
    from scipy import sparse
    old=prior("oral")
    rows,ids,titles,src=base.geo_metadata("GSE322620")
    base.require(src["matrix"]["sha256"]==old["source"]["matrix"]["sha256"],"GEO metadata drift")
    base.require(titles==["PT%d"%i for i in range(1,6)]+["MET%d"%i for i in range(1,6)],"exact title order drift")
    files=base.sample_files(rows,ids)
    records=[]
    for i,(sid,title) in enumerate(zip(ids,titles)):
        if not title.startswith(group):continue
        p=old["records"][i]
        result={"specimen":title,"group":group,"status":"HOLD_FULL_MATRIX"}
        try:
            base.require(p["status"]=="DIMENSION_PREFIX_VERIFIED","prior dimension gate")
            f,fm=base.get(files[sid]["features"],3_000_000)
            b,bm=base.get(files[sid]["barcodes"],3_000_000)
            base.require(fm["sha256"]==p["feature_sha256"] and bm["sha256"]==p["barcode_sha256"],"feature/barcode source drift")
            features=[v.split("\t") for v in gzip.decompress(f).decode().splitlines() if v.strip()]
            bars=gzip.decompress(b).decode().splitlines()
            base.require(len(features)==p["genes"] and len(bars)==p["cells"],"dimension drift")
            base.require(all(len(x)==3 and x[2]=="Gene Expression" for x in features),"not a gene expression only feature universe")
            syms=[x[1] for x in features]; genes=collections.Counter(syms)
            raw,source=base.get(files[sid]["matrix"],150_000_000)
            base.require(hashlib.sha256(raw[:65536]).hexdigest()==p["matrix_prefix_sha256"],"matrix prefix drift")
            x=scipy.io.mmread(gzip.GzipFile(fileobj=io.BytesIO(raw)))
            base.require(sparse.issparse(x) and x.shape==(p["genes"],p["cells"]) and x.nnz==p["declared_nonzero"],"full matrix dimensions/nnz mismatch")
            base.require(np.issubdtype(x.data.dtype,np.integer) and np.all(x.data>0),"nonpositive/noninteger count")
            m=x.tocsc();base.require(m.nnz==x.nnz,"duplicate coordinates")
            del x,raw
            totals=np.asarray(m.sum(axis=0)).ravel()
            detected=np.diff(m.indptr)
            mt=np.array([s.upper().startswith("MT-") for s in syms])
            base.require(int(mt.sum())>0,"mitochondrial gene universe missing")
            mito=np.asarray(m[mt,:].sum(axis=0)).ravel()
            good=(detected>=200)&(totals>=500)&(mito<=.2*totals)
            base.require(int(good.sum())>0,"no QC eligible cells")
            counts=np.asarray(m[:,good].sum(axis=1)).ravel()
            denom=int(counts.sum());base.require(denom>0,"empty pseudobulk")
            values={}
            for family,markers in PANEL.items():
                vals={}
                for g in markers:
                    if genes[g]==1:
                        v=int(counts[syms.index(g)])
                        vals[g]={"counts":v,"cpm":v*1e6/denom,"log2_cpm_plus1":math.log2(v*1e6/denom+1)}
                    else:vals[g]={"status":"HOLD_NONUNIQUE_OR_ABSENT"}
                values[family]=vals
            result.update(status="FULL_MATRIX_QC_AND_SPECIMEN_PSEUDOBULK",source=source,genes=p["genes"],cells=p["cells"],
                          declared_nonzero=p["declared_nonzero"],qc_cells=int(good.sum()),qc_umis=denom,
                          mitochondrial_features=int(mt.sum()),panel=values)
        except (OSError,ValueError,TypeError,IndexError,KeyError,UnicodeError) as e:
            result["reason"]=str(e)[:150]
        records.append(result);time.sleep(1)
    return {"status":"ORAL_DESCRIPTIVE_WHOLE_SPECIMEN_CONTEXT","group":group,"records":records,
            "qc":{"min_detected_genes":200,"min_UMI":500,"max_mito_fraction":.2},
            "panel":PANEL,"limits":["Specimens, not known independent/pair-matched donors. No cell-level output, no tests, no malignant-cell attribution, no cell type deconvolution or clinical prediction. Whole-specimen CPM is confounded by tissue composition, tumor purity and lymph node background. No doublet/ambient RNA correction; simple QC only. Full file SHA recorded, prefix parity is not publisher full checksum. Known SPP1/CXCR4 axis (DOI 10.1002/advs.202400524) is context, not novel."]}

def hcc():
    import openpyxl
    old=prior("hcc"); xml=old["xml"]
    raw,src=base.get(xml["url"],xml["bytes"])
    base.require(src["sha256"]==xml["sha256"],"publisher XML drift")
    prefix="https://pmc-oa-opendata.s3.amazonaws.com/PMC9068765.1/"
    vocab=("gene","protein","uniprot","ensembl","symbol","rna","mrna","transcript","proteom","phospho","logfc","log2fc","fold change",
           "p.value","p-value","p value","adjusted","fdr","sample","patient","tumor","tumour","survival","bclc","stage","event","follow",
           "clinical","etiology","grade","biopsy","id","correlation")
    results=[]
    for item in old["inventory"]:
        n=item["file"];size=item["bytes"]
        row={"file":n,"caption":item["caption"],"status":"HOLD_WORKBOOK"}
        try:
            base.require(re.fullmatch(r"41467_2022_29960_MOESM\d+_ESM\.xlsx",n) and 0<size<=12_000_000,"unexpected publisher inventory")
            data,meta=base.get(prefix+n,size)
            etag=(meta.get("etag") or "").strip('"')
            base.require(len(data)==size and re.fullmatch("[0-9a-f]{32}",etag) and hashlib.md5(data).hexdigest()==etag,"object size/MD5 mismatch")
            wb=openpyxl.load_workbook(io.BytesIO(data),read_only=True,data_only=True)
            sheets=[]
            for j,ws in enumerate(wb,1):
                rows=[]
                for k,rr in enumerate(ws.iter_rows(min_row=1,max_row=min(ws.max_row or 1,8),values_only=True),1):
                    tokens=collections.Counter()
                    for v in rr:
                        if not isinstance(v,str) or len(v)>150:continue
                        v=v.lower()
                        for term in vocab:
                            if re.search(r"\b"+re.escape(term)+r"\b",v):tokens[term]+=1
                    rows.append({"row":k,"nonempty":sum(v is not None for v in rr),"recognized_terms":dict(tokens)})
                sheets.append({"sheet_index":j,"dimensions":[ws.max_row,ws.max_column],"first_rows":rows})
            wb.close();row.update(status="PUBLIC_WORKBOOK_MD5_AND_SCHEMA_READ",source=meta,sheets=sheets)
        except (OSError,ValueError,TypeError,KeyError,IndexError) as e:row["reason"]=str(e)[:150]
        results.append(row);time.sleep(1)
    return {"status":"HCC_FULL_PUBLIC_SUPPLEMENT_SCHEMA_AUDIT","workbooks":results,
            "verified":sum(x["status"]=="PUBLIC_WORKBOOK_MD5_AND_SCHEMA_READ" for x in results),
            "limits":["Header concepts and dimensions only; no clinical values, no personal identifiers, no raw EGA access, no expression/patient join. Next gate: exact Data 2 gene-level matched mRNA/protein effect schema and genuinely independent comparison, or approved RNA matrix access. Ng et al. already performed discordance/clinical clustering."]}

def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["PT","MET","hcc","selftest"]);p.add_argument("--out",type=Path);a=p.parse_args()
    if a.mode=="selftest":
        assert len(set(g for panel in PANEL.values() for g in panel))==sum(map(len,PANEL.values()))
        assert len(prior("oral")["records"])==10 and len(prior("hcc")["inventory"])==15
        print("Focused fixed-panel and archive checks passed");return
    base.require(a.out is not None,"output path required")
    try:x=hcc() if a.mode=="hcc" else oral(a.mode)
    except (OSError,ValueError,TypeError,KeyError,IndexError,UnicodeError) as e:x={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:200]}
    x["scientific_status"]="HOLD_NO_INDEPENDENT_INFERENCE"
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
    print(json.dumps(x,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
