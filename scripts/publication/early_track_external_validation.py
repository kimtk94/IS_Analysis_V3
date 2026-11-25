#!/usr/bin/env python3
"""Fixed HCC external transcript validation and published protein source qualification."""
import argparse,collections,csv,gzip,hashlib,io,json,math,re,statistics,sys,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"prior_helpers/scripts/publication"))
import early_track_acceleration as base
require=base.require
PANEL=("SOAT1","HMGCR","SQLE","LDLR","FASN","ACACA","ACLY","SREBF2")
ARCH=ROOT/"docs/publication/evidence/early_track_matrix_readout_20261125.json"
base.HOSTS.update({"media.springernature.com","ftp.pride.ebi.ac.uk"})
def prior(mode):
    a=json.loads(ARCH.read_text());require(a["run_id"]==45607454164,"archive identity")
    x=[j["aggregate"] for j in a["records"] if j["aggregate"].get("mode")==mode]
    require(len(x)==1,"missing prior");return x[0]
def one_field(meta,key,n):
    x=meta.get(key,[]);require(len(x)==1 and len(x[0])==n,"missing field "+key);return x[0]
def characteristic(meta,name,n):
    pat=re.compile(re.escape(name)+r":\s*(.+)",re.I)
    x=[[pat.fullmatch(v.strip())[1].strip() for v in rr] for rr in meta.get("!Sample_characteristics_ch1",[]) if len(rr)==n and all(pat.fullmatch(v.strip()) for v in rr)]
    require(len(x)==1,"missing characteristic "+name);return x[0]
def decode_table(raw,begin,end):
    stream=io.TextIOWrapper(gzip.GzipFile(fileobj=io.BytesIO(raw)),encoding="utf-8-sig",errors="strict")
    for line in stream:
        if line.strip()==begin:break
    else:raise ValueError("table begin absent")
    rr=csv.reader(stream,delimiter="\t")
    header=next(rr);require(len(set(header))==len(header),"duplicate table header")
    def rows():
        for row in rr:
            if row==[end]:return
            require(len(row)==len(header),"table width drift")
            yield row
        raise ValueError("table end absent")
    return header,rows()
def sized_gzip(url,limit):
    # Exact remote size from a one-byte Range, then bounded complete transfer.
    # Full gzip reads in decode_table validate the end marker and CRC; no arbitrary cap escalation.
    import urllib.error
    error=None
    for attempt in range(2):
        try:
            raw,first=base.get(url,1,(0,0))
            m=re.fullmatch(r"bytes 0-0/(\d+)",first["content_range"] or "")
            require(m is not None,"missing exact remote size")
            n=int(m[1]);require(0<n<=limit,f"source too large: {n} bytes exceeds {limit}")
            raw,src=base.get(url,n)
            require(len(raw)==n and raw[:2]==b"\x1f\x8b","size/gzip identity mismatch")
            if first.get("etag") and src.get("etag"):
                require(first["etag"]==src["etag"],"remote entity changed during transfer")
            # Exhaust gzip in bounded chunks to verify CRC without expanding in memory.
            total=0
            with gzip.GzipFile(fileobj=io.BytesIO(raw)) as f:
                for chunk in iter(lambda:f.read(1_048_576),b""):
                    total+=len(chunk);require(total<=2_000_000_000,"gzip expansion cap")
            src.update(remote_size=n,uncompressed_bytes=total,range_probe=first)
            return raw,src
        except (OSError,EOFError) as e:
            error=e
            if attempt==0:time.sleep(3)
    raise error

def platform():
    url="https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/soft/GPL6244_family.soft.gz"
    raw,src=sized_gzip(url,180_000_000)
    h,rows=decode_table(raw,"!platform_table_begin","!platform_table_end")
    require("ID" in h and "gene_assignment" in h,"official GPL6244 annotation schema")
    i=h.index("ID");j=h.index("gene_assignment");m=collections.defaultdict(list);seen=set()
    for row in rows:
        key=row[i];require(key not in seen,"duplicate platform probe ID");seen.add(key)
        # Affymetrix format: accession // symbol // description // cytoband // Entrez ID.
        parts=[v.strip().split(" // ") for v in row[j].split(" /// ") if v.strip()]
        if not parts or any(len(v)<2 for v in parts):continue
        symbols={v[1].strip() for v in parts}
        if len(symbols)==1 and next(iter(symbols)) in PANEL:m[next(iter(symbols))].append(key)
    require(len(seen)>=30000,"incomplete GPL6244")
    return {g:sorted(set(m[g])) for g in PANEL},src,len(seen)
def sign_summary(delta,scale):
    from scipy.stats import binomtest
    n=len(delta);nonzero=[x for x in delta if x!=0];positive=sum(x>0 for x in nonzero)
    require(n>=5,"insufficient paired patients")
    med=statistics.median(delta);loo=[statistics.median(delta[:i]+delta[i+1:]) for i in range(n)]
    p=float(binomtest(positive,len(nonzero),.5).pvalue) if nonzero else 1.0
    return {"pairs":n,"positive":positive,"negative":len(nonzero)-positive,"ties":n-len(nonzero),
            "median_delta":med,"scale":scale,"loo_median_range":[min(loo),max(loo)],"p_exact_sign":p,"bonferroni_8":min(1,8*p)}
def swiss():
    import numpy as np
    acc="GSE64041";folder=base.geo_folder(acc)+"matrix/"
    raw,listing=base.get(folder,2_000_000)
    files=[v for v in base.links(raw) if re.fullmatch(r"GSE64041(?:-GPL6244)?_series_matrix\.txt\.gz",v)]
    require(len(files)==1,"exact series matrix unavailable")
    raw,src=sized_gzip(folder+files[0],180_000_000)
    # Metadata prefix is parsed from the verified full gzip, not a separate request.
    with gzip.GzipFile(fileobj=io.BytesIO(raw)) as f:
        lines=[];metadata_bytes=0
        for line in f:
            t=line.decode("utf-8-sig")
            if t.startswith("!series_matrix_table_begin"):break
            lines.append(t);metadata_bytes+=len(line)
            require(metadata_bytes<8_000_000,"metadata cap")
        else:raise ValueError("series matrix table missing")
    meta=base.meta_rows("".join(lines));ids=one_field(meta,"!Sample_geo_accession",125)
    titles=one_field(meta,"!Sample_title",125);tissue=characteristic(meta,"tissue",125)
    disease=characteristic(meta,"disease state",125)
    plat=one_field(meta,"!Sample_platform_id",125)
    require(set(plat)=={"GPL6244"} and len(set(ids))==125,"sample/platform identity")
    pairs={};normal=0
    for i,(title,t,d) in enumerate(zip(titles,tissue,disease)):
        m=re.fullmatch(r"patient (\d+) (tumor|non-tumor)",title.strip(),re.I)
        if m:
            group=m[2].lower();require(t.lower()==("tumor" if group=="tumor" else "liver parenchyma") and d.lower()=="hepatocellular carcinoma patient","title/tissue/disease disagreement")
            key=(m[1],group);require(key not in pairs,"duplicate pair site");pairs[key]=ids[i]
        else:
            # Five healthy references are explicitly withheld, never treated as adjacent tissue.
            require(re.fullmatch(r"normal[1-5]",title.strip(),re.I) and t.lower()=="liver" and d.lower()=="normal","unresolved healthy specimen");normal+=1
    patients=sorted({p for p,g in pairs});require(len(patients)==60 and len(pairs)==120 and normal==5 and all((p,g) in pairs for p in patients for g in ("tumor","non-tumor")),"60 pair/5 healthy denominator drift")
    probes,ann,platform_n=platform()
    h,rows=decode_table(raw,"!series_matrix_table_begin","!series_matrix_table_end")
    require(h[0]=="ID_REF" and len(h)==126 and set(h[1:])==set(ids),"GSM matrix columns mismatch")
    index={v:i for i,v in enumerate(h)}
    selected=set(sum(probes.values(),[]));vectors={};seen=set()
    for row in rows:
        k=row[0];require(k not in seen,"duplicate matrix feature");seen.add(k)
        if k not in selected:continue
        v=np.asarray(row[1:],dtype=float);require(np.isfinite(v).all(),"nonfinite probe value")
        vectors[k]=v
    require(len(seen)==33297,"reported matrix feature count drift")
    results={}
    for gene in PANEL:
        available=[k for k in probes[gene] if k in vectors]
        if not available or len(available)!=len(probes[gene]):
            results[gene]={"status":"HOLD_INCOMPLETE_EXACT_PROBES","mapped":len(available),"annotated":len(probes[gene])};continue
        # Freeze equal-probe average, not strongest differential probe. Sensitivity: all probe directions.
        ds=[]
        for k in available:
            v=vectors[k];ds.append([float(v[index[pairs[p,"tumor"]]-1]-v[index[pairs[p,"non-tumor"]]-1]) for p in patients])
        delta=np.mean(ds,axis=0).tolist();z=sign_summary(delta,"publisher log2 RMA tumor minus non-tumor")
        z.update(status="EXTERNAL_PAIRED_RNA_CONTEXT",probes=len(available),probe_median_range=[min(map(statistics.median,ds)),max(map(statistics.median,ds))])
        results[gene]=z
    return {"status":"HCC_SWISS_PAIRED_FIXED_PANEL","source":src,"directory":listing,"annotation":ann,"platform_rows":platform_n,"feature_rows":len(seen),"pairs":60,"healthy_withheld":5,"results":results,
            "limits":["Independent geographic cohort relative to Jiang HBV cohort; different etiology/platform/stage composition. Fixed eight genes, no post-result selection. GPL6244 build35 publisher annotation; equal uniquely assigned probe mean. RMA and FPKM effects not directly comparable. No claims of novel mechanism or clinical utility."]}
def korea():
    meta,ids,titles,src=base.geo_metadata("GSE36376")
    require(len(ids)==433 and len(set(ids))==433,"433 sample denominator")
    tissue=characteristic(meta,"tissue",433);counts=dict(collections.Counter(v.lower() for v in tissue))
    return {"status":"KOREAN_HCC_METADATA_QUALIFICATION","source":src,"samples":433,"tissue_counts":counts,"characteristic_fields":base.field_profile(meta,ids),
            "limits":["Samsung Medical Center GEO GSE36376; 240 tumor/193 adjacent reported, check actual count. Illumina chip identifiers and SAMPLE numbers do not prove patient linkage. No unpaired significance tests with possible within-patient overlap."]}
def protein_source():
    import openpyxl
    # Exact publisher link for Extended Data Figure 9 (protein and mRNA cholesterol-metabolism plots).
    url="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-019-0987-8/MediaObjects/41586_2019_987_MOESM17_ESM.xlsx"
    raw,src=base.get(url,30_000_000)
    wb=openpyxl.load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
    inventory=[]
    for i,ws in enumerate(wb,1):
        labels=collections.Counter();gene_counts=collections.Counter();headers=[]
        for row in ws.iter_rows(min_row=1,max_row=min(ws.max_row or 1,20),values_only=True):
            terms=[]
            for v in row:
                if not isinstance(v,str) or len(v)>100:continue
                s=v.strip()
                if s in PANEL:gene_counts[s]+=1
                # Record only known analyte/site/measurement concepts, no source IDs.
                for k in ("protein","mRNA","iBAQ","FPKM","intensity","abundance","tumor","tumour","non-tumor","non-tumour","cholesterol"):
                    if re.search(r"\b"+re.escape(k)+r"\b",s,re.I):labels[k]+=1
            for v in row:
                if isinstance(v,str) and v.strip() in PANEL:terms.append(v.strip())
            if terms:headers.append({"fixed_gene_labels":terms})
        inventory.append({"sheet_index":i,"dimensions":[ws.max_row,ws.max_column],"known_terms":dict(labels),"fixed_gene_labels":dict(gene_counts),"gene_header_rows":headers})
    wb.close()
    return {"status":"HCC_EXTENDED_FIG9_SOURCE_READ","source":src,"sheets":inventory,
            "limits":["Exact original Extended Data Fig9 source; already published cholesterol protein/mRNA comparisons. Gene labels and measurement schema only; no inferred paired identifiers or reuse as independent validation."]}
def supplement_effects():
    import openpyxl
    p=prior("hcc_protein")["source"];raw,src=base.get(p["url"],p["bytes"]);require(src["sha256"]==p["sha256"],"supplement archive drift")
    out=[]
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        for n in (6,7,9):
            hits=[f for f in z.infolist() if Path(f.filename).name==f"2017-04-04932E-Supplementary Table {n}.xlsx"]
            require(len(hits)==1,"exact supplement member")
            w=openpyxl.load_workbook(io.BytesIO(z.read(hits[0])),read_only=True,data_only=True);s=[]
            for i,ws in enumerate(w,1):
                head=[]
                for j,row in enumerate(ws.iter_rows(min_row=1,max_row=min(ws.max_row,8),values_only=True),1):
                    # Explicit public gene-statistic headers; suppress unknown cells and all patient values.
                    h=[str(v).strip() if v is not None else "" for v in row]
                    ix=[k for k,v in enumerate(h) if re.fullmatch(r"gene\s*symbol|symbol",v,re.I)]
                    if len(ix)!=1:continue
                    stats=[{"column":k+1,"label":v} for k,v in enumerate(h) if len(v)<=100 and re.fullmatch(r"(?:log2?\s*)?(?:fold.change|fc|ratio)(?:\s*\([^)]{1,60}\))?|(?:adjusted\s*)?p[- ]?value|fdr|q[- ]?value|uniprot\s*id|gene\s*symbol",v,re.I)]
                    matches=collections.Counter()
                    for rr in ws.iter_rows(min_row=j+1,values_only=True):
                        if len(rr)>ix[0] and isinstance(rr[ix[0]],str) and rr[ix[0]].strip() in PANEL:matches[rr[ix[0]].strip()]+=1
                    head.append({"header_row":j,"statistic_labels":stats,"fixed_gene_row_counts":dict(matches)})
                s.append({"sheet":i,"dimensions":[ws.max_row,ws.max_column],"profiles":head})
            w.close();out.append({"table":n,"sheets":s})
    return {"status":"HCC_PUBLIC_PROTEIN_EFFECT_TABLE_GATE","source":src,"tables":out,
            "limits":["Original publisher tables only; no inferred signed logFC or direction. Missing fixed genes may reflect selected lists, not no effect. Do not treat original effect tables as independent validation."]}
def selftest():
    assert len(PANEL)==8
    assert sign_summary([1]*5,"test")["p_exact_sign"]==.0625
    assert sign_summary([1]*5,"test")["bonferroni_8"]==.5
    assert prior("hcc_rna")["pairs"]==35 and prior("sc_headers")["verified"]==20
    assert prior("sc_pilot")["full_pair_pass"] is True
    print("Focused frozen-family, exact sign, and source archive assertions passed.")
def protein_effect_readout():
    """Original selected gene list only; direction and analyte identity must be explicitly documented."""
    import openpyxl
    arch=json.loads((ROOT/"docs/publication/evidence/early_track_external_first_readback_20261125.json").read_text())
    require(arch["run_id"]==45607898177,"first external run mismatch")
    p=next(x["aggregate"] for x in arch["records"] if x["aggregate"].get("mode")=="supplement_effects")["source"]
    raw,src=base.get(p["url"],p["bytes"]);require(src["sha256"]==p["sha256"],"publisher ZIP drift")
    result={}
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        members=[f for f in z.infolist() if Path(f.filename).name=="2017-04-04932E-Supplementary Table 6.xlsx"]
        require(len(members)==1,"exact Table 6 missing")
        wb=openpyxl.load_workbook(io.BytesIO(z.read(members[0])),read_only=True,data_only=True)
        require(len(wb.worksheets)==2 and wb.worksheets[0].max_row==2009,"Table 6 shape drift")
        # Second sheet is a small glossary; emit only source-authored definition fields
        # that explicitly refer to the known columns, not other records or clinical data.
        glossary={}
        for rr in wb.worksheets[1].iter_rows(values_only=True):
            if len(rr)<2 or not isinstance(rr[0],str) or not isinstance(rr[1],str):continue
            key=rr[0].strip().lower()
            if key in {"uniprot id","gene symbol","logfc","p-value","p value","fdr","adjusted p-value","regulation","direction"}:
                glossary[key]=rr[1][:200]
        ws=wb.worksheets[0];it=ws.iter_rows(values_only=True);h=[str(x).strip() for x in next(it)]
        require(h[:4]==["UniProt ID","Gene symbol","logFC","P-value"] or [x.lower() for x in h[:4]]==["uniprot id","gene symbol","logfc","p-value"],"exact source effect header mismatch")
        rows=collections.defaultdict(list);n=0
        for rr in it:
            n+=1;g=rr[1]
            if g not in PANEL:continue
            a,b=rr[2:4]
            require(isinstance(a,(int,float)) and isinstance(b,(int,float)) and math.isfinite(a) and math.isfinite(b) and 0<=b<=1,"invalid published gene effect")
            rows[g].append({"published_logFC":float(a),"published_p_value":float(b)})
        require(n==2008,"table row count drift")
        for g in PANEL:
            q=rows[g]
            result[g]={"status":"PUBLISHED_TABLE_EFFECT_ONLY",**q[0]} if len(q)==1 else {"status":"HOLD_ABSENT_OR_NONUNIQUE_IN_SELECTED_TABLE","rows":len(q)}
        wb.close()
    return {"status":"HCC_ORIGINAL_EFFECT_READOUT","source":src,"table":6,"rows":n,"glossary":glossary,"results":result,
            "limits":["Selected original published table, no independent protein validation. Published logFC and p-values exactly as stored, not re-estimated; no inference about numerator, denominator, paired sample identity, or missing candidates. Interpret only if source legend and article confirm direction/analyte. No clinical or patient-level values."]}
def source_repair_selftest():
    selftest()
    import tempfile
    a=json.loads((ROOT/"docs/publication/evidence/early_track_external_first_readback_20261125.json").read_text())
    require(a["run_id"]==45607898177 and len(a["records"])==4,"archived first external run")
    # Explicit complete-pair validation and two-sided test consistency.
    require(sign_summary([1]*8,"test")["bonferroni_8"]==0.0625,"sign test correction")
    print("First external readback provenance and paired correction checks passed.")

def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["swiss","korea","protein_source","supplement_effects","protein_effect_readout","selftest"]);p.add_argument("--out",type=Path);a=p.parse_args()
    if a.mode=="selftest":source_repair_selftest();return
    try:z=globals()[a.mode]()
    except (OSError,EOFError,ValueError,KeyError,TypeError,IndexError,UnicodeError,StopIteration,zipfile.BadZipFile) as e:z={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:220]}
    z.update(mode=a.mode,scientific_status="HOLD_NOVELTY_AND_INDEPENDENT_VALIDATION")
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n");print(json.dumps(z,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
