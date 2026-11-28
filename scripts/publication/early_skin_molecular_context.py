#!/usr/bin/env python3
"""Frozen Skin nominal cis-eQTL inventory, and narrow nonidentifying source diagnostics."""
import argparse,collections,csv,gzip,hashlib,io,json,math,re,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
for d,paths in (("prior",["scripts/publication"]),("indexed",["scripts","scripts/is"]),("orthogonal",["scripts/publication"])):
    for sub in paths:sys.path.insert(0,str(ROOT/d/sub))
import early_track_acceleration as base
import fixed_subtype_indexed as idx
import early_track_orthogonal_sources as orth
require=base.require
TISSUES={"sun":"Skin_Sun_Exposed_Lower_leg","protected":"Skin_Not_Sun_Exposed_Suprapubic","fibroblast":"Cells_Cultured_fibroblasts"}
FILEBASE="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/"
READBACK=ROOT/"docs/publication/evidence/early_track_latest_readbacks_20261125.json"
CROW=ROOT/"docs/publication/evidence/skin_crows_feet_verified_readback_20261125.json"
EXPECTED={"GCST90094903":"0c2e1c6cb4411f8dca4d69fa65671e9c68711d1ac6fcc660625371af316f99af",
          "GCST90094904":"567e55db53cb9528680073dd122190528bd640a0706c286d40999f4b283e5927"}
base.HOSTS.update({"zenodo.org","media.springernature.com"})
def leads():
    a=json.loads(READBACK.read_text());require(len(a["records"])==9,"readback archive count")
    x=[v["aggregate"] for v in a["records"] if v["aggregate"].get("accession")=="GCST90094903"]
    require(len(x)==1,"under eye source missing")
    y=json.loads(CROW.read_text());require(y["job_id"]==139193308088 and y["run_id"]==45610752424,"crow source provenance")
    out={}
    for d in (x[0],y["aggregate"]):
        acc=d["accession"];q=d["full_source_readout"]
        require(q["status"]=="COMPLETE_MD5_VERIFIED_LEAD_SCAN" and q["full_source"]["sha256"]==EXPECTED[acc] and q["rows_scanned"]==4950621,"full source mismatch")
        require(q["metadata"]["genome_assembly"]=="GRCh38","source assembly")
        for vid,z in q["lead_rows"].items():
            if z["status"]!="MATCHED_FULL_SOURCE_LEAD":continue
            require(re.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2])_\d+_[ACGT]_[ACGT]",vid),"invalid source variant")
            ch,p,r,b=vid.split("_")
            require(r!=b and {r,b} not in ({"A","T"},{"C","G"}),"ambiguous allele")
            v=out.setdefault(vid,{"rsid":z["rsid"],"traits":{}})
            require(v["rsid"]==z["rsid"] and acc not in v["traits"],"source variant identity")
            v["traits"][acc]={k:z[k] for k in ("beta_alt","standard_error","p_value")}
    require(len(out)==16 and sum(len(v["traits"]) for v in out.values())==17,"frozen panel denominator")
    return out
def validate_header(h,tbi):
    required={"variant","chromosome","position","ref","alt","molecular_trait_id","gene_id","beta","se","pvalue","maf","an"}
    require(required<=set(h) and len(h)==len(set(h)),"official nominal column contract")
    require(h[tbi["seq"]-1]=="chromosome" and h[tbi["begin"]-1]=="position" and tbi["end"] in (0,tbi["begin"]),"index coordinate column contract")
    require(tbi["format"]&0x10000==0,"unexpected zero-based index")
def exact_row(row,vid):
    ch,pos,ref,alt=vid.split("_")
    require(row["variant"]=="chr"+vid and row["chromosome"].removeprefix("chr")==ch and row["position"]==pos and row["ref"]==ref and row["alt"]==alt,"variant identity")
    g=row["gene_id"];mt=row["molecular_trait_id"]
    require(re.fullmatch(r"ENSG\d+(?:\.\d+)?",g) is not None and re.fullmatch(r"ENSG\d+(?:\.\d+)?",mt) is not None and g.split(".")[0]==mt.split(".")[0],"not gene expression trait")
    b,se,p,maf,an=(float(row[k]) for k in ("beta","se","pvalue","maf","an"))
    require(all(math.isfinite(v) for v in (b,se,p,maf,an)) and se>0 and 0<=p<=1 and 0<maf<=.5 and an>=2 and an.is_integer(),"nominal statistic contract")
    return {"gene":g,"molecular_trait":mt,"beta_alt":b,"se":se,"p":p,"maf":maf,"allele_number":int(an)}
def qtl(mode):
    panel=leads();tissue=TISSUES[mode];url=FILEBASE+tissue+".tsv.gz"
    out={"status":"HOLD_NOMINAL_SOURCE","tissue":tissue,"source_url":url,"panel_variants":len(panel),"rows":{}}
    try:
        raw,ir=idx.get(url+".tbi",12_000_000);tbi=idx.parse_index(raw);out["index_sha256"]=idx.sha(raw)
        first,cr=idx.get(url,65536,(0,65535));total=int(cr.rsplit("/",1)[1]);out.update(header_sha256=idx.sha(first),remote_bytes=total)
        h=idx.parse_bgzf_range(first,0,0)["data_lines"].readline().lstrip("#").rstrip("\r\n").split("\t")
        out["header"]=h;validate_header(h,tbi)
        # Reuse same indexed range for nearby markers, never download complete multi-GB GTEx files.
        cache={}
        for vid,v in sorted(panel.items()):
            ch,p,ref,alt=vid.split("_");p=int(p);z={"rsid":v["rsid"],"status":"HOLD_QUERY"}
            try:
                span=idx.indexed_span(tbi,ch,p,p,total)
                if span is None:z.update(status="NO_INDEX_CHUNK");out["rows"][vid]=z;continue
                a,b,x,y=span
                if span not in cache:
                    time.sleep(1.0)
                    rr,rc=idx.get(url,b-a+1,(a,b));cache[span]=(rr,rc)
                    require(sum(len(t[0]) for t in cache.values())<=120_000_000,"total range byte cap")
                rr,rc=cache[span];stream=idx.parse_bgzf_range(rr,x,y)["data_lines"]
                matches=collections.defaultdict(list);rejected=0;coord=0
                for row in csv.reader(stream,delimiter="\t"):
                    if len(row)!=len(h):continue
                    d=dict(zip(h,row))
                    if d["chromosome"].removeprefix("chr")!=ch or d["position"]!=str(p):continue
                    coord+=1
                    if d["ref"]!=ref or d["alt"]!=alt:continue
                    try:q=exact_row(d,vid);matches[q["gene"]].append(q)
                    except (ValueError,KeyError,TypeError,IndexError):rejected+=1
                kept=[v[0] for v in matches.values() if len(v)==1]
                z.update(status="NOMINAL_EXACT_VARIANT_READBACK",coordinate_rows=coord,exact_gene_rows=len(kept),rejected_rows=rejected,
                         duplicate_gene_keys=sum(len(v)>1 for v in matches.values()),range=rc,range_sha256=idx.sha(rr),
                         nominal_rows=sorted(kept,key=lambda x:(x["p"],x["gene"])))
            except (OSError,EOFError,ValueError,KeyError,TypeError,IndexError,UnicodeError) as e:z["reason"]=str(e)[:160]
            out["rows"][vid]=z
        out.update(status="SKIN_GTEX_V8_NOMINAL_POINT_READBACK",queries_read=sum(v["status"]=="NOMINAL_EXACT_VARIANT_READBACK" for v in out["rows"].values()),
                   variants_with_qtl_rows=sum(v.get("exact_gene_rows",0)>0 for v in out["rows"].values()),nominal_tests=sum(v.get("exact_gene_rows",0) for v in out["rows"].values()))
    except (OSError,EOFError,ValueError,KeyError,TypeError,IndexError,UnicodeError) as e:out["reason"]=str(e)[:200]
    out["limits"]=["Official eQTL Catalogue GTEx v8 imported nominal gene-expression statistics; point-indexed ranges, not full-file checksum. Exact GRCh38 ref/alt match, not LD proxies. GTEx mostly European American; no Korean-specific tissue QTL validation. Three tissue contexts share donors, not independent. No MR, conditional effects, colocalization or tissue-specific interaction inference."]
    return out
def combine(root):
    panel=leads();out={};tests=0;observed=0;missing=0;by=collections.defaultdict(dict)
    for mode,t in TISSUES.items():
        p=root/(mode+".json");d=json.loads(p.read_text()) if p.exists() else {}
        require(not d or d.get("tissue")==t,"tissue artifact identity")
        if d.get("status")=="SKIN_GTEX_V8_NOMINAL_POINT_READBACK":
            require(set(d["rows"])==set(panel),"panel artifact drift")
        out[mode]={"status":d.get("status","HOLD_MISSING"),"queries_read":d.get("queries_read",0),"nominal_tests":d.get("nominal_tests",0),"variants_with_rows":d.get("variants_with_qtl_rows",0)}
        for vid in panel:
            z=d.get("rows",{}).get(vid,{})
            if z.get("status")!="NOMINAL_EXACT_VARIANT_READBACK":missing+=1;continue
            observed+=1
            for q in z["nominal_rows"]:
                tests+=1;k=(vid,q["gene"].split(".")[0]);require(mode not in by[k],"duplicate stable gene key")
                by[k][mode]=q
    # Family is all observed exact nominal gene tests across 3 contexts, not just selected hits.
    # With incomplete queries, adjusted values are NOT computed; strong nominal flag only.
    full=(missing==0)
    rows=[]
    for (vid,g),a in sorted(by.items()):
        if not any(z["p"]<=5e-8 for z in a.values()):continue
        q={k:{**v,**({"bonferroni_observed_family":min(1,tests*v["p"])} if full else {})} for k,v in a.items()}
        rows.append({"variant":vid,"rsid":panel[vid]["rsid"],"gene":g,"GWAS":panel[vid]["traits"],"tissues":q})
    return {"status":"SKIN_FIXED_LEAD_MOLECULAR_CONTEXT","tissues":out,"unique_variants":16,"target_variant_tissue_queries":48,"completed_queries":observed,"missing_queries":missing,
            "observed_nominal_tests":tests,"full_query_family":full,"strong_nominal_pairs":len(rows),"strong_nominal_rows":rows,
            "limits":["Preselected nominal P<=5e-8 display threshold, not a colocalization finding. Bonferroni observed test family only if all 48 queries completed, not all genome-wide tests. Same variant-associated expression and wrinkle association does not prove shared causal mechanism. Full regional GWAS/QTL overlap, credible set and ancestry-aware LD needed. Pigmentation loci OCA2 rs74653330 and near KITLG rs72620727 known prior; examine pleiotropy not new gene discovery."]}
def shape(v):
    if v is None:return "empty"
    s=str(v).strip()
    if re.fullmatch(r"ENSG\d+(?:\.\d+)?",s):return "ensembl_gene"
    if re.fullmatch(r"rs\d+",s):return "rsid"
    if re.fullmatch(r"[ACGT]{14,20}(?:-\d+)?(?:_\d+)?",s):return "10x_barcode_like"
    if re.fullmatch(r"[A-Za-z][A-Za-z0-9.-]{0,19}",s):return "short_alphanumeric"
    if re.fullmatch(r"\d+",s):return "integer_like"
    return "other"
def oral_labels():
    result=[]
    # Diagnose only metadata and barcodes. No matrix redownload.
    for lib,sid in enumerate(orth.GSMS,1):
        z={"gsm":sid}
        try:
            folder,names,_=orth.filelist(sid);src={}
            for k,p,cap in (("barcode",r"barcodes\.tsv\.gz$",15_000_000),("label",r"celltype_metadata.*\.csv\.gz$",5_000_000)):
                data,src[k]=base.get(folder+orth.one(names,p),cap)
                if k=="barcode":bars=[v[0] for v in csv.reader(io.StringIO(orth.text_gz(data)),delimiter="\t") if v]
                else:rr=list(csv.reader(io.StringIO(orth.text_gz(data))))
            require(rr and all(len(v)==3 for v in rr) and len(set(bars))==len(bars),"prior three-column/unique barcode contract")
            b=set(bars);profiles=[]
            # Exact raw key, documented merge suffix; additionally report core barcode diagnostic
            # overlap only (not an authorized join) to identify the actual mismatch.
            core=collections.Counter(x.split("-")[0] for x in bars if re.fullmatch(r"[ACGT]{14,20}-\d+",x))
            for i in range(3):
                v=[row[i] for row in rr]
                shapes=collections.Counter(shape(x) for x in v)
                exact=sum(x in b for x in v);merge=sum(x.endswith("_"+str(lib)) and x.rsplit("_",1)[0] in b for x in v)
                cores=[re.match(r"^([ACGT]{14,20})(?=[^ACGT]|$)",x) for x in v]
                core_matches=sum(m is not None and core[m[1]]==1 for m in cores)
                cats=collections.Counter(v)
                # Only exact known safe cell-type names; no arbitrary cell or patient keys.
                safe={k:n for k,n in cats.items() if re.fullmatch(r"(?:Macro(?:phages?)?|Monocytes?|Myeloid|B[ -]?cells?|T[ -]?cells?|NK[ -]?cells?|Fibroblasts?|Endothelial(?: cells)?|Epithelial(?: cells)?|Mast(?: cells)?|DC|Dendritic(?: cells)?|Tumou?r(?: cells)?|Cancer(?: cells)?|Keratinocytes?|Unknown|Unassigned|Doublet|Negative|cell|barcode|celltype(?:_classification)?)",k,re.I)}
                profiles.append({"column":i+1,"unique_values":len(cats),"shape_counts":dict(shapes),"exact_barcode_matches":exact,"author_merge_suffix_matches":merge,"unique_core_overlap_diagnostic_only":core_matches,"recognized_cell_labels":safe})
            z.update(status="ACTUAL_LABEL_COLUMN_DIAGNOSTICS",sources=src,rows=len(rr),raw_barcode_columns=len(bars),columns=profiles)
        except (OSError,EOFError,ValueError,KeyError,TypeError,IndexError,UnicodeError) as e:z.update(status="HOLD_SOURCE",reason=str(e)[:180])
        result.append(z);time.sleep(.3)
    return {"status":"ORAL_LABEL_COLUMN_DIAGNOSTICS","records":result,"limits":["No cell identifiers emitted. Core overlap is diagnostic only, not accepted matching; no transformation without author evidence. No HTO assignment or expression test; 3 patients maximum with both LN types."]}
def hcc_schema():
    import openpyxl
    out={}
    for name in orth.ZEN:
        raw,src=orth.zen_download(name);w=openpyxl.load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
        sheets=[]
        for wi,ws in enumerate(w,1):
            if ws.max_row<10:
                # Small publisher notes: count explicit measurement words without outputting arbitrary values.
                concepts=collections.Counter()
                for row in ws.iter_rows(values_only=True):
                    for v in row:
                        if isinstance(v,str):
                            for term in ("RNA","protein","transcript","FPKM","TPM","log2","normalized","normalised","LFQ","iBAQ","tumor","tumour","normal","unifocal","multifocal","non-tumor","non-tumour","gene symbol","UniProt","Ensembl"):
                                if re.search(r"(?<!\w)"+re.escape(term)+r"(?!\w)",v,re.I):concepts[term]+=1
                sheets.append({"sheet":wi,"dimensions":[ws.max_row,ws.max_column],"note_concepts":dict(concepts)});continue
            header=next(ws.iter_rows(min_row=1,max_row=1,values_only=True));first_col=[]
            for row in ws.iter_rows(min_row=2,min_col=1,max_col=1,values_only=True):first_col.append(row[0])
            counts=collections.Counter(v.strip() for v in first_col if isinstance(v,str) and v.strip() in orth.PANEL)
            sheets.append({"sheet":wi,"dimensions":[ws.max_row,ws.max_column],"first_header_is_blank":header[0] is None or header[0]=="",
                           "first_column_shapes":dict(collections.Counter(shape(v) for v in first_col)),"first_column_unique":len(set(first_col)),
                           "fixed_gene_exact_label_counts":{g:counts[g] for g in orth.PANEL}})
        w.close();out[name]={"source":src,"sheets":sheets}
    supurl="https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12943-026-02604-x/MediaObjects/12943_2026_2604_MOESM3_ESM.xlsx"
    raw,src=base.get(supurl,15_000_000);require(hashlib.sha256(raw).hexdigest()=="1b7aa84bc224adb25cba9f01f628f8e6bf29b383b0891817bc1dabd555e23989","prior metadata SHA")
    w=openpyxl.load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
    s=w.worksheets[0];header=next(s.iter_rows(min_row=2,max_row=2,values_only=True))
    # Header labels only, reject any long/freeform values. This is public schema not participant rows.
    allowed=[]
    for i,v in enumerate(header,1):
        if isinstance(v,str) and len(v)<=120 and not re.search(r"[\r\n]|@|https?://",v):
            allowed.append({"column":i,"published_header":v})
    w.close()
    return {"status":"HCC_PUBLISHER_SCHEMA_DIAGNOSTICS","workbooks":out,"manifest":{"source":src,"public_headers":allowed},
            "limits":["Exact first-column fixed gene keys and public header labels only. No raw expression or participant data, no assumed tissue from sheet order, no tests. Previously 47 and 41 exact label intersections, 22 with manifest, but patient/tissue/feature match unverified."]}
def selftest():
    p=leads();assert len(p)==16 and len(p["15_27983407_C_T"]["traits"])==2
    row={"variant":"chr15_27983407_C_T","chromosome":"15","position":"27983407","ref":"C","alt":"T","molecular_trait_id":"ENSG00000000001","gene_id":"ENSG00000000001","beta":".5","se":".1","pvalue":"1e-6","maf":".1","an":"100"}
    assert exact_row(row,"15_27983407_C_T")["beta_alt"]==.5
    row["alt"]="A"
    try:exact_row(row,"15_27983407_C_T")
    except ValueError:pass
    else:raise AssertionError("wrong allele accepted")
    assert shape("ACGTACGTACGTACGT-1_2")=="10x_barcode_like"
    print("Frozen panel/source, exact allele, expression gene, and barcode schema checks passed.")
def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=[*TISSUES,"combine","oral_labels","hcc_schema","selftest"]);p.add_argument("--out",type=Path);p.add_argument("--root",type=Path,default=Path("inputs"));a=p.parse_args()
    if a.mode=="selftest":selftest();return
    try:
        z=qtl(a.mode) if a.mode in TISSUES else combine(a.root) if a.mode=="combine" else globals()[a.mode]()
    except (OSError,EOFError,ValueError,KeyError,TypeError,IndexError,UnicodeError) as e:z={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:220]}
    z.update(mode=a.mode,scientific_status="HOLD_NOVELTY_AND_CAUSAL_VALIDATION")
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n");print(json.dumps(z,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
