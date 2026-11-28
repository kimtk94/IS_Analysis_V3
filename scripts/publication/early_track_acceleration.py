#!/usr/bin/env python3
"""Bounded public-source follow-up for six early tracks; no private data or automatic promotion."""
import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
import math
import re
import sys
import tempfile
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path

HOSTS={"ftp.ncbi.nlm.nih.gov","www.ncbi.nlm.nih.gov","pmc-oa-opendata.s3.amazonaws.com",
       "ftp.ebi.ac.uk","rest.ensembl.org","grch37.rest.ensembl.org","charls.pku.edu.cn"}
TRACKS=("oral","muscle","hcc","skin","metabolic","ckm")
GEO="https://ftp.ncbi.nlm.nih.gov/geo/"
S3="https://pmc-oa-opendata.s3.amazonaws.com/"
SKIN_BASE="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90094001-GCST90095000/"

def require(ok,msg):
    if not ok: raise ValueError(msg)

def get(url,cap,span=None,body=None):
    u=urllib.parse.urlsplit(url)
    require(u.scheme=="https" and u.hostname in HOSTS,"unapproved host")
    headers={"User-Agent":"IS-early-track-public-source/1.0","Accept-Encoding":"identity"}
    if span: headers["Range"]="bytes=%d-%d"%span
    if body is not None: headers.update({"Content-Type":"application/json","Accept":"application/json"})
    req=urllib.request.Request(url,data=body,headers=headers)
    with urllib.request.urlopen(req,timeout=45) as f:
        require(urllib.parse.urlsplit(f.geturl()).hostname==u.hostname,"redirect host changed")
        require(f.status==(206 if span else 200),"unexpected status")
        cr=f.headers.get("Content-Range")
        if span:
            m=re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",cr or "")
            require(m and (int(m[1]),int(m[2]))==span and int(m[3])>span[1],"incorrect byte range")
        raw=f.read(cap+1)
        require(0<len(raw)<=cap,"empty/oversize transfer")
        if span: require(len(raw)==span[1]-span[0]+1,"truncated range")
        return raw,{"url":url,"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest(),"content_range":cr,"etag":f.headers.get("ETag")}

def gzip_prefix(raw,cap=16_000_000):
    out=zlib.decompressobj(16+zlib.MAX_WBITS).decompress(raw,cap+1)
    require(len(out)<=cap,"gzip prefix expanded beyond cap")
    return out.decode("utf-8","strict")

def geo_folder(acc):
    require(re.fullmatch(r"GSE\d+",acc),"bad accession")
    return GEO+"series/"+acc[:-3]+"nnn/"+acc+"/"

def links(raw):
    return sorted(set(urllib.parse.unquote(a) for a in re.findall(r"""href\s*=\s*["']([^"']+)["']""",raw.decode("utf-8","replace"),re.I)))

def meta_rows(txt):
    # Only quoted GEO metadata fields; no expression rows are processed.
    rows=collections.defaultdict(list)
    for line in txt.splitlines():
        if line.startswith("!series_matrix_table_begin"): break
        if not line.startswith(("!Sample_","!Series_")): continue
        a=next(csv.reader([line],delimiter="\t"))
        rows[a[0]].append(a[1:])
    return rows

def geo_metadata(acc):
    folder=geo_folder(acc)+"matrix/"
    raw,listing=get(folder,2_000_000)
    files=[x for x in links(raw) if re.fullmatch(re.escape(acc)+r"(?:-[A-Za-z0-9_]+)?_series_matrix\.txt\.gz",x)]
    require(len(files)==1,"one exact GEO series matrix not established")
    url=folder+files[0]
    try:
        raw,src=get(url,4_000_000)
        txt=gzip_prefix(raw)
        mode="bounded_complete_gzip"
    except (OSError,ValueError):
        raw,src=get(url,1_048_576,(0,1_048_575))
        txt=gzip_prefix(raw)
        mode="bounded_compressed_prefix"
    require("!series_matrix_table_begin" in txt,"metadata header incomplete")
    rows=meta_rows(txt)
    ids=rows.get("!Sample_geo_accession",[[]])[0]
    require(ids and len(ids)==len(set(ids)) and all(re.fullmatch("GSM[0-9]+",s) for s in ids),"invalid sample IDs")
    title=rows.get("!Sample_title",[[]])[0]
    require(len(title)==len(ids),"sample title mismatch")
    series_ids=[v for rr in rows.get("!Series_geo_accession",[]) for v in rr]
    require(acc in series_ids,"series accession mismatch")
    return rows,ids,title,{"directory":listing,"matrix":src,"mode":mode}

def field_profile(rows,ids):
    n=len(ids); result=collections.defaultdict(lambda:{"covered":0,"numeric":0,"unique":set()})
    for rr in rows.get("!Sample_characteristics_ch1",[]):
        require(len(rr)==n,"misaligned GEO characteristic")
        for i,x in enumerate(rr):
            k,sep,v=x.partition(":")
            if not sep or not v.strip():continue
            k=k.strip().lower();v=v.strip()
            z=result[k];z["covered"]+=1;z["unique"].add(v)
            if re.fullmatch(r"[+-]?\d+(?:\.\d+)?",v):z["numeric"]+=1
    return {k:{"covered":v["covered"],"numeric":v["numeric"],"distinct_values":len(v["unique"])} for k,v in result.items()}

def sample_files(rows,ids):
    out={i:{} for i in ids}
    for k,rrs in rows.items():
        if not k.startswith("!Sample_supplementary_file"):continue
        for rr in rrs:
            require(len(rr)==len(ids),"supplementary file metadata alignment")
            for sid,v in zip(ids,rr):
                u=urllib.parse.urlsplit(v.replace("ftp://","https://",1))
                if u.hostname!="ftp.ncbi.nlm.nih.gov" or "/"+sid+"/" not in u.path:continue
                for name,pat in (("features",r"features\.tsv\.gz$"),("barcodes",r"barcodes\.tsv\.gz$"),("matrix",r"matrix\.mtx\.gz$")):
                    if re.search(pat,u.path,re.I):
                        require(name not in out[sid],"duplicate file type")
                        out[sid][name]=urllib.parse.urlunsplit(u)
    return out

def matrix_dimensions(txt):
    a=txt.splitlines()
    require(a and a[0].lower().startswith("%%matrixmarket matrix coordinate integer"),"not integer coordinate matrix")
    for line in a[1:]:
        if line.startswith("%"):continue
        d=line.split()
        require(len(d)==3 and all(v.isdigit() for v in d),"invalid matrix dimensions")
        g,c,nnz=map(int,d);require(g>0 and c>0 and 0<nnz<=g*c,"invalid matrix size")
        return g,c,nnz
    raise ValueError("matrix dimensions missing")

def oral():
    rows,ids,titles,src=geo_metadata("GSE322620")
    require(len(ids)==10,"expected ten original specimens")
    groups=collections.Counter()
    for title in titles:
        m=re.fullmatch(r"(PT|MET)[ _-]*[1-5]",title.strip(),re.I)
        groups[m.group(1).upper() if m else "unresolved"]+=1
    # PT/MET names do not establish pairing; original summary says five per group.
    files=sample_files(rows,ids)
    records=[]
    for sid,title in zip(ids,titles):
        z={"status":"HOLD_FILES","group":"PT" if re.match(r"^PT",title,re.I) else "MET" if re.match(r"^MET",title,re.I) else "unresolved","file_types":sorted(files[sid])}
        try:
            require(set(files[sid])=={"features","barcodes","matrix"},"missing exact official 10x triplet")
            feat,fm=get(files[sid]["features"],3_000_000)
            bar,bm=get(files[sid]["barcodes"],3_000_000)
            mat,mm=get(files[sid]["matrix"],65536,(0,65535))
            frows=[v.split("\t") for v in gzip.decompress(feat).decode().splitlines() if v.strip()]
            brows=[v for v in gzip.decompress(bar).decode().splitlines() if v.strip()]
            d=matrix_dimensions(gzip_prefix(mat))
            require(d[:2]==(len(frows),len(brows)) and all(len(v)>=2 for v in frows),"feature/barcode dimension mismatch")
            require(len({v[0] for v in frows})==len(frows) and len(set(brows))==len(brows),"duplicate feature/barcode")
            # Hash public files, never emit barcode or cell-level data; no matrix-body QC.
            z.update(status="DIMENSION_PREFIX_VERIFIED",genes=d[0],cells=d[1],declared_nonzero=d[2],
                     feature_sha256=fm["sha256"],barcode_sha256=bm["sha256"],matrix_prefix_sha256=mm["sha256"])
        except (OSError,ValueError,KeyError,IndexError,UnicodeError) as e:z["reason"]=type(e).__name__
        records.append(z);time.sleep(1)
    return {"status":"ORAL_NEW_COHORT_QUALIFICATION","source":src,"specimens":len(ids),"title_group_counts":dict(groups),
            "characteristic_fields":field_profile(rows,ids),"matrix_prefix_verified":sum(x["status"]=="DIMENSION_PREFIX_VERIFIED" for x in records),
            "records":records,"limits":["Unpaired unless documented same-donor identity, no inference from PT/MET number; no full matrix body checksum or cell type assignment; no DEG or patient-level inference."]}

def muscle():
    results={}
    for acc in ("GSE270823","GSE163434","GSE157585"):
        try:
            rows,ids,titles,src=geo_metadata(acc)
            fields=field_profile(rows,ids)
            endpoint={k:v for k,v in fields.items() if re.search(r"lean mass|hypertroph|cross.section|qCSA|muscle (?:size|mass)|change in|response",k,re.I)}
            identity={k:v for k,v in fields.items() if re.search(r"subject|participant|patient|donor|study|time.point|limb|leg",k,re.I)}
            results[acc]={"status":"METADATA_READ","source":src,"samples":len(ids),"fields":fields,
                          "endpoint_candidates":endpoint,"linkage_candidates":identity}
        except (OSError,ValueError,KeyError,IndexError,UnicodeError) as e:results[acc]={"status":"HOLD_SOURCE","reason":str(e)[:150]}
    return {"status":"MUSCLE_NEW_ENDPOINT_SOURCE_AUDIT","sources":results,
            "limits":["No numeric field is automatically a validated outcome; confirm participant/limb, units, baseline/post, overlap and selection. GSE163434 reuses GSE157585; GSE270823 has seven shared baseline specimens. Preserve #343 negative result, no panel retuning or prediction."]}

def hcc():
    # Public OA supplement inventory; do not access EGA or infer sample matches.
    import openpyxl
    raw,listing=get(S3+"?list-type=2&prefix=PMC9068765.&delimiter=%2F&max-keys=20",200000)
    root=ET.fromstring(raw)
    prefixes=[e.text for e in root.findall(".//{*}CommonPrefixes/{*}Prefix") if re.fullmatch(r"PMC9068765\.\d+/",e.text or "")]
    require(prefixes,"no OA article version")
    ver=max(prefixes,key=lambda x:int(x.split(".")[1].rstrip("/")))
    raw,meta=get(S3+"metadata/"+ver.rstrip("/")+".json",300000)
    obj=json.loads(raw)
    require(str(obj.get("pmcid"))=="PMC9068765" and obj.get("doi")=="10.1038/s41467-022-29960-8" and obj.get("is_pmc_openaccess") is True,"article identity/OA gate")
    raw,xsrc=get(S3+ver+ver.rstrip("/")+".xml",2_000_000)
    xml=ET.fromstring(raw)
    require("10.1038/s41467-022-29960-8" in raw.decode(),"XML DOI mismatch")
    names={}
    for node in xml.iter():
        if node.tag.split("}")[-1] not in ("supplementary-material","media"):continue
        caption=" ".join(" ".join(n.itertext()) for n in node if n.tag.split("}")[-1] in ("caption","label"))
        for v in node.iter():
            href=v.get("{http://www.w3.org/1999/xlink}href")
            if href and re.fullmatch(r"[A-Za-z0-9_.-]+\.xlsx",href):
                names[href]=caption[:200]
    raw,objs=get(S3+"?list-type=2&prefix="+ver+"&max-keys=100",500000)
    lr=ET.fromstring(raw)
    sizes={x.findtext("{*}Key"):int(x.findtext("{*}Size")) for x in lr.findall(".//{*}Contents")}
    inventory=[{"file":n,"caption":v,"bytes":sizes.get(ver+n)} for n,v in names.items()]
    inspected=[]
    # Exact labeled Data 1 first, not guessed MOESM file number.
    for n,caption in names.items():
        if not re.search(r"\bSupplementary Data\s*1\b",caption,re.I):continue
        size=sizes.get(ver+n)
        if size is None or not 0<size<=12_000_000:continue
        b,m=get(S3+ver+n,size)
        require(len(b)==size,"workbook size mismatch")
        etag=(m.get("etag") or "").strip('"')
        require(re.fullmatch("[0-9a-f]{32}",etag) and hashlib.md5(b).hexdigest()==etag,"simple object MD5 unavailable/mismatch")
        wb=openpyxl.load_workbook(io.BytesIO(b),read_only=True,data_only=True)
        schemas=[]
        for ws in wb:
            labels=collections.Counter()
            for row in ws.iter_rows(min_row=1,max_row=min(ws.max_row or 1,8),values_only=True):
                for v in row:
                    if not isinstance(v,str):continue
                    for concept in ("patient","sample","tumor","tumour","RNA","proteom","survival","follow","BCLC","stage","etiology","grade","event"):
                        if re.search(concept,v,re.I):labels[concept]+=1
            schemas.append({"rows":ws.max_row,"columns":ws.max_column,"recognized_header_terms":dict(labels)})
        wb.close()
        inspected.append({"file":n,"bytes":size,"sha256":m["sha256"],"sheet_schemas":schemas})
        break
    return {"status":"HCC_ALTERNATE_PUBLIC_SUPPLEMENT_AUDIT","article":meta,"xml":xsrc,"inventory":inventory,"inspected":inspected,
            "limits":["122 tumors/114 patients and 51 multi-omic tumors/49 patients reported by paper; not observed matched rows. EGA raw sequencing restricted. No individual records emitted, no RNA-protein join or survival model."]}

def skin():
    import yaml
    out={}
    for acc in ("GCST90094903","GCST90094904"):
        z={"status":"HOLD_SOURCE"}
        try:
            folder=SKIN_BASE+acc+"/harmonised/"
            raw,listing=get(folder,2_000_000)
            names=[n for n in links(raw) if re.fullmatch(re.escape(acc)+r"[^/]*\.h\.tsv\.gz",n)]
            require(len(names)==1,"one harmonised file not established")
            url=folder+names[0]
            raw,mdsrc=get(url+"-meta.yaml",300000)
            md=yaml.safe_load(raw)
            require(md.get("gwas_id")==acc and md.get("data_file_name")==names[0] and md.get("genome_assembly")=="GRCh38" and md.get("is_harmonised") is True,"harmonised metadata mismatch")
            raw,idxsrc=get(url+".tbi",8_000_000)
            require(gzip.decompress(raw).startswith(b"TBI\x01"),"TBI magic")
            raw,hsrc=get(url,65536,(0,65535))
            header=gzip_prefix(raw).splitlines()[0].lstrip("#").split("\t")
            require(len(header)==len(set(header)),"duplicate header")
            required=[("hm_chrom","chromosome"),("hm_pos","base_pair_location"),("hm_effect_allele","effect_allele"),
                      ("hm_other_allele","other_allele"),("hm_beta","beta"),("standard_error","hm_standard_error"),
                      ("p_value","hm_p_value"),("hm_rsid","rsid","rs_id"),("hm_code","harmonisation_code")]
            require(all(any(v in header for v in alternatives) for alternatives in required),"missing allele/effect/index fields")
            z.update(status="INDEXED_SOURCE_PREFIX_QUALIFIED",directory=listing,metadata=mdsrc,index=idxsrc,header=hsrc,
                     columns=header,trait=md.get("trait_description"),samples=md.get("samples"))
        except (OSError,ValueError,KeyError,TypeError,IndexError) as e:z["reason"]=str(e)[:200]
        out[acc]=z
        time.sleep(3)
    return {"status":"SKIN_INDEXED_RECOVERY_GATE","studies":out,
            "limits":["Different route from failed full original TSV transfer; index magic/header only, not full-file integrity or effect extraction. No MR or novelty claim. Keep phenotypes separate."]}

def metabolic():
    raw,src=get("https://charls.pku.edu.cn/en/News.htm",2_000_000)
    text=re.sub(r"<[^>]+>"," ",raw.decode("utf-8","replace"))
    observed={k:bool(re.search(pattern,text,re.I)) for k,pattern in {
        "wave_3_blood_2015":r"Blood Data of CHARLS 2015",
        "wave_4_2018":r"Wave 4 \(2018\)",
        "wave_5_2020":r"Wave 5 \(2020\)"}.items()}
    # Two observed clinical waves can establish a transition but cannot establish relapse.
    return {"status":"METABOLIC_LONGITUDINAL_RELEASE_AUDIT","source":src,"official_release_mentions":observed,
            "endpoint_gate":"HOLD_3_COMPLETE_BIOMARKER_WAVES_NOT_ESTABLISHED",
            "limits":["2011/2015 blood panels reported in release documentation; later survey release is not equivalent to repeated fasting blood assays. No individual data, genetics, recovery/relapse counts or models. KoGES remains primary candidate for three-wave genetic design."]}

def reference_maps(panel):
    # Both assembly reference positions must be independently resolved; use no proxy.
    import publication_common_instrument_pilot as old
    refs={}
    for build,host in (("GRCh37","grch37.rest.ensembl.org"),("GRCh38","rest.ensembl.org")):
        raw,source=get("https://"+host+"/variation/homo_sapiens",10_000_000,body=json.dumps({"ids":[r["rsid"] for r in panel]}).encode())
        refs[build]=(json.loads(raw),source)
    out={}
    for r in panel:
        rs=r["rsid"];maps={}
        for build,(data,source) in refs.items():
            valid=[]
            for m in data.get(rs,{}).get("mappings",[]):
                alleles=str(m.get("allele_string","")).split("/")
                if m.get("assembly_name")!=build or str(m.get("seq_region_name"))!=str(r["chrom"]) or m.get("start")!=m.get("end") or m.get("strand")!=1 or len(alleles)!=2:continue
                if old.orientation(r["ea"],r["oa"],*alleles) is None:continue
                valid.append((int(m["start"]),alleles))
            if len(valid)==1:maps[build]=valid[0]
        if set(maps)!={"GRCh37","GRCh38"} or maps["GRCh37"][0]!=r["pos"] or maps["GRCh37"][1]!=maps["GRCh38"][1]:
            out[rs]={"status":"HOLD_REFERENCE"};continue
        pos,(ref,alt)=maps["GRCh38"]
        sign=old.orientation(r["ea"],r["oa"],alt,ref)[0]
        out[rs]={"status":"MATCHED_REFERENCE","chr":str(r["chrom"]),"pos":pos,"ref":ref,"alt":alt,"exposure_sign":sign}
    return out,{k:v[1] for k,v in refs.items()}

def ckm():
    import publication_common_instrument_pilot as old
    import fixed_subtype_indexed as idx
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)/"ldl.gz"
        source=old.download(old.LDL_URL,p,60_000_000,180)
        require(source["sha256"]=="ddf75ed591ffcce8c3e81ae2f4c77f4428cd7dbf4628b24467aa8cf6e85e3dce","GLGC 2013 source SHA drift")
        panel,selection=old.load_ldl(p)
    require(len(panel)==19 and selection["source_rows"]==2437751 and selection["eligible_genomewide_significant_rows"]==2511,"frozen exposure selection changed")
    ref,ref_sources=reference_maps(panel)
    idx.PANEL=tuple(ref)
    outcomes={}; studies={}
    for name,study,trait in (("LAA","EUR_LAS","Large artery stroke"),("CE","EUR_CES","Cardioembolic stroke"),("SVD","EUR_SVS","Small vessel stroke")):
        x=idx.extract(study,ref);md=x.get("metadata",{})
        verified=(md.get("trait_description")==[trait] and {v for row in md.get("samples",[]) for v in row.get("sample_ancestry_category",[])}=={"European"} and bool(x.get("fields",{}).get("rsid")))
        hits={}
        if verified:
            for rs,v in x["rows"].items():
                if v["status"]=="MATCHED_HARMONISED_ROW":
                    hits[rs]={"by":v["beta_alt"]*ref[rs]["exposure_sign"],"sy":v["se"]}
        outcomes[name]=hits
        studies[name]={"source_metadata_qualified":verified,"matched":len(hits),"readback":x}
    ans=old.analyze(panel,outcomes)
    return {"status":"CKM_LDL_INDEXED_PILOT","exposure_source":source,"selection":selection,"reference_sources":ref_sources,
            "reference_verified":sum(v["status"]=="MATCHED_REFERENCE" for v in ref.values()),"studies":studies,"analysis":ans,
            "limits":["Historical 19 SNP sparse panel, same source and rules as #334. Harmonised indexed ranges, not original full file checksum; no SBP retuning. Direct difference regression conditional on scales, covariance bounds and independent instruments. No new causal claim."]}

def selftest():
    assert matrix_dimensions("%%MatrixMarket matrix coordinate integer general\n% x\n3 4 5\n")== (3,4,5)
    try:matrix_dimensions("%%MatrixMarket matrix coordinate real general\n3 4 5\n")
    except ValueError:pass
    else:raise AssertionError("noninteger matrix accepted")
    x=meta_rows('!Sample_geo_accession\t"GSM1"\t"GSM2"\n!Sample_characteristics_ch1\t"lean mass: 1.2"\t"lean mass: 2"\n!series_matrix_table_begin\n')
    assert field_profile(x,["GSM1","GSM2"])["lean mass"]=={"covered":2,"numeric":2,"distinct_values":2}
    assert geo_folder("GSE322620").endswith("GSE322nnn/GSE322620/")
    assert len(TRACKS)==6
    print("Five focused source-contract checks passed; no biological result.")

def main():
    p=argparse.ArgumentParser();p.add_argument("track",choices=TRACKS+("selftest",));p.add_argument("--out",type=Path)
    a=p.parse_args()
    if a.track=="selftest":selftest();return
    require(a.out is not None,"output path required")
    try:result=globals()[a.track]()
    except (OSError,ValueError,KeyError,TypeError,IndexError,ET.ParseError,UnicodeError) as e:
        result={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:250]}
    result["track"]=a.track
    result["scientific_status"]="HOLD_NOT_INDEPENDENTLY_VALIDATED"
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(result,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
