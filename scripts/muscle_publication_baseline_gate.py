#!/usr/bin/env python3
"""Bounded public muscle source/design audit. No phenotype prediction or individual exports."""
import argparse,collections,csv,gzip,hashlib,json,math,re,urllib.request,urllib.parse,zipfile
from pathlib import Path

SERIES={"GSE277819":81,"GSE42507":44}
FIGSHARE="29330168"
def require(ok,msg):
    if not ok:raise ValueError(msg)
def fetch(url,path,cap):
    req=urllib.request.Request(url,headers={"User-Agent":"Public-muscle-design-audit/1.0"})
    with urllib.request.urlopen(req,timeout=35) as f,path.open("wb") as w:
        h=hashlib.sha256();n=0
        while b:=f.read(1024*1024):
            n+=len(b);require(n<=cap,"Source byte cap exceeded")
            h.update(b);w.write(b)
    require(n>0,"Empty source")
    return {"url":url,"bytes":n,"sha256":h.hexdigest()}
def geo_base(acc):
    return "https://ftp.ncbi.nlm.nih.gov/geo/series/"+acc[:-3]+"nnn/"+acc
def category(v):
    # Report short categorical source labels only, not numeric trait or sample IDs.
    s=v.strip().casefold()
    if s in {"pre","post","baseline","before","after","male","female","m","f","lr","mr","hr","non","mod","xtr","1 set","4 sets","1set","4set","4sets","4 set","1-set","4-set","low","medium","high","nonresponder","responder","non-responder","low responder","medium responder","high responder","modest responder","extreme responder","low responders","medium responders","high responders"}:return s
    if re.fullmatch(r"(?:low|medium|high|non|modest|extreme)[ -]?(?:responder|response)s?",s):return s
    return "OTHER"
def phenotype_fields(k):
    s=re.sub(r"[^a-z0-9]","",k.lower())
    return bool(re.search(r"qcsa|crosssectionalarea|leanmass|leanbodymass|hypertrophy|musclethickness|myofibercsa|fiberarea|strengthchange|1rmchange|csachange|deltacsa|musclemasschange",s))
def numeric(v):
    s=v.strip().removesuffix("%").strip()
    try:return math.isfinite(float(s))
    except ValueError:return False
def geo_matrix(acc,root):
    path=root/(acc+".gz")
    src=fetch(geo_base(acc)+"/matrix/"+acc+"_series_matrix.txt.gz",path,60_000_000)
    fields=collections.defaultdict(list);series_fields={};header=None;n=0;gene_ids=set();dup=0;malformed=0;start=False;end=False
    with gzip.open(path,"rt",encoding="utf-8-sig",newline="") as f:
        for line in f:
            line=line.rstrip("\r\n")
            if line=="!series_matrix_table_begin":start=True;continue
            if line=="!series_matrix_table_end":end=True;continue
            if start and not end:
                a=next(csv.reader([line],delimiter="\t"))
                if header is None:header=a;continue
                n+=1
                require(n<=250000,"Expression row cap")
                if len(a)!=len(header):malformed+=1;continue
                if a[0] in gene_ids:dup+=1
                gene_ids.add(a[0])
                continue
            if line.startswith("!Sample_"):
                a=next(csv.reader([line],delimiter="\t"))
                fields[a[0]].append(a[1:])
            elif line.startswith(("!Series_geo_accession","!Series_platform_id","!Series_supplementary_file")):
                a=next(csv.reader([line],delimiter="\t"));series_fields.setdefault(a[0],[]).extend(a[1:])
    require(series_fields.get("!Series_geo_accession")==[acc],"Wrong GEO series")
    access=fields["!Sample_geo_accession"]
    require(len(access)==1 and len(access[0])==SERIES[acc] and len(set(access[0]))==len(access[0]),"Source accession count")
    samples=access[0];require(all(re.fullmatch(r"GSM\d+",s) for s in samples),"Invalid sample accessions")
    require(all(len(r)==len(samples) for a in fields.values() for r in a),"Metadata row width")
    if header is not None:require(header[1:]==samples,"Expression header sample order differs")
    chars=collections.defaultdict(list)
    for row in fields.get("!Sample_characteristics_ch1",[]):
        for i,v in enumerate(row):
            if ": " not in v:continue
            k,val=v.split(": ",1);chars[k.strip().casefold()].append((i,val.strip()))
    descriptions=[]
    for k,v in sorted(chars.items()):
        distinct=len({x[1] for x in v})
        d={"key":k,"nonempty":sum(bool(x[1]) for x in v),"distinct_values":distinct,
           "numeric_values":sum(numeric(x[1]) for x in v),"possible_growth_outcome_field":phenotype_fields(k)}
        if re.search(r"time|condition|group|response|treatment|set|sex",k,re.I):
            d["recognized_category_counts"]=dict(collections.Counter(category(x[1]) for x in v))
        descriptions.append(d)
    # Titles/sources are masked and reduced to counts, not exported sample labels.
    patterns={}
    for key in ("!Sample_title","!Sample_source_name_ch1"):
        if key in fields:
            patterns[key]=dict(collections.Counter(re.sub(r"\d+","#",x)[:100] for x in fields[key][0]))
    id_fields=[d for d in descriptions if re.search(r"subject|participant|donor|patient|individual|sample.?id",d["key"],re.I)]
    pairing={}
    counts_audit={}
    if acc=="GSE277819":
        require(set(("participant id","condition"))<=set(chars),"Observed participant/condition fields missing")
        by_id=dict(chars["participant id"]);by_arm=dict(chars["condition"])
        require(set(by_id)==set(by_arm)==set(range(81)),"Incomplete participant/condition metadata")
        groups=collections.defaultdict(list)
        for i in range(81):groups[by_id[i]].append(category(by_arm[i]))
        pairing={"participants":len(groups),"condition_counts":dict(collections.Counter(category(v) for v in by_arm.values())),
            "observations_per_participant":dict(collections.Counter(len(v) for v in groups.values())),
            "complete_baseline_and_two_distinct_followups":sum(sorted(v)==["1set","4set","baseline"] or sorted(v)==["1set","4sets","baseline"] or sorted(v)==["1set","4 set","baseline"] for v in groups.values())}
        require(len(groups)==27 and all(len(v)==3 for v in groups.values()),"Participant repeated-measure structure")
        try:
            file="GSE277819_raw_counts.txt.gz"
            require(file in [urllib.parse.urlsplit(u).path.rsplit("/",1)[-1] for u in series_fields.get("!Series_supplementary_file",[])],"Exact public count source not listed")
            q=root/"counts.gz";qsrc=fetch(geo_base(acc)+"/suppl/"+file,q,30000000)
            with gzip.open(q,"rt",encoding="utf-8-sig",newline="") as f:
                rd=csv.reader(f,delimiter="\t")
                h=next(rd);require(len(h)>1 and len(set(h))==len(h),"Count header schema")
                expected=set(samples);titles=set(fields["!Sample_title"][0])
                columns=h[1:]
                exact_type="GSM" if set(columns)==expected and len(columns)==81 else "title" if set(columns)==titles and len(columns)==81 else "UNRESOLVED"
                # Preserve count-header patterns, never raw individual sample IDs.
                counts_audit={"source":qsrc,"column_count":len(columns),"first_column":h[0][:100],
                  "column_patterns":dict(collections.Counter(re.sub(r"\d+","#",x)[:100] for x in columns)),"sample_identity":exact_type}
                require(exact_type!="UNRESOLVED","No exact count-column identity")
                features=set();totals=[0]*len(columns);rows=0;stable_ensembl=0;nonnegative_integer_rows=0
                for rr in rd:
                    rows+=1;require(rows<=250000 and len(rr)==len(h),"Count row bound or width")
                    require(rr[0] not in features,"Duplicate count feature")
                    features.add(rr[0]);stable_ensembl+=bool(re.fullmatch(r"ENSG\d+(?:\.\d+)?",rr[0]))
                    require(all(re.fullmatch(r"\d+(?:\.0+)?",x) for x in rr[1:]),"Noninteger count")
                    vals=[int(x.split(".")[0]) for x in rr[1:]]
                    totals=[x+y for x,y in zip(totals,vals)]
                    nonnegative_integer_rows+=1
                require(rows>1000 and all(t>0 for t in totals),"Empty gene or sample universe")
                counts_audit.update(status="COUNT_MATRIX_EXACT_COLUMN_AND_INTEGER_QC",features=rows,ensembl_features=stable_ensembl,
                  nonnegative_integer_rows=nonnegative_integer_rows,nonempty_libraries=len(totals))
        except (OSError,ValueError,EOFError,UnicodeError,csv.Error,StopIteration) as e:
            counts_audit.update(status="HOLD",reason=str(e)[:240])
        finally:(root/"counts.gz").unlink(missing_ok=True)
    return {"status":"PUBLIC_GEO_METADATA_AUDITED","series":acc,"source":src,"samples":len(samples),
       "platforms":series_fields.get("!Series_platform_id",[]),"sample_metadata_fields":sorted(fields),
       "masked_title_source_patterns":patterns,"characteristics":descriptions,"possible_id_fields":id_fields,
       "possible_growth_outcome_fields":[d for d in descriptions if d["possible_growth_outcome_field"]],
       "expression_matrix":{"table_begin":start,"table_end":end,"features":n,"duplicate_ids":dup,"malformed_rows":malformed,
                           "sample_order_matches":header is not None and header[1:]==samples},
       "series_supplementary_files":[urllib.parse.urlsplit(u).path.rsplit("/",1)[-1] for u in series_fields.get("!Series_supplementary_file",[])],"participant_pairing":pairing,"supplement_count_audit":counts_audit,
       "limits":["No individual outcomes, identifiers or expression values exported.",
                 "Exact metadata field discovery only; no inference of identity from column order or category label. A field name is not proof of a continuous outcome.",
                 "GSE277819 is outcome-selected nine participants per response group; GSE42507 is outcome-clustered. Neither can yield population-calibrated prediction without appropriate selection correction and independent test."]}
def ena(root):
    fields="run_accession,experiment_accession,sample_accession,secondary_study_accession,study_accession,experiment_title,library_name,sample_alias"
    url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1164267&result=read_run&fields="+fields+"&format=tsv"
    p=root/"ena.tsv";src=fetch(url,p,500000)
    with p.open(encoding="utf-8",newline="") as f:
        rd=csv.DictReader(f,delimiter="\t");require(set(fields.split(","))<=set(rd.fieldnames or []),"ENA header")
        rows=sorted(rd,key=lambda x:x["run_accession"])
    require(0<len(rows)<=200 and all(x["run_accession"] for x in rows),"Nonempty bounded ENA inventory")
    for k in ("run_accession","experiment_accession","sample_accession"):
        require(all(re.fullmatch(r"(?:SRR|ERR|DRR|SRX|ERX|DRX|SAMN|SAMEA|SAMD|SRS|ERS|DRS)\d+",x[k]) for x in rows),"Invalid accession")
    raw="".join("\t".join(x[k] for k in fields.split(","))+"\n" for x in rows)
    return {"status":"ENA_PUBLIC_RUN_INVENTORY","source":src,"sorted_projection_sha256":hashlib.sha256(raw.encode()).hexdigest(),
            "counts":{k:len({x[k] for x in rows}) for k in ("run_accession","experiment_accession","sample_accession")},"rows":len(rows),"study_accession_counts":dict(collections.Counter(x["study_accession"] for x in rows)),"secondary_study_counts":dict(collections.Counter(x["secondary_study_accession"] for x in rows)),"expected_81_rows":len(rows)==81,
            "masked_patterns":{k:dict(collections.Counter(re.sub(r"\d+","#",x[k])[:160] for x in rows)) for k in ("library_name","sample_alias","experiment_title")},
            "limit":"Technical sample accession is not an independent participant; verify three observations per donor, baseline and limb labels."}
def header_tokens(v):
    if not isinstance(v,str):return None
    s=v.strip()
    if len(s)>100 or not re.search("[A-Za-z]",s):return None
    # Only method/phenotype/ID column headings. Never print arbitrary row values.
    if re.search(r"qCSA|cross.sectional.area|lean.body.mass|lean.mass|muscle.mass|hypertrophy|muscle.thickness|participant|subject|sample.?id|age|sex|baseline|change|treatment|response|training.volume|log.?FC|adj.?P|FDR|p.value|gene.?id|gene.?symbol|pathway",s,re.I) and not re.search(r"\d{3,}|GSM\d+|ENSG\d+",s,re.I):
        return re.sub(r"\d+","#",s)[:100]
    return None
def xlsx_schema(path):
    from openpyxl import load_workbook
    with zipfile.ZipFile(path) as z:
        require(sum(x.file_size for x in z.infolist())<=250_000_000,"Workbook expanded size cap")
    w=load_workbook(path,read_only=True,data_only=True);out=[]
    for s in w.worksheets:
        require((s.max_row or 0)<=250000 and (s.max_column or 0)<=2000,"Workbook dimensions cap")
        # Inspect only candidate header rows; no export of actual rows.
        hits=[]
        for i,row in enumerate(s.iter_rows(min_row=1,max_row=min(8,s.max_row or 8),max_col=min(200,s.max_column or 200),values_only=True),1):
            a=sorted({x for v in row if (x:=header_tokens(v))})
            if a:hits.append({"row":i,"recognized_header_tokens":a})
        out.append({"sheet":s.title,"declared_rows":s.max_row,"declared_columns":s.max_column,"possible_schema":hits})
    w.close();return out
def supplement(root):
    p=root/"figshare.json";source=fetch("https://api.figshare.com/v2/articles/"+FIGSHARE,p,2000000)
    d=json.loads(p.read_text());require(str(d["id"])==FIGSHARE and re.fullmatch(r"10\.6084/m9\.figshare\.29330168(?:\.v\d+)?",str(d.get("doi",""))),"Exact Figshare record/DOI identity")
    require(not d.get("is_embargoed",False),"Embargoed supplementary record")
    files=[];total=0
    for x in d["files"]:
        name=x["name"];size=int(x["size"]);fid=int(x["id"])
        r={"id":fid,"name":name,"declared_bytes":size,"published_md5":x.get("computed_md5") or x.get("supplied_md5"),"status":"INVENTORY_ONLY"}
        files.append(r)
        if not name.lower().endswith(".xlsx") or size>15000000 or total+size>30000000:continue
        total+=size;path=root/(str(fid)+".xlsx")
        try:
            r["source"]=fetch("https://ndownloader.figshare.com/files/"+str(fid),path,15000000)
            require(path.stat().st_size==size,"Figshare declared size mismatch")
            with path.open("rb") as f:md5=hashlib.file_digest(f,"md5").hexdigest()
            require(bool(r["published_md5"]) and md5==r["published_md5"],"Figshare MD5 mismatch/unavailable")
            r["workbook_schema"]=xlsx_schema(path);r["status"]="PUBLIC_WORKBOOK_SCHEMA_AUDITED"
        except (OSError,ValueError,KeyError,zipfile.BadZipFile) as e:r["status"]="HOLD";r["reason"]=str(e)[:240]
        finally:path.unlink(missing_ok=True)
    return {"status":"PUBLIC_SUPPLEMENT_INVENTORY_AND_SCHEMA","record":FIGSHARE,"title":d["title"],"doi":d.get("doi"),"source":source,"files":files,
        "limit":"Workbook schema only, no individual phenotype export, no inferred matching between tables, no outcome-derived modeling. Other file types inventoried, not analyzed."}
def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["GSE277819","GSE42507","ena","supplement"]);p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    a.root.mkdir(parents=True,exist_ok=True)
    try:
        if a.mode in SERIES:r=geo_matrix(a.mode,a.root)
        elif a.mode=="ena":r=ena(a.root)
        else:r=supplement(a.root)
    except (OSError,ValueError,KeyError,TypeError,UnicodeError,EOFError,csv.Error,zipfile.BadZipFile) as e:
        r={"status":"HOLD_SOURCE_OR_SCHEMA","mode":a.mode,"reason":str(e)[:300]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2,allow_nan=False)+"\n");print(json.dumps(r,indent=2,allow_nan=False))
if __name__=="__main__":main()
