#!/usr/bin/env python3
"""Fixed five SNPs, public indexed GWAS only; conditional reported-beta comparisons."""
import argparse, collections, csv, gzip, hashlib, io, json, math, re, struct, time
from pathlib import Path
from urllib.parse import urljoin,urlsplit
from urllib.request import Request,urlopen
import yaml
from qtl_original_indexed_http_range import CheckedReader,tbi_bins_for_region,parse_bgzf_range

PANEL=("rs9515201","rs4771674","rs4918058","rs55983834","rs12415501")
STUDIES={"EUR_AIS":"GCST90104540","EUR_CES":"GCST90104541","EUR_LAS":"GCST90104542","EUR_SVS":"GCST90104543",
         "EAS_AIS":"GCST90104545","EAS_CES":"GCST90104546","EAS_LAS":"GCST90104547","EAS_SVS":"GCST90104548"}
BASE="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/"
MAX_RANGE=12000000
def require(v,msg):
    if not v:raise ValueError(msg)
def sha(x):return hashlib.sha256(x).hexdigest()
def get(url,cap,span=None):
    host=urlsplit(url).hostname
    require(host in {"ftp.ebi.ac.uk","rest.ensembl.org"} and urlsplit(url).scheme=="https","Unapproved URL")
    h={"Accept-Encoding":"identity","User-Agent":"IS-research-indexed-fixed-panel/1.0"}
    if host=="rest.ensembl.org":h["Accept"]="application/json";h["Content-Type"]="application/json"
    if span:h["Range"]=f"bytes={span[0]}-{span[1]}"
    with urlopen(Request(url,headers=h),timeout=40) as f:
        require(urlsplit(f.geturl()).hostname==host,"Unexpected redirect host")
        require(f.status==(206 if span else 200),"HTTP status not expected")
        content_range=f.headers.get("Content-Range","")
        if span:
            m=re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",content_range)
            require(m is not None and (int(m[1]),int(m[2]))==span and int(m[3])>span[1],"Exact range not returned")
        raw=f.read(cap+1)
    require(0<len(raw)<=cap,"Byte bound exceeded or empty")
    if span:require(len(raw)==span[1]-span[0]+1,"Truncated indexed range")
    return raw,content_range
def write(p,x):
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
    print(json.dumps(x,indent=2,allow_nan=False),flush=True)
def reference():
    out={}
    for rs in PANEL:
        try:
            raw,_=get("https://rest.ensembl.org/variation/human/"+rs+"?content-type=application/json",2000000)
            obj=json.loads(raw);require(obj["name"]==rs,"Variant identity changed")
            m=[x for x in obj["mappings"] if x.get("assembly_name")=="GRCh38" and re.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2])",str(x.get("seq_region_name"))) and x.get("strand")==1 and x.get("start")==x.get("end")]
            require(len(m)==1,"Nonunique primary GRCh38 mapping")
            x=m[0];a=x["allele_string"].split("/")
            require(len(a)==2 and all(re.fullmatch("[ACGT]",v) for v in a) and a[0]!=a[1],"Non-biallelic SNP")
            require(set(a) not in (set("AT"),set("CG")),"Palindromic SNP not eligible")
            ch=str(x["seq_region_name"]);pos=int(x["start"])
            if rs in {"rs9515201","rs55983834"}:
                require((ch,pos)==({"rs9515201":("13",110388451),"rs55983834":("10",103841452)}[rs]),"Prior publisher position disagreement")
            out[rs]={"status":"MATCHED_REFERENCE","chr":ch,"pos":pos,"ref":a[0],"alt":a[1],"ensembl_json_sha256":sha(raw)}
        except (OSError,ValueError,KeyError,TypeError) as e:out[rs]={"status":"HOLD_REFERENCE","reason":str(e)[:180]}
        time.sleep(1)
    return out
def parse_index(raw):
    # TBI permits col_end=0 for point records. Preserve exact original hash.
    payload=gzip.decompress(raw);require(len(payload)<80000000,"Index expansion bound")
    rd=CheckedReader(payload);require(rd.take(4).tobytes()==b"TBI\x01","Wrong TBI magic")
    n,fmt,cs,cb,ce,meta,skip,ln=rd.unpack("<8i")
    require(1<=n<=80 and 0<ln<10000 and 1<=cs<=100 and 1<=cb<=100 and 0<=ce<=100,"Invalid TBI header")
    names=rd.take(ln).tobytes();require(names.endswith(b"\0"),"TBI names")
    names=names[:-1].decode("ascii").split("\0");require(len(names)==n and len(set(names))==n,"TBI names")
    refs={}
    for name in names:
        nb,=rd.unpack("<i");require(0<=nb<=500000,"TBI bins")
        bins={}
        for _ in range(nb):
            b,=rd.unpack("<I");nc,=rd.unpack("<i");require(0<=nc<=1000000,"TBI chunks")
            chunks=[rd.unpack("<QQ") for _ in range(nc)]
            if b==37450:require(nc==2,"TBI pseudo-bin");continue
            require(b not in bins and all(y>=x for x,y in chunks),"TBI chunk bounds");bins[b]=chunks
        ni,=rd.unpack("<i");require(0<=ni<=300000,"TBI linear index")
        refs[name]=(bins,[rd.unpack("<Q")[0] for _ in range(ni)])
    require(len(payload)-rd.offset in (0,8),"TBI trailing bytes")
    return {"refs":refs,"seq":cs,"begin":cb,"end":ce,"format":fmt,"meta":meta}
def indexed_span(tbi,ch,start,end,total):
    names=[v for v in (ch,"chr"+ch) if v in tbi["refs"]]
    require(len(names)==1,"Missing/ambiguous index chromosome")
    bins,lin=tbi["refs"][names[0]]
    off=lin[(start-1)>>14] if (start-1)>>14<len(lin) else 0
    arr=[(x,y) for k in tbi_bins_for_region(start,end) for x,y in bins.get(k,[]) if y>off]
    if not arr:return None
    x=max(min(x for x,y in arr),off);y=max(y for x,y in arr)
    a=x>>16;b=min(total-1,(y>>16)+65535)
    require(0<=a<=b and b-a+1<=MAX_RANGE,"Indexed range exceeds budget")
    return a,b,x&65535,(y>>16)-a
def choose(h,*options):
    for x in options:
        if x in h:return x
    raise ValueError("Missing field: "+"/".join(options))
def source_fields(h,tbi):
    fields={"chr":choose(h,"hm_chrom","chromosome"),"pos":choose(h,"hm_pos","base_pair_location"),
        "ea":choose(h,"hm_effect_allele","effect_allele"),"oa":choose(h,"hm_other_allele","other_allele"),
        "beta":choose(h,"hm_beta","beta"),"se":choose(h,"standard_error","hm_standard_error"),
        "p":choose(h,"p_value","hm_p_value"),"eaf":choose(h,"hm_effect_allele_frequency","effect_allele_frequency"),
        "code":choose(h,"hm_code","harmonisation_code")}
    require(h[tbi["seq"]-1]==fields["chr"] and h[tbi["begin"]-1]==fields["pos"] and tbi["end"] in (0,tbi["begin"]),"Index/header coordinate contract")
    fields["rsid"]=next((x for x in ("hm_rsid","rsid","rs_id") if x in h),None)
    return fields
def orient(row,f,r,rs):
    if f["rsid"] and row[f["rsid"]].strip().lower()!=rs:raise ValueError("rsID differs at exact position")
    ea,oa=row[f["ea"]].upper(),row[f["oa"]].upper()
    require({ea,oa}=={r["ref"],r["alt"]},"Alleles do not match reference")
    require(int(row[f["code"]]) in {10,11,12,13},"Not a resolved nonpalindromic harmonisation code")
    b,se,p,eaf=[float(row[f[k]]) for k in ("beta","se","p","eaf")]
    require(all(math.isfinite(v) for v in (b,se,p,eaf)) and se>0 and 0<=p<=1 and 0<eaf<1,"Invalid numeric fields")
    sign=1 if ea==r["alt"] else -1
    return {"status":"MATCHED_HARMONISED_ROW","beta_alt":sign*b,"se":se,"p":p,"alt_eaf":eaf if sign==1 else 1-eaf,"hm_code":int(row[f["code"]]),"ci95":[sign*b-1.959964*se,sign*b+1.959964*se]}
def extract(study,ref):
    acc=STUDIES[study];out={"study":study,"accession":acc,"status":"HOLD_SOURCE","rows":{rs:{"status":"HOLD_SOURCE"} for rs in PANEL}}
    try:
        directory=BASE+acc+"/harmonised/"
        html,_=get(directory,2000000)
        names=sorted(set(re.findall(r'href=["\']([^"\']+\.h\.tsv\.gz)["\']',html.decode("utf-8"))))
        require(len(names)==1 and re.fullmatch(re.escape(acc)+r"[^/]*\.h\.tsv\.gz",names[0]),"One exact harmonised file not established")
        url=urljoin(directory,names[0]);require(url.startswith(directory),"Directory escape")
        meta,_=get(url+"-meta.yaml",200000)
        md=yaml.safe_load(meta);require(isinstance(md,dict),"Metadata schema")
        require(md.get("gwas_id")==acc and md.get("data_file_name")==names[0],"Harmonised metadata accession/filename")
        require(str(md.get("genome_assembly","")).upper()=="GRCH38" and md.get("is_harmonised") is True,"GRCh38 harmonisation not established")
        out.update(source=url,metadata_sha256=sha(meta),publisher_full_file_md5=md.get("data_file_md5sum"),metadata={k:md.get(k) for k in ("genome_assembly","trait_description","samples","sample_ancestry","is_harmonised")})
        time.sleep(3)
        idx,_=get(url+".tbi",8000000);tbi=parse_index(idx);out["index_sha256"]=sha(idx)
        time.sleep(3)
        first,cr=get(url,65536,(0,65535));total=int(cr.rsplit("/",1)[1]);out["header_range_sha256"]=sha(first);out["source_total_bytes"]=total
        hdr=parse_bgzf_range(first,0,0)["data_lines"].readline().lstrip("#").rstrip("\r\n").split("\t")
        require(len(set(hdr))==len(hdr),"Duplicate header");f=source_fields(hdr,tbi);out["fields"]=f
        by=collections.defaultdict(list)
        for rs,r in ref.items():
            if r["status"]=="MATCHED_REFERENCE":by[r["chr"]].append((rs,r))
            else:out["rows"][rs]=r
        out["ranges"]=[]
        for ch,rr in sorted(by.items()):
            for rs,r in rr:
                try:
                    span=indexed_span(tbi,ch,r["pos"],r["pos"],total)
                    if span is None:out["rows"][rs]={"status":"HOLD_NO_INDEX_CHUNK"};continue
                    a,b,x,y=span;time.sleep(3)
                    raw,cr=get(url,b-a+1,(a,b));stream=parse_bgzf_range(raw,x,y)["data_lines"]
                    out["ranges"].append({"rsid":rs,"range":cr,"sha256":sha(raw),"bytes":len(raw)})
                    matches=[]
                    for row in csv.reader(stream,delimiter="\t"):
                        if len(row)!=len(hdr):continue
                        d=dict(zip(hdr,row))
                        if d[f["chr"]].removeprefix("chr")==ch and d[f["pos"]]==str(r["pos"]):
                            if {d[f["ea"]].upper(),d[f["oa"]].upper()}=={r["ref"],r["alt"]}:matches.append(d)
                    require(len(matches)==1,"No unique allele-compatible exact coordinate row")
                    out["rows"][rs]=orient(matches[0],f,r,rs)
                except (OSError,ValueError,KeyError,TypeError,IndexError) as e:out["rows"][rs]={"status":"HOLD_VARIANT","reason":str(e)[:180]}
        out["status"]="INDEXED_ROWS_READBACK"
        out["matched"]=sum(r["status"]=="MATCHED_HARMONISED_ROW" for r in out["rows"].values())
    except (OSError,ValueError,KeyError,TypeError,IndexError) as e:out["reason"]=str(e)[:240]
    return out
def contrast(a,b):
    delta=a["beta_alt"]-b["beta_alt"];s0=math.hypot(a["se"],b["se"]);sw=a["se"]+b["se"]
    p0=math.erfc(abs(delta)/(s0*math.sqrt(2)));pw=math.erfc(abs(delta)/(sw*math.sqrt(2)))
    return {"status":"CONDITIONAL_REPORTED_BETA_CONTRAST","delta":delta,"ci95_max_variance":[delta-1.959964*sw,delta+1.959964*sw],"p_zero_covariance":p0,"p_max_variance":pw,"bonferroni_40_zero":min(1,40*p0),"bonferroni_40_max":min(1,40*pw)}
def combine(root):
    d={}
    for k,acc in STUDIES.items():
        p=root/k/(k+".json")
        if p.exists():
            x=json.loads(p.read_text());require(x["study"]==k and x["accession"]==acc and set(x["rows"])==set(PANEL),"Study artifact identity");d[k]=x
        else:d[k]={"study":k,"accession":acc,"status":"HOLD_NO_ARTIFACT","rows":{rs:{"status":"HOLD_NO_ARTIFACT"} for rs in PANEL}}
    pairs=[("EAS_"+t,"EUR_"+t) for t in ("AIS","CES","LAS","SVS")]+[(a+"_SVS",a+"_"+t) for a in ("EUR","EAS") for t in ("CES","LAS")]
    res=[]
    for rs in PANEL:
        for a,b in pairs:
            x=d[a]["rows"][rs];y=d[b]["rows"][rs]
            row={"rsid":rs,"first":a,"second":b,"status":"HOLD_UNMATCHED"}
            if x["status"]==y["status"]=="MATCHED_HARMONISED_ROW":row.update(contrast(x,y))
            res.append(row)
    require(len(res)==40,"Fixed contrast family")
    return {"status":"EXPLORATORY_INDEXED_COMPARISONS","family":40,"matched":sum(r["status"]=="CONDITIONAL_REPORTED_BETA_CONTRAST" for r in res),"studies":d,"contrasts":res,
        "limits":["Indexed harmonised rows, not full original-file checksum validation. Exact source phenotype metadata must be reviewed against publisher study mapping.","Reported beta comparisons conditional on compatible effect models/scales; covariance bounds address unknown estimator covariance only. Not case-case regression or evidence of ancestry-specific mechanism.","Five preselected known variants, no independent replication, no new MR or colocalization. Missing is not zero."]}
def main():
    ap=argparse.ArgumentParser();ap.add_argument("mode",choices=["reference","extract","combine","selftest"]);ap.add_argument("--study",choices=STUDIES);ap.add_argument("--ref",type=Path);ap.add_argument("--root",type=Path);ap.add_argument("--out",type=Path);a=ap.parse_args()
    if a.mode=="selftest":
        f={k:k for k in ("chr","pos","ea","oa","beta","se","p","eaf","code","rsid")}
        r={"ref":"C","alt":"A"};x={"ea":"C","oa":"A","beta":"0.2","se":"0.1","p":"0.05","eaf":"0.3","code":"10","rsid":"rs9515201"}
        y=orient(x,f,r,"rs9515201");assert y["beta_alt"]==-.2 and y["alt_eaf"]==.7
        assert contrast(y,y)["p_max_variance"]==1
        assert len(STUDIES)==8 and len(PANEL)==5
        print("Focused tests passed; no biological result.");return
    result=reference() if a.mode=="reference" else extract(a.study,json.loads(a.ref.read_text())) if a.mode=="extract" else combine(a.root)
    write(a.out,result)
if __name__=="__main__":main()
