#!/usr/bin/env python3
"""Recover rs2767060 single-index-variant LD API schema and audit exact proxy tokens.
Post-readout reference-only lookup, no statistical association or novelty claim.
"""
import argparse,collections,hashlib,importlib,json,math,re
from pathlib import Path
from openpyxl import load_workbook
m=importlib.import_module("scripts.is.rs2767060_regional_annotation")
PIN="7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4"
RS="rs2767060"
def require(ok,msg):
    if not ok:raise ValueError(msg)
def read_ld(v,pop):
    require(isinstance(v,list),"Expected LD array")
    if not v:return {"status":"HOLD_EMPTY_LD_RESPONSE","returned_rows":0,"proxies":[]}
    proxies={};schemas=collections.Counter()
    for row in v:
        require(isinstance(row,dict),"LD row type")
        schemas[",".join(sorted(row))]+=1
        # The single-ID endpoint can report the partner as variation;
        # pairwise-shaped payloads must explicitly include the query variant.
        if {"variation1","variation2"}<=set(row):
            a,b=row["variation1"],row["variation2"]
            require(RS in (a,b),"Query variant absent from pair")
            other=b if a==RS else a
        elif "variation" in row and isinstance(row["variation"],str):
            other=row["variation"]
        else:raise ValueError("Unknown LD partner schema: "+",".join(sorted(row)))
        require(isinstance(other,str) and re.fullmatch(r"rs[0-9]+",other),"Proxy rsID")
        require(row.get("population_name")==pop,"Population identity")
        r2,d=float(row["r2"]),float(row["d_prime"])
        require(math.isfinite(r2) and math.isfinite(d) and .8-1e-8<=r2<=1 and abs(d)<=1,"LD statistic")
        if other==RS:continue
        item={"rsid":other,"r2":r2,"d_prime":d}
        if "chr" in row or "start" in row or "end" in row:
            # API JSON may encode integer coordinate fields as decimal strings.
            # Validate syntax before conversion; never round or infer coordinates.
            start,end=row.get("start"),row.get("end")
            require(re.fullmatch(r"[0-9]+",str(start)) is not None and re.fullmatch(r"[0-9]+",str(end)) is not None,"Proxy coordinate format: "+str(type(start).__name__)+"/"+str(type(end).__name__))
            start,end=int(start),int(end)
            # Ensembl insertion convention: start = end + 1, not an invalid interval.
            # https://www.ensembl.org/info/genome/variation/prediction/classification.html
            require(str(row.get("chr")).removeprefix("chr")=="10" and end>0 and 0<start<=end+1 and max(abs(start-35833363),abs(end-35833363))<=500000,"Proxy coordinate identity: "+str([row.get("chr"),start,end]))
            item.update(chr="10",start=start,end=end,coordinate_class=("INSERTION" if start==end+1 else "POINT" if start==end else "SPAN"),consequence_type=row.get("consequence_type"))
        require(other not in proxies or proxies[other]==item,"Duplicate LD conflict")
        proxies[other]=item
    return {"status":"REFERENCE_LD_PARSED" if proxies else "HOLD_NO_PARSED_PROXIES","returned_rows":len(v),"schemas":dict(schemas),"proxies":sorted(proxies.values(),key=lambda x:(-x["r2"],x["rsid"]))}
def audit(evidence,workbook):
    p=json.loads(evidence.read_text());x=p["result"]
    require(p["provenance"]["run_id"]==45539259361 and x["identity_verified"] and x["rsid"]==RS and x["ensembl_allele_string"]=="C/T" and x["GRCh38_pos"]==35833363,"Prior allele identity")
    out={"status":"POST_READOUT_LD_PROXY_PUBLICATION_AUDIT","rsid":RS,"populations":{},"publication_matches":[],"limits":["Single-ID Ensembl 1000 Genomes Phase 3 LD context, not GWAS-matched individual-level LD. Exact variant C/T identity verified in prior reference run. Empty or unparsed responses remain HOLD.","Exact rsID tokens in the pinned one-workbook audit do not prove phenotype association or novelty. LD proxy evidence does not identify a causal gene or causal variant."]}
    for pop in m.POPULATIONS:
        v,source=m.request("/ld/human/"+RS+"/"+pop,{"window_size":500,"r2":.8,"attribs":1})
        expected_sha={"1000GENOMES:phase_3:CHB":"c5a134d1e802955b41a8f28dfa6e262758e11ab2ba096a53f1ef3ae09c3e3835","1000GENOMES:phase_3:JPT":"8d711ef0f479b4c46545bb10092fa9cb9f25972fbcb6ef33d501e4d99787c528","1000GENOMES:phase_3:CEU":"c642d54161a64e4146ccb43276c1fdfcb7eed3b0ea63ae5a001226d0b18a0007"}
        try:
            require(source.get("sha256")==expected_sha[pop],"LD response changed from initial verified raw response digest")
            d=read_ld(v,pop)
        except (ValueError,KeyError,TypeError) as e:d={"status":"HOLD_SCHEMA","returned_rows":len(v) if isinstance(v,list) else None,"observed_keys":sorted(v[0]) if isinstance(v,list) and v and isinstance(v[0],dict) else [],"reason":str(e)[:180],"proxies":[]}
        d["source"]=source;out["populations"][pop]=d
    sets={k:{v["rsid"] for v in d["proxies"]} for k,d in out["populations"].items() if d["status"]=="REFERENCE_LD_PARSED"}
    out["reference_populations_parsed"]=len(sets)
    out["proxy_union"]=len(set().union(*sets.values())) if sets else 0
    out["proxy_shared_all_three"]=len(set.intersection(*sets.values())) if len(sets)==3 else None
    require(workbook.stat().st_size==756910 and hashlib.sha256(workbook.read_bytes()).hexdigest()==PIN,"Workbook pin")
    if sets:
        allrs=set().union(*sets.values())
        w=load_workbook(workbook,read_only=True,data_only=True,keep_links=False)
        try:
            for s in w.worksheets:
                require((s.max_row or 0)<=100000 and (s.max_column or 0)<=400,"Workbook bounds")
                for i,row in enumerate(s.iter_rows(values_only=True),1):
                    vals=["" if v is None else str(v) for v in row]
                    matches=set()
                    for v in vals:matches.update(rs.lower() for rs in re.findall(r"(?<![A-Za-z0-9])rs[0-9]+(?![A-Za-z0-9])",v,re.I) if rs.lower() in allrs)
                    for rs in sorted(matches):
                        out["publication_matches"].append({"rsid":rs,"sheet":s.title,"row":i,"populations":[k for k,ss in sets.items() if rs in ss],"values":vals[:40]})
        finally:w.close()
    out["workbook_sha256"]=PIN
    return out
def main():
    # Focused checks of documented insertion coordinates and strict invalid interval rejection.
    pop="1000GENOMES:phase_3:CHB";base={"variation":"rs123","population_name":pop,"r2":"0.9","d_prime":"1","chr":"10","start":"35832725","end":"35832724"}
    assert read_ld([base],pop)["proxies"][0]["coordinate_class"]=="INSERTION"
    try:read_ld([{**base,"end":"35832723"}],pop)
    except ValueError:pass
    else:raise AssertionError("Invalid reversed interval accepted")
    p=argparse.ArgumentParser();p.add_argument("--evidence",type=Path,required=True);p.add_argument("--workbook",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    x=audit(a.evidence,a.workbook);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
    short={**x,"populations":{k:{**{a:b for a,b in v.items() if a!="proxies"},"proxy_count":len(v["proxies"]),"coordinate_classes":dict(collections.Counter(z.get("coordinate_class","UNSPECIFIED") for z in v["proxies"])),"first_5":v["proxies"][:5]} for k,v in x["populations"].items()}}
    print(json.dumps(short,indent=2,allow_nan=False))
if __name__=="__main__":main()
