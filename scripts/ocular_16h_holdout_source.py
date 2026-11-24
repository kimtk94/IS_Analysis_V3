#!/usr/bin/env python3
"""Bounded public 16h bulk TM holdout source qualification, before scoring expression."""
import argparse,collections,csv,hashlib,io,json,re,urllib.request,zipfile
from pathlib import Path
PROJECT="PRJNA781014"
FIELDS=("run_accession","secondary_study_accession","study_accession","experiment_accession","sample_accession","sample_alias","library_name","experiment_title")
def require(ok,msg):
    if not ok:raise ValueError(msg)
def download(url,path,cap):
    req=urllib.request.Request(url,headers={"User-Agent":"Public-TM-research-source-audit/1.0"})
    with urllib.request.urlopen(req,timeout=40) as f,path.open("wb") as w:
        h=hashlib.sha256();n=0
        while b:=f.read(1024*1024):
            n+=len(b);require(n<=cap,"Bounded source cap")
            h.update(b);w.write(b)
    require(n>0,"Empty public response")
    return {"url":url,"bytes":n,"sha256":h.hexdigest()}
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession="+PROJECT+"&result=read_run&fields="+",".join(FIELDS)+"&format=tsv"
    path=a.root/"ena.tsv";out={"status":"HOLD_SOURCE","project":PROJECT}
    try:
        out["ena"]=download(url,path,200000)
        with path.open(encoding="utf-8",newline="") as f:
            r=csv.DictReader(f,delimiter="\t")
            require(set(FIELDS)<=set(r.fieldnames or []),"ENA exact field schema")
            rows=sorted(r,key=lambda x:x["run_accession"])
        require(len(rows)==12 and len({x["run_accession"] for x in rows})==12 and len({x["experiment_accession"] for x in rows})==12,"Expected twelve distinct study runs/experiments")
        require(all(x["study_accession"]==PROJECT for x in rows),"Unexpected study accession")
        ids={x["secondary_study_accession"] for x in rows}
        require(len(ids)==1 and re.fullmatch(r"SRP\d+",next(iter(ids))),"Secondary study ID")
        srp=next(iter(ids));out["secondary_study"]=srp
        projection="".join("\t".join(x[k] for k in FIELDS)+"\n" for x in rows)
        out["ena_projection_sha256"]=hashlib.sha256(projection.encode()).hexdigest()
        # Public run/experiment identifiers may be retained, but not sample IDs
        # or donor labels. Record normalized label patterns and a strict pairing
        # candidate only when an arm suffix is explicit and corroborated.
        patterns={}
        for field in ("sample_alias","library_name","experiment_title"):
            counts=collections.Counter(re.sub(r"\d+","#",x[field].strip()) for x in rows)
            patterns[field]=dict(counts)
        out["label_patterns"]=patterns
        # Explicit suffix only; no matching by run order. Further source review
        # required before interpreting candidate pairs or scoring expression.
        candidate={}
        for field in ("sample_alias","library_name"):
            groups=collections.defaultdict(list);bad=0
            for x in rows:
                m=re.fullmatch(r"(.+?)[ _-]?(DEX|ETH|VEH|VEHICLE|CONTROL|CTRL)",x[field].strip(),re.I)
                if not m:bad+=1;continue
                arm="DEX" if m.group(2).upper()=="DEX" else "VEH"
                groups[m.group(1).rstrip(" _-").casefold()].append(arm)
            candidate[field]={"unresolved_rows":bad,"candidate_pair_count":sum(sorted(v)==["DEX","VEH"] for v in groups.values()),"groups":len(groups)}
        out["pairing_candidates_not_yet_verified"]=candidate
        bundle=a.root/"counts.zip"
        try:
            out["dee2"]=download("https://dee2.io/huge/hsapiens/"+srp+"_NA.zip",bundle,50000000)
            with zipfile.ZipFile(bundle) as z:
                require(z.testzip() is None,"DEE2 CRC")
                def one(name):
                    members=[x for x in z.namelist() if x.rsplit("/",1)[-1]==name]
                    require(len(members)==1,"DEE2 member "+name)
                    return members[0]
                with z.open(one("MetadataSummary.tsv")) as f:
                    rr=list(csv.DictReader(io.TextIOWrapper(f,encoding="utf-8"),delimiter="\t"))
                by={x["run_accession"]:x for x in rows}
                require(len(rr)==len({x["SRR_accession"] for x in rr}),"Duplicate DEE2 run")
                require(all(x["SRR_accession"] in by and x["SRX_accession"]==by[x["SRR_accession"]]["experiment_accession"] and x["SRP_accession"]==srp and x["QC_summary"] in ("PASS","WARN","FAIL") for x in rr),"DEE2 source/run/QC mismatch")
                with z.open(one("GeneCountMatrix.tsv")) as f:
                    header=f.readline().decode("utf-8").rstrip("\r\n").split("\t")
                require(header[0]=="GeneID" and len(header[1:])==len(set(header[1:])) and set(header[1:])=={x["SRR_accession"] for x in rr},"Count-metadata column mismatch")
                out["count_release"]={"counted_runs":len(rr),"expected_runs":len(rows),"qc_counts":dict(collections.Counter(x["QC_summary"] for x in rr)),"missing_run_count":12-len(rr)}
                out["status"]="16H_PUBLIC_SOURCE_QUALIFIED_PENDING_EXACT_PAIR_REVIEW"
        except (OSError,ValueError,KeyError,zipfile.BadZipFile,UnicodeError) as e:
            out["count_release"]={"status":"HOLD","reason":str(e)[:240]}
            out["status"]="16H_ENA_ONLY"
    except (OSError,ValueError,KeyError,UnicodeError,csv.Error) as e:
        out["reason"]=str(e)[:240]
    out["limits"]=["No expression scores calculated; no donor-level identifiers or expression exported.",
        "16h and 7d cohorts in Kathirvel 2022 used different eye banks and sequencing centers, not randomized time contrasts.",
        "Candidate suffix groups are not approved sample pairing. Review original design and metadata before scoring.",
        "Existing 241 panel was selected on GSE313736 and 100nM 48h external statistics; this study is a distinct candidate validation source."]
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps(out,indent=2,allow_nan=False))
if __name__=="__main__":main()
