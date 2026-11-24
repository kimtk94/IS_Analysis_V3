#!/usr/bin/env python3
"""PRJNA781014 exact public-source audit; score only explicit complete pairs."""
import argparse,collections,csv,hashlib,io,json,math,re,statistics,urllib.request,zipfile
from pathlib import Path
from ocular_frozen_241_external import PIN,read_panel
from ocular_paired_panel_validation import gene_counts,summarize_scores
PROJECT="PRJNA781014"
FIELDS=("run_accession","secondary_study_accession","study_accession","experiment_accession","sample_accession","sample_alias","library_name","experiment_title")
TREATMENTS={"DEX":"DEX","DEXAMETHASONE":"DEX","ETH":"VEH","ETHANOL":"VEH","VEH":"VEH","VEHICLE":"VEH","CONTROL":"VEH","CTRL":"VEH"}
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
def parsed(label):
    # Only explicit treatment tokens; a numeric common stem is required.
    s=label.strip().upper()
    hits=list(re.finditer(r"(?<![A-Z])(?:DEXAMETHASONE|DEX|ETHANOL|ETH|VEHICLE|VEH|CONTROL|CTRL)(?![A-Z])",s))
    if len(hits)!=1:return None
    m=hits[0];stem=re.sub(r"[^A-Z0-9]","",s[:m.start()]+s[m.end():])
    if not stem or not re.search(r"\d",stem):return None
    return stem,TREATMENTS[m.group()]
def candidates(rows,out):
    options=[];observed={}
    for field in ("sample_alias","library_name","experiment_title"):
        groups=collections.defaultdict(dict);bad=dup=0
        per_run={}
        for r in rows:
            v=parsed(r[field])
            if v is None:bad+=1;continue
            stem,arm=v
            if arm in groups[stem]:dup+=1
            groups[stem][arm]=r["run_accession"]
            per_run[r["run_accession"]]=arm
        observed[field]=per_run
        pairs=sorted((v["DEX"],v["VEH"]) for v in groups.values() if set(v)=={"DEX","VEH"})
        complete=bad==dup==0 and len(groups)==5 and len(pairs)==5
        out[field]={"unparsed":bad,"duplicate_arm_in_group":dup,"groups":len(groups),
                    "paired_groups":len(pairs),"complete_5_pair_label_contract":complete,
                    "masked_patterns":dict(collections.Counter(re.sub(r"\d+","#",r[field])[:180] for r in rows))}
        if complete:options.append((field,pairs))
    if not options:return None
    require(all(p==options[0][1] for _,p in options),"Conflicting complete pairing label sources")
    pairs=options[0][1]
    truearm={r:arm for a,b in pairs for r,arm in ((a,"DEX"),(b,"VEH"))}
    require(all(truearm[r]==a for d in observed.values() for r,a in d.items()),"Treatment label contradiction")
    return {"fields":[f for f,_ in options],"pairs":pairs}
def one(z,name):
    hits=[x for x in z.namelist() if x.rsplit("/",1)[-1]==name]
    require(len(hits)==1,"DEE2 member "+name)
    return hits[0]
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    root=a.root;out={"status":"HOLD_SOURCE","project":PROJECT}
    try:
        path=root/"ena.tsv"
        url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession="+PROJECT+"&result=read_run&fields="+",".join(FIELDS)+"&format=tsv"
        out["ena"]=download(url,path,200000)
        with path.open(encoding="utf-8",newline="") as f:
            r=csv.DictReader(f,delimiter="\t")
            require(set(FIELDS)<=set(r.fieldnames or []),"ENA field contract")
            rows=sorted(r,key=lambda x:x["run_accession"])
        sizes={k:len({x[k] for x in rows}) for k in ("run_accession","experiment_accession","sample_accession")}
        out["source_inventory"]={"rows":len(rows),**sizes}
        require(len(rows)==10 and all(n==10 for n in sizes.values()),"Ten distinct source sample/experiment/run contract")
        require(all(x["study_accession"]==PROJECT for x in rows),"Project mismatch")
        ids={x["secondary_study_accession"] for x in rows}
        require(len(ids)==1 and re.fullmatch(r"SRP\d+",next(iter(ids))),"Secondary study ID")
        srp=next(iter(ids));out["secondary_study"]=srp
        projection="".join("\t".join(x[k] for k in FIELDS)+"\n" for x in rows)
        out["ena_sorted_projection_sha256"]=hashlib.sha256(projection.encode()).hexdigest()
        out["pairing_label_audit"]={}
        design=candidates(rows,out["pairing_label_audit"])
        out["pairing"]={"status":"EXPLICIT_LABEL_PAIR_CANDIDATE","fields":design["fields"],"pairs":5} if design else {"status":"HOLD_PAIRING"}
        out["status"]="16H_SOURCE_INVENTORY_COMPLETE"
        bundle=root/"counts.zip"
        out["dee2"]=download("https://dee2.io/huge/hsapiens/"+srp+"_NA.zip",bundle,50000000)
        with zipfile.ZipFile(bundle) as z:
            require(z.testzip() is None,"DEE2 CRC")
            with z.open(one(z,"MetadataSummary.tsv")) as f:
                rr=list(csv.DictReader(io.TextIOWrapper(f,encoding="utf-8"),delimiter="\t"))
            by={x["run_accession"]:x for x in rows}
            runs={x["SRR_accession"]:x for x in rr}
            require(len(rr)==len(runs)>0,"Duplicate or empty DEE2 run list")
            require(all(k in by and x["SRX_accession"]==by[k]["experiment_accession"] and x["SRP_accession"]==srp and x["QC_summary"] in ("PASS","WARN","FAIL") for k,x in runs.items()),"DEE2 source/run/QC mismatch")
            with z.open(one(z,"GeneCountMatrix.tsv")) as f:
                header=f.readline().decode("utf-8").rstrip("\r\n").split("\t")
            require(header[0]=="GeneID" and len(header[1:])==len(set(header[1:])) and set(header[1:])==set(runs),"Count-metadata column mismatch")
            out["count_release"]={"counted_runs":len(rr),"expected_runs":10,"qc_counts":dict(collections.Counter(x["QC_summary"] for x in rr)),"missing_run_count":10-len(rr)}
            out["status"]="16H_COUNTS_QUALIFIED"
            if design is None:
                out["score_gate"]="HOLD_PAIRING"
            else:
                panel=read_panel(root/"panel.zip")
                with z.open(one(z,"GeneCountMatrix.tsv")) as f:
                    expr,n_genes=gene_counts(io.TextIOWrapper(f,encoding="utf-8"),panel,runs)
                groups={"all_241":sorted(expr),"ECM_16":sorted(g for g in expr if panel[g]["ecm_member"])}
                require(len(groups["all_241"])>=math.ceil(.9*241) and len(groups["ECM_16"])==16,"Fixed panel mapping gate")
                scores={k:[] for k in groups};primary={k:[] for k in groups}
                n_complete=n_primary=0
                for dex,veh in design["pairs"]:
                    if dex not in runs or veh not in runs:continue
                    n_complete+=1
                    good=all(runs[x]["QC_summary"] in ("PASS","WARN") for x in (dex,veh))
                    n_primary+=int(good)
                    for name,genes in groups.items():
                        val=statistics.mean((expr[g][dex]-expr[g][veh])*(1 if panel[g]["acute_log2fc"]>0 else -1) for g in genes)
                        scores[name].append(val)
                        if good:primary[name].append(val)
                out["score_gate"]="DESCRIPTIVE_EXPLICIT_LABEL_PAIRED_SCORE" if n_primary else "HOLD_NO_QC_COMPLETE_PAIRS"
                out["panel_sha256"]=PIN
                out["mapped_genes"]={k:len(v) for k,v in groups.items()}
                out["count_genes"]=n_genes
                out["complete_public_pairs"]=n_complete
                out["primary_qc_pairs"]=n_primary
                out["primary"]=summarize_scores(primary)
                out["all_complete_sensitivity"]=summarize_scores(scores)
                out["status"]="16H_PAIRED_SCORE_COMPUTED" if n_primary else "16H_COUNTS_QUALIFIED"
    except (OSError,ValueError,KeyError,zipfile.BadZipFile,UnicodeError,csv.Error) as e:
        out["hold_reason"]=str(e)[:240]
    out["limits"]=["Original paper reports n=6 at 16h; ENA returns ten distinct public runs/samples. Unreleased/missing source samples and donor identity not independently resolved.",
        "Pair candidates use explicit treatment labels and matching numeric stems only; no matching by run order. No new clinical endpoints.",
        "No donor IDs, donor-level values, or expression matrix exported. Score based on exact frozen signs and assigned-gene library totals.",
        "Canonical sorted projection recorded rather than pinning order-sensitive raw bytes, because the two earlier 1,647-byte responses had different hashes.",
        "Two-panel Holm controls these two descriptive sign tests, not prior selection or cohort overlap. Five complete pairs cannot give two-sided P<.05.",
        "Different duration, culture and sequencing sites preclude causal time comparison. No negative inference from absent source records."]
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps(out,indent=2,allow_nan=False))
if __name__=="__main__":main()
