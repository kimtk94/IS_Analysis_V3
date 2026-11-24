#!/usr/bin/env python3
"""Public GSE313736 metadata-only design gate. No genotype/clinical donor fields retrieved or emitted."""
import argparse,collections,csv,gzip,hashlib,json,re,urllib.request
from pathlib import Path

BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313736"
FILES=("GSE313736_HTMC_scRNA_cell_metadata.tsv.gz","GSE313736_HTMC_scRNA_demultiplexing_and_cell_classification.txt.gz")
def norm(s):return re.sub(r"[^a-z0-9]","",s.lower())
def fetch(url,p,cap):
    req=urllib.request.Request(url,headers={"User-Agent":"Public-academic-expression-research/1.0"})
    with urllib.request.urlopen(req,timeout=40) as f,p.open("wb") as w:
        h=hashlib.sha256();n=0
        while b:=f.read(1024*1024):
            n+=len(b)
            if n>cap:raise ValueError("Bounded transfer exceeded")
            h.update(b);w.write(b)
    if n==0:raise ValueError("Empty source")
    return {"url":url,"bytes":n,"sha256":h.hexdigest()}
def inspect_table(p):
    with gzip.open(p,"rt",encoding="utf-8-sig",newline="") as f:
        head=f.readline().rstrip("\r\n")
        if "\t" not in head:raise ValueError("Unexpected source delimiter")
        h=head.split("\t");keys=[norm(x) for x in h]
        if len(set(h))!=len(h):raise ValueError("Duplicate source columns")
        # Output schema and aggregate counts only. Never export cell barcodes,
        # donor IDs, age/sex/ancestry/disease/PRS fields or genotype data.
        donor_names={"donor","donorid","donoridentity","donorname","individual","individualid","subject","subjectid","sampleid","best","demuxletbest","demuxletdonor","demuxletassignment"}
        arm_names={"condition","treatment","treated","dex","dexstatus","group","treatmentgroup","drug","stim","stimulation","experimentalcondition","library","libraryid","libraryname","origident","pool","sample","samplename"}
        donor_cols=[i for i,x in enumerate(keys) if x in donor_names]
        arm_cols=[i for i,x in enumerate(keys) if x in arm_names]
        # Strict source token checks, not substring extraction from donor genotype or other fields.
        donor_patterns=collections.Counter();arm_counts=collections.defaultdict(collections.Counter)
        joints=collections.defaultdict(collections.Counter)
        # Include no field values in output except normalized treatment counts.
        rd=csv.reader(f,delimiter="\t");n=bad=0
        def donor(v):
            v=v.strip()
            m=re.fullmatch(r"(?:SNG-)?HTMC[-_]?(\d+)",v,re.I)
            return int(m.group(1)) if m else None
        def arm(v):
            v=v.strip().casefold().replace(" ","").replace("_","-")
            if v in ("dex","dexamethasone","treated","dex-treated","500nm","500nmdex","dex500nm"):return "DEX"
            if v in ("vehicle","control","untreated","ctrl","dmso","eth","no-dex"):return "VEH"
            if re.fullmatch(r"pool-?dex[a-z]",v):return "DEX"
            if re.fullmatch(r"pool-?[a-z]",v):return "VEH"
            return None
        for row in rd:
            if not row:continue
            n+=1
            if n>2000000:raise ValueError("Metadata row cap")
            if len(row)!=len(h):bad+=1;continue
            ds={i:donor(row[i]) for i in donor_cols};aa={i:arm(row[i]) for i in arm_cols}
            for i,d in ds.items():
                if d is not None:donor_patterns[i]+=1
            for i,a in aa.items():
                arm_counts[i][a or "UNRESOLVED"]+=1
            for i,d in ds.items():
                if d is None:continue
                for j,a in aa.items():
                    if a is not None:joints[i,j][d,a]+=1
        if n<1000 or bad:raise ValueError("Metadata incomplete/schema failure "+str((n,bad)))
        # No field choice based on expression; report every exact candidate pair.
        summaries=[]
        for (i,j),q in sorted(joints.items()):
            dex={d for d,a in q if a=="DEX"};veh={d for d,a in q if a=="VEH"}
            paired=dex&veh
            summaries.append({"donor_column":h[i],"arm_column":h[j],"donors_dex":len(dex),"donors_vehicle":len(veh),
                  "donors_with_both":len(paired),"recognized_donor_arm_rows":sum(q.values()),
                  "paired_with_at_least_50_rows_per_arm":sum(min(q[d,"DEX"],q[d,"VEH"])>=50 for d in paired),
                  "paired_with_at_least_100_rows_per_arm":sum(min(q[d,"DEX"],q[d,"VEH"])>=100 for d in paired)})
        return {"header":h,"rows":n,"malformed":bad,"donor_id_pattern_counts":{h[i]:donor_patterns[i] for i in donor_cols},
                "treatment_label_counts":{h[i]:dict(arm_counts[i]) for i in arm_cols},"candidate_pair_coverage":summaries,
                "status":"AGGREGATE_METADATA_DESIGN_ONLY",
                "limit":"Row counts not validated unique cells or QC singlets. Column/label interpretation must be confirmed against original study; no expression or new P value."}
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    results=[]
    for name in FILES:
        dest=a.root/name;url=BASE+"/suppl/"+name
        r={"file":name}
        try:r.update(fetch(url,dest,80000000));r["analysis"]=inspect_table(dest)
        except (OSError,ValueError,EOFError,UnicodeError,csv.Error) as e:r["analysis"]={"status":"HOLD","reason":str(e)[:250]}
        finally:dest.unlink(missing_ok=True)
        results.append(r)
    out={"status":"PUBLIC_SINGLE_CELL_DESIGN_GATE_NOT_EXPRESSION","source":"GSE313736","files":results,
         "limits":["No donor clinical metadata, genotype VCF, individual identifiers, expression matrices or raw reads exported.",
                   "Exact donor/condition definitions and QC singlets must be checked before pseudobulk; 74 reported donors are not 74 proven usable pairs.",
                   "This source gate is not a replication result. Full frozen-panel scoring protocol must be fixed before accessing expression."]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n");print(json.dumps(out,indent=2,allow_nan=False))
if __name__=="__main__":main()
