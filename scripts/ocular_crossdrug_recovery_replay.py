#!/usr/bin/env python3
"""Focused ocular replay: compact GEO platform, paired FASTQ layout, alternate count coverage."""
import argparse,csv,gzip,hashlib,json,math,re,shutil,urllib.request,zipfile
from pathlib import Path
from publication_ocular_ckd_sources import ocular,SAMPLES,digest
from ocular_paired_panel_validation import ena_design,BUNDLE_SHA

def write(out,r):
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(r,indent=2,allow_nan=False)+"\n")
    print(json.dumps(r,indent=2,allow_nan=False),flush=True)

def crossdrug(root,out):
    p=root/"platform.txt"
    # Platform-only SOFT response, not the full family archive; ensure exact accession and table.
    with p.open("rt",encoding="utf-8-sig",errors="strict") as f:
        head=[]
        for line in f:
            if line.strip().lower()=="!platform_table_begin":break
            head.append(line.strip())
            if len(head)>10000:raise ValueError("Platform preamble exceeds bound")
        else:raise ValueError("No platform table")
        if "^PLATFORM = GPL6480" not in head:raise ValueError("Wrong platform accession")
        h=f.readline().rstrip("\r\n").split("\t")
        if not {"ID","CONTROL_TYPE","GENE_SYMBOL"}<=set(h):raise ValueError("Unexpected platform source columns "+repr(h))
    # Reuse the already frozen scoring function unchanged.
    with p.open("rb") as f,gzip.open(root/"platform.gz","wb") as w:shutil.copyfileobj(f,w)
    r=ocular(root)
    r["platform_text_sha256"]=digest(p)
    r["platform_source"]="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6480&targ=self&form=text&view=full"
    r["limitations"]=["Two biological TM strains, not 23 independent donors. No inferential p value or clinical response claim.",
        "2020 integrated analysis reported a GSE16643 donor excluded as QC outlier; both strains retained and shown separately here.",
        "Cross-drug response already studied in original paper. This is frozen-panel portability, not novel discovery or proof of drug equivalence.",
        "Percentile shifts depend on all complete non-control probes; relative ranks not log2 fold changes."]
    for panel,d in r["panels"].items():
        d["same_direction_positive_strain_drug_combinations"]=sum(x["mean_signed_percentile_shift"]>0 for x in d["scores"].values())
        d["all_six_positive"]=d["same_direction_positive_strain_drug_combinations"]==6
    write(out,r)

def recovery(root,out):
    design=ena_design(root/"ena.tsv")
    b=root/"counts.zip"
    if digest(b)!=BUNDLE_SHA:raise ValueError("DEE2 bundle changed")
    with zipfile.ZipFile(b) as z:
        if z.testzip():raise ValueError("DEE2 archive CRC")
        import io
        with z.open("SRP319815_NA/MetadataSummary.tsv") as f:
            md=list(csv.DictReader(io.TextIOWrapper(f),delimiter="\t"))
    available={r["SRR_accession"]:r for r in md}
    if len(available)!=len(md):raise ValueError("Duplicate DEE2 run")
    with (root/"fastq.tsv").open() as f:
        raw=list(csv.DictReader(f,delimiter="\t"))
    entries={r["run_accession"]:r for r in raw}
    expected={r["run_accession"] for arms in design.values() for r in arms.values()}
    if set(entries)!=expected or len(entries)!=len(raw) or not set(available)<=expected:raise ValueError("Manifest identities changed")
    runs=[]
    for donor,arms in sorted(design.items()):
        for arm,a in sorted(arms.items()):
            run=a["run_accession"];r=entries[run]
            urls=r["fastq_ftp"].split(";");md5=r["fastq_md5"].split(";");size=r["fastq_bytes"].split(";")
            if r["library_layout"]!="PAIRED" or len(urls)!=len(md5) or len(urls)!=len(size) or len(set(urls))!=len(urls):
                raise ValueError("Invalid declared run manifest "+run)
            files={}
            for u,m,n in zip(urls,md5,size):
                if not u.startswith("ftp.sra.ebi.ac.uk/") or not re.fullmatch("[a-fA-F0-9]{32}",m) or not n.isdecimal() or int(n)<=0:
                    raise ValueError("Unexpected ENA file fields")
                name=u.rsplit("/",1)[-1]
                if name in files:raise ValueError("Repeated filename")
                files[name]={"file":name,"url":"https://"+u,"declared_md5":m,"declared_bytes":int(n)}
            pair=[run+"_1.fastq.gz",run+"_2.fastq.gz"]
            singles=sorted(set(files)-set(pair))
            status="EXACT_PAIR_PLUS_UNSUFFIXED" if set(pair)<=set(files) and singles==[run+".fastq.gz"] else "EXACT_PAIR_ONLY" if set(files)==set(pair) else "HOLD_LAYOUT"
            runs.append({"donor":donor,"arm":arm,"source_group":a["experiment_title"],"run":run,"present_in_DEE2":run in available,
                "status":status,"paired_files":[files[k] for k in pair if k in files],"additional_files":[files[k] for k in singles],
                "all_declared_bytes":sum(x["declared_bytes"] for x in files.values())})
    summary={"status":"ENA_DECLARED_LAYOUT_NOT_FASTQ_CONTENT_VERIFIED","runs":runs,
      "source_manifest_sha256":digest(root/"fastq.tsv"),"source_runs":len(runs),
      "declared_bytes_all_runs":sum(r["all_declared_bytes"] for r in runs),
      "declared_bytes_missing_runs":sum(r["all_declared_bytes"] for r in runs if not r["present_in_DEE2"]),
      "exact_pair_plus_unsuffixed":sum(r["status"]=="EXACT_PAIR_PLUS_UNSUFFIXED" for r in runs),
      "exact_pair_only":sum(r["status"]=="EXACT_PAIR_ONLY" for r in runs),
      "limits":["Declared filenames and MD5 are not content or paired-read integrity validation.",
          "Unsuffixed file may contain unpaired reads; do not call it a third mate or discard without checking library layout and read IDs.",
          "No raw reads transferred, no new quantification. All 16 runs must use a uniform pipeline if recovered."]}
    write(out,summary)

def recount(out):
    url="https://duffel.rail.bio/recount3/human/data_sources/sra/metadata/sra.recount_project.MD.gz"
    try:
        req=urllib.request.Request(url,headers={"User-Agent":"Public-ocular-research/1.0"})
        with urllib.request.urlopen(req,timeout=30) as f:raw=f.read(20000001)
        if len(raw)>20000000:raise ValueError("Public project catalog cap")
        rows=gzip.decompress(raw).decode("utf-8").splitlines()
        if not rows:raise ValueError("Empty catalog")
        header=rows[0]
        found=[r for r in rows[1:] if "SRP319815" in r.split("\t")]
        result={"status":"PUBLIC_RECOUNT3_PROJECT_CATALOG_CHECK","url":url,"source_sha256":hashlib.sha256(raw).hexdigest(),"header":header,"rows":len(rows)-1,"exact_project_rows":found,
            "limit":"Catalog coverage only. Project absence is not absence of raw RNA-seq. If present, verify all 16 runs and coverage units before scoring."}
    except (OSError,ValueError,EOFError,UnicodeError) as exc:
        result={"status":"HOLD_CATALOG_ACCESS","url":url,"reason":str(exc)[:300]}
    write(out,result)

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["crossdrug","recovery","recount"]);p.add_argument("--root",type=Path,default=Path("."));p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    if a.mode=="crossdrug":crossdrug(a.root,a.out)
    elif a.mode=="recovery":recovery(a.root,a.out)
    else:recount(a.out)
