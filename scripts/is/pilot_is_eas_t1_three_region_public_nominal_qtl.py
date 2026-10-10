#!/usr/bin/env python3
"""Public eQTL Catalogue nominal QTL pilot for three strong AIS/BBJ regions.

Downloads three high-support source groups (G0015, G0023, G0007), three public
nominal eQTL datasets, +/-500kb GRCh37 lead windows converted to GRCh38 using
UCSC hg19ToHg38 chain with coordinate/strand checks. Does not do coloc, does
not label absent QTL rows as biologically negative, never touch prior 646.
"""
import argparse,csv,gzip,hashlib,json,math,subprocess,sys,time
from collections import Counter,defaultdict
from pathlib import Path
DATASETS={
 "QTD000131":("GTEx_Artery_Aorta",387,"QTS000015"),
 "QTD000171":("GTEx_Brain_Cortex",205,"QTS000015"),
 "QTD000434":("ROSMAP_Brain_DLPFC",560,"QTS000025"),
}
REGIONS=("IS_XDATA_G0015","IS_XDATA_G0023","IS_XDATA_G0007")
def read(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(4194304),b""):h.update(b)
 return h.hexdigest()
def write(p,rows):
 with p.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
  w.writeheader();w.writerows(rows)
def lift_interval(lo,r):
 chrom="chr"+str(r["chr"])
 positions=[max(1,int(r["start"])-500000),int(r["start"]),int(r["end"]),int(r["end"])+500000]
 hits=[]
 for bp in positions:
  got=lo.convert_coordinate(chrom,bp-1)
  if len(got)!=1 or got[0][0]!=chrom or got[0][2]!="+":
   raise ValueError("Ambiguous GRCh37/38 boundary "+r["group_id"]+":"+str(bp))
  hits.append(int(got[0][1])+1)
 if hits!=sorted(hits):
  raise ValueError("Nonmonotone liftOver interval; do not query as contiguous")
 return chrom.replace("chr",""),hits[0],hits[-1],hits
def get_url(dataset):
 q=DATASETS[dataset][2]
 return f"https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/{q}/{dataset}/{dataset}.all.tsv.gz"
def fetch(dataset,region,path):
 url=get_url(dataset)
 if path.exists():raise FileExistsError("Do not overwrite existing QTL slice")
 p=subprocess.run(["tabix",url,region],capture_output=True,timeout=95)
 if p.returncode:
  raise RuntimeError("Remote tabix failed: "+p.stderr.decode("utf-8","replace")[-250:])
 if len(p.stdout)>100*1024*1024:raise ValueError("QTL region over 100MB; restrict interval")
 if len(p.stdout)==0:return None
 temp=path.with_name(path.name+".part")
 if temp.exists():raise FileExistsError("Partial output from old run: "+str(temp))
 try:
  with gzip.GzipFile(filename="",mode="wb",fileobj=temp.open("wb"),mtime=0,compresslevel=1) as f:
   f.write(p.stdout)
  temp.replace(path)
 finally:
  if temp.exists():temp.unlink()
 return path
def audit_slice(path,chrom,start,end,genes):
 source=Counter();hit=defaultdict(lambda:{"records":0,"nominal_sig":0,"minp":1.0})
 with gzip.open(path,"rt") as f:
  for ln in f:
   if not ln.strip():continue
   col=ln.rstrip("\r\n").split("\t")
   source["rows"]+=1
   if len(col)!=19:raise ValueError("eQTL Catalogue public nominal schema not 19 columns")
   if col[1]!=chrom or not start<=int(col[2])<=end:
    raise ValueError("Retrieved variant outside requested GRCh38 region")
   if col[5]!=f"chr{chrom}_{col[2]}_{col[3]}_{col[4]}":
    raise ValueError("QTL variant/ref/alt inconsistency")
   source["unique_variant_gene_records"]+=1
   gene=col[16].split(".")[0]
   if gene not in genes:continue
   source["in_positional_genes"]+=1
   target=hit[gene];target["records"]+=1
   try:p=float(col[8])
   except ValueError:continue
   if math.isfinite(p) and 0<=p<=1:
    target["minp"]=min(p,target["minp"])
    target["nominal_sig"]+=p<.05
    source["target_valid_p"]+=1
   else:source["target_invalid_p"]+=1
 return dict(source),dict(hit)
def main():
 ap=argparse.ArgumentParser()
 for opt in ("regions","gene-universe","chain","vendor-pyliftover","raw-out","audit-out"):
  ap.add_argument("--"+opt,required=True,type=Path)
 a=ap.parse_args()
 for p in (a.regions,a.gene_universe,a.chain):
  if not p.is_file():raise FileNotFoundError(p)
 if a.raw_out.exists() and any(a.raw_out.iterdir()):
  raise FileExistsError("No overwrite of QTL pilot raw directory")
 if a.audit_out.exists() and any(a.audit_out.iterdir()):
  raise FileExistsError("No overwrite of QTL pilot audit directory")
 sys.path.insert(0,str(a.vendor_pyliftover))
 from pyliftover import LiftOver
 lo=LiftOver(str(a.chain))
 region_records={r["group_id"]:r for r in read(a.regions)}
 if not all(r in region_records for r in REGIONS):raise ValueError("Expected two pilot regions missing")
 positional=read(a.gene_universe)
 if len(positional)!=2225:raise ValueError("Original 2225 gene universe changed")
 genes_byregion={}
 for rid in REGIONS:
  genes_byregion[rid]={g["gene_id_stable"] for g in positional
   if rid in g["region_ids"].split(";")}
 if {k:len(v) for k,v in genes_byregion.items()}!={"IS_XDATA_G0015":71,"IS_XDATA_G0023":14,"IS_XDATA_G0007":10}:
  raise ValueError("Prespecified region gene universe changed")
 a.raw_out.mkdir(parents=True,exist_ok=True)
 a.audit_out.mkdir(parents=True,exist_ok=True)
 summaries=[];genecover=[]
 for rid in REGIONS:
  rr=region_records[rid]
  chrom,s,e,mapped=lift_interval(lo,rr)
  for did,(label,n,_) in DATASETS.items():
   qregion=f"{chrom}:{s}-{e}"
   path=a.raw_out/f"{rid}__{did}__GRCh38_{chrom}_{s}_{e}.tsv.gz"
   summary={"region_id":rid,"dataset_id":did,"study":label,
    "catalogue_n":n,"original_region_GRCh37":f'{rr["chr"]}:{rr["start"]}-{rr["end"]}',
    "input_GRCh37_flank_bp":500000,
    "queried_region_GRCh38":qregion,
    "liftOver_anchor_checks":";".join(map(str,mapped)),
    "n_region_positional_genes":len(genes_byregion[rid]),
    "url":get_url(did),"status":"","regional_qtl_rows":0,
    "target_gene_qtl_rows":0,"target_genes_with_any_rows":0,
    "valid_target_QTL_p_rows":0,"raw_size_bytes":0,"raw_sha256":"",
    "error":"","scientific_status":"NOMINAL_QTL_SOURCE_COVERAGE_NOT_COLOC"}
   try:
    got=fetch(did,qregion,path)
    if got is None:
     summary["status"]="EMPTY_QTL_REGION_NOT_A_NEGATIVE_GENE_RESULT"
    else:
     counters,per_gene=audit_slice(path,chrom,s,e,genes_byregion[rid])
     summary.update({
      "status":"PUBLIC_NOMINAL_QTL_REGION_DOWNLOADED",
      "regional_qtl_rows":counters.get("rows",0),
      "target_gene_qtl_rows":counters.get("in_positional_genes",0),
      "target_genes_with_any_rows":len(per_gene),
      "valid_target_QTL_p_rows":counters.get("target_valid_p",0),
      "raw_size_bytes":path.stat().st_size,
      "raw_sha256":sha(path)})
     for gene in sorted(genes_byregion[rid]):
      vals=per_gene.get(gene,{"records":0,"nominal_sig":0,"minp":None})
      genecover.append({
       "region_id":rid,"dataset_id":did,"gene_id":gene,
       "n_raw_QTL_records":vals["records"],
       "n_nominal_p_lt_0p05":vals["nominal_sig"],
       "min_QTL_p":vals["minp"] if vals["records"] else "",
       "coverage_status":"NOMINAL_RECORDS_IN_RETRIEVED_WINDOW" if vals["records"] else
                         "NO_NOMINAL_RECORD_IN_THIS_SOURCE_WINDOW_NOT_NEGATIVE",
       "coloc_status":"NOT_RUN_NO_FULL_ALLELE_HARMONIZATION"
      })
   except Exception as exc:
    summary["status"]="RETRIEVAL_OR_SCHEMA_FAILURE_NOT_DATA_ABSENCE"
    summary["error"]=str(exc)[:500]
   summaries.append(summary)
   print("QTL_PILOT_SOURCE",rid,did,summary["status"],
    "allrows",summary["regional_qtl_rows"],"targetgenes",
    summary["target_genes_with_any_rows"],"error",summary["error"][:150],flush=True)
   time.sleep(3)
 write(a.audit_out/"IS_EAS_T1_THREE_REGION_THREE_PUBLIC_QTL_SOURCE_COVERAGE.tsv",summaries)
 if genecover:write(a.audit_out/"IS_EAS_T1_THREE_REGION_GENE_NOMINAL_QTL_COVERAGE.tsv",genecover)
 stats={
  "status":"PUBLIC_NOMINAL_QTL_PILOT_NO_COLOC",
  "regions":list(REGIONS),"QTL_datasets":list(DATASETS),
  "expected_slices":9,"source_successes":sum(r["status"]=="PUBLIC_NOMINAL_QTL_REGION_DOWNLOADED" for r in summaries),
  "source_failures":sum(r["status"]=="RETRIEVAL_OR_SCHEMA_FAILURE_NOT_DATA_ABSENCE" for r in summaries),
  "per_gene_nominal_counts_not_tissue_specific_functional_validation":True,
  "QTL_ancestry_verified_for_BBJ":False,"source_GWAS_cohort_independence_verified":False,
  "genome_build":"GRCh38","no_variant_harmonization_or_coloc":True,
  "input_sha256":{str(p):sha(p) for p in (a.regions,a.gene_universe,a.chain)}
 }
 (a.audit_out/"IS_EAS_T1_QTL_PILOT_STATUS.json").write_text(json.dumps(stats,indent=2)+"\n")
 print("IS_EAS_T1_THREE_REGION_QTL_SOURCE_PILOT_COMPLETE",json.dumps(stats),flush=True)
 if stats["source_successes"]<1:sys.exit(2)
if __name__=="__main__":main()
