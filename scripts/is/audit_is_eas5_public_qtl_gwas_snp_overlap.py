#!/usr/bin/env python3
"""Read-only 5-region nominal QTL -> BBJ/GIGASTROKE EAS AIS GRCh37 SNP overlap.

GRCh38 QTL source -> chain position/strand -> independent GRCh37 FASTA
REF check -> exact, unordered allele-aware GWAS summary intersection.
GWAS beta/eaf oriented to the GRCh37 ALT; QTL beta effect allele provenance
NOT verified and is NOT reoriented or qualified for coloc. Context +/-500kb
is a discovery window, NOT an independent LD-defined causal locus.
"""
import argparse,csv,gzip,hashlib,json,math,sys
from collections import Counter,defaultdict,OrderedDict
from pathlib import Path
TRAITS={"BBJ_IS","EAS_AIS"}
STUDIES={"QTD000131","QTD000171","QTD000434"}
REGIONS=("IS_XDATA_G0015","IS_XDATA_G0023","IS_XDATA_G0007","IS_XDATA_G0025","IS_XDATA_G0018")
COMP=str.maketrans("ACGT","TGCA")
def read(path):
 with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def save(path,rows):
 with path.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
  w.writeheader();w.writerows(rows)
def sha(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for b in iter(lambda:f.read(4194304),b""):h.update(b)
 return h.hexdigest()
def number(v):
 try:
  x=float(v)
  return x if math.isfinite(x) else None
 except (ValueError,TypeError):return None
class Fasta37:
 def __init__(self,path):
  self.fh=path.open("rb")
  with Path(str(path)+".fai").open() as f:
   self.index={}
   for line in f:
    cols=line.rstrip("\n").split("\t")
    self.index[cols[0]]=tuple(map(int,cols[1:5]))
  self.cache=OrderedDict()
 def get(self,chrom,pos):
  if chrom not in self.index:return None
  length,offset,bases,bytewidth=self.index[chrom]
  if not 1<=pos<=length:return None
  first=((pos-1)//65536)*65536+1
  key=(chrom,first)
  if key not in self.cache:
   end=min(first+65535,length)
   def loc(p):
    p-=1
    return offset+(p//bases)*bytewidth+(p%bases)
   self.fh.seek(loc(first))
   raw=self.fh.read(loc(end)-loc(first)+1)
   data=raw.replace(b"\n",b"").replace(b"\r",b"").upper()
   if len(data)!=end-first+1:raise ValueError("FASTA FAI offset or line-width drift")
   self.cache[key]=data
   if len(self.cache)>96:self.cache.popitem(last=False)
  else:self.cache.move_to_end(key)
  return chr(self.cache[key][pos-first])
 def close(self):self.fh.close()
def canonical_qtl(chrom,pos,ref,alt,lo,fa,cache):
 key=(chrom,pos,ref,alt)
 if key in cache:return cache[key]
 if not (ref in "ACGT" and alt in "ACGT" and len(ref)==len(alt)==1 and ref!=alt):
  return "QTL_NOT_BIALLELIC_SNP","",None
 hit=lo.convert_coordinate("chr"+chrom,int(pos)-1)
 if len(hit)!=1:
  result=("LIFTOVER_NOT_UNIQUE","",None)
 elif hit[0][0]!="chr"+chrom or hit[0][1]!=int(hit[0][1]) or hit[0][1]<0:
  result=("LIFTOVER_NON_INTEGER_OR_CHR_SHIFT","",None)
 else:
  chrom37,point,strand,_=hit[0]
  r,a=ref,alt
  if strand=="-":r,a=ref.translate(COMP),alt.translate(COMP)
  elif strand!="+":return "BAD_LIFTOVER_STRAND","",None
  refbase=fa.get(chrom,int(point)+1)
  if refbase is None:result=("FASTA37_REFERENCE_NOT_FOUND","",None)
  elif refbase!=r:
   result=("GRCH37_FASTA_REF_DISAGREES",f"{chrom}:{int(point)+1}:{r}:{a}",None)
  else:
   genomic=f"{chrom}:{int(point)+1}:{r}:{a}"
   result=("GRCH37_FASTA_REF_VERIFIED",genomic,strand)
 cache[key]=result
 return result
def checked_gwas_row(row,trait,fa):
 try:
  chrom=row["chr"].replace("chr","")
  pos=int(row["pos"])
  variant=row["variant_id"].upper()
  parts=variant.split(":")
  if len(parts)!=4 or parts[0]!=chrom or int(parts[1])!=pos:
   return "ID_CHR_POS_INVALID",None
  ref,alt=parts[2:]
  if len(ref)!=1 or len(alt)!=1 or ref not in "ACGT" or alt not in "ACGT" or ref==alt:
   return "GWAS_NOT_BIALLELIC_SNP",None
  if trait=="BBJ_IS" and (row["ref"]!=ref or row["alt"]!=alt):
   return "GWAS_SOURCE_REF_ALT_MISMATCH",None
  if fa.get(chrom,pos)!=ref:return "GWAS_GRCH37_FASTA_REF_MISMATCH",None
  ea=row["effect_allele"].upper()
  oa=row["other_allele"].upper()
  if not ((ea==alt and oa==ref) or (ea==ref and oa==alt)):
   return "GWAS_EFFECT_ALLELE_SET_MISMATCH",None
  b,s,p,f=(number(row[x]) for x in ("beta","se","p","eaf"))
  if None in (b,s,p,f) or s<=0 or not 0<=p<=1 or not 0<=f<=1:
   return "GWAS_EFFECT_P_OR_AF_INVALID",None
  oriented=(b if ea==alt else -b)
  aligned_freq=f if ea==alt else 1-f
  return "VALID_GWAS_SOURCE_REF_ALT_EFFECT",{
   "variant":variant,"beta_alt":oriented,"se":s,"p":p,
   "af_alt":aligned_freq,"palindrome":frozenset((ref,alt)) in (
   frozenset(("A","T")),frozenset(("C","G")))
  }
 except (ValueError,IndexError,KeyError):return "GWAS_ROW_INVALID",None
def load_gwas(path,trait,groups,fa):
 index=defaultdict(dict)
 stats=Counter()
 lookup=defaultdict(list)
 for rid,r in groups.items():lookup[r["chr"]].append((rid,
  max(1,int(r["start"])-500000),int(r["end"])+500000))
 with gzip.open(path,"rt") as f:
  rdr=csv.DictReader(f,delimiter="\t")
  if not {"chr","pos","variant_id","effect_allele","other_allele","beta","se","p","eaf"}.issubset(rdr.fieldnames):
   raise ValueError("GWAS source canonical header incomplete")
  for row in rdr:
   stats["all_records_scanned"]+=1
   chrom=row["chr"].replace("chr","")
   if chrom not in lookup:continue
   try:p=int(row["pos"])
   except ValueError:stats["position_invalid"]+=1;continue
   regions=[r for r,lo,hi in lookup[chrom] if lo<=p<=hi]
   if not regions:continue
   stats["context_rows"]+=1
   status,record=checked_gwas_row(row,trait,fa)
   stats[status]+=1
   if record:
    for rid in regions:
     prev=index[rid].get(record["variant"])
     if prev:
      if prev!=record:raise ValueError("Inconsistent duplicate source GWAS variant")
      stats["identical_duplicate_variant"]+=1
     else:index[rid][record["variant"]]=record
   if stats["all_records_scanned"]%5000000==0:
    print("SOURCE_GWAS_PROGRESS",trait,stats["all_records_scanned"],flush=True)
 if stats["all_records_scanned"]<5000000:raise ValueError("Suspiciously incomplete GWAS source")
 print("SOURCE_GWAS_LOADED",trait,"all",stats["all_records_scanned"],
       "context",stats["context_rows"],"valid",stats["VALID_GWAS_SOURCE_REF_ALT_EFFECT"],
       "indexed",sum(map(len,index.values())),flush=True)
 return index,dict(stats)
def qtl_file(source,rootT1,rootT4):
 rid=source["region_id"];did=source["dataset_id"]
 if rid not in REGIONS or did not in STUDIES:raise ValueError("Unexpected pilot source")
 chrom,coords=source["queried_region_GRCh38"].split(":")
 start,end=coords.split("-")
 root=rootT1 if rid in REGIONS[:3] else rootT4
 return root/f"{rid}__{did}__GRCh38_{chrom}_{start}_{end}.tsv.gz"
def validate_reverse_chain(path):
 if "hg38ToHg19" not in path.name:
  raise ValueError("QTL GRCh38 to GWAS GRCh37 requires hg38ToHg19 chain")
 with path.open("rt") as f:
  hdr=f.readline().split()
 # UCSC chain header target chromosome sizes distinguish GRCh38 from GRCh37.
 if len(hdr)!=13 or hdr[0]!="chain" or hdr[2]!="chr1" or hdr[3]!="248956422" or hdr[7]!="chr1" or hdr[8]!="249250621":
  raise ValueError("Chain target GRCh38/query GRCh37 assembly headers not validated")
 return True
def run(a):
 validate_reverse_chain(a.chain)
 from pyliftover import LiftOver
 if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError("No existing output")
 source=read(a.source_manifest)
 genes=read(a.gene_universe)
 regions=read(a.regions)
 if len(source)!=15 or len(genes)!=2225 or len(regions)!=30:
  raise ValueError("Source pilot or positional gene universe count drift")
 region={r["group_id"]:r for r in regions if r["group_id"] in REGIONS}
 if set(region)!=set(REGIONS):raise ValueError("5 source regions changed")
 gene_lookup={x["gene_id_stable"]:x for x in genes}
 source_checks={}
 for x in source:
  p=qtl_file(x,a.t1_raw,a.technical_raw)
  if not p.is_file() or sha(p)!=x["raw_sha256"]:
   raise ValueError("Raw public QTL source changed or missing "+str(p))
  source_checks[str(p)]=x["raw_sha256"]
 lo=LiftOver(str(a.chain));fa=Fasta37(a.fasta)
 print("GWAS_LOADING_START",flush=True)
 bindex,bstats=load_gwas(a.bbj,"BBJ_IS",region,fa)
 aindex,astats=load_gwas(a.ais,"EAS_AIS",region,fa)
 a.out_dir.mkdir(parents=True,exist_ok=True)
 cache={};summary=[];gene_stats=defaultdict(Counter)
 matchpath=a.out_dir/"IS_QTL_5REGION_EXACT_GWAS_VARIANT_MATCHES.tsv.gz"
 fields=["region_id","QTL_dataset","gene_id","source_variant_GRCh38",
  "mapped_variant_GRCh37","mapped_strand","QTL_p","QTL_beta_SOURCE_ORIENTATION_UNVERIFIED",
  "QTL_se","QTL_maf","QTL_ac","QTL_an",
  "BBJ_is_present","BBJ_beta_ALT","BBJ_p","BBJ_ALT_EAF",
  "EAS_AIS_is_present","EAS_AIS_beta_ALT","EAS_AIS_p","EAS_AIS_ALT_EAF",
  "palindromic_SNP_review","QTL_effect_allele_provenance",
  "science_gate"]
 with gzip.open(matchpath,"wt",newline="") as h:
  w=csv.DictWriter(h,fieldnames=fields,delimiter="\t");w.writeheader()
  for s in source:
   rid,did=s["region_id"],s["dataset_id"]
   path=qtl_file(s,a.t1_raw,a.technical_raw)
   c=Counter()
   with gzip.open(path,"rt") as f:
    for line in f:
     if not line.strip():continue
     c["source_rows"]+=1
     parts=line.rstrip("\r\n").split("\t")
     if len(parts)!=19:raise ValueError("QTL nominal 19-column schema drift")
     gene=parts[16].split(".")[0]
     if gene not in gene_lookup or rid not in gene_lookup[gene]["region_ids"].split(";"):
      c["non_target_gene"]+=1;continue
     g=gene_stats[(rid,did,gene)]
     g["QTL_raw_gene_variant_records"]+=1
     if parts[1]!=region[rid]["chr"] or parts[5]!=f"chr{parts[1]}_{parts[2]}_{parts[3]}_{parts[4]}":
      raise ValueError("Original QTL variant tuple and chromosome drift")
     if parts[11]!="SNP" or len(parts[3])!=1 or len(parts[4])!=1:
      g["QTL_non_SNP"]+=1;c["non_SNP"]+=1;continue
     g["QTL_SNP_rows"]+=1
     status,mapped,strand=canonical_qtl(parts[1],int(parts[2]),
       parts[3].upper(),parts[4].upper(),lo,fa,cache)
     g[status]+=1;c[status]+=1
     if status!="GRCH37_FASTA_REF_VERIFIED":continue
     qbeta=number(parts[9]);qse=number(parts[10]);qp=number(parts[8])
     qmaf=number(parts[7]);qac=number(parts[12]);qan=number(parts[13])
     if None in (qbeta,qse,qp,qmaf,qac,qan) or qse<=0 or not 0<=qp<=1 or not 0<=qmaf<=.5 or qan<=0 or not 0<=qac<=qan:
      g["QTL_STATS_INVALID"]+=1;continue
     if abs(min(qac/qan,1-qac/qan)-qmaf)>2e-4:
      g["QTL_AC_AN_MAF_INCONSISTENT"]+=1;continue
     g["QTL_STATS_VALID"]+=1
     bb=bindex[rid].get(mapped);ais=aindex[rid].get(mapped)
     # Both source SNP records are already independently verified against GRCh37 FASTA.
     if bb:g["EXACT_REF_ALT_BBJ"]+=1
     if ais:g["EXACT_REF_ALT_EAS_AIS"]+=1
     if not bb and not ais:continue
     g["EXACT_EITHER_GWAS"]+=1;c["EXACT_EITHER_GWAS"]+=1
     if bb and ais:g["EXACT_BOTH_GWAS"]+=1
     pal=frozenset(mapped.split(":")[2:]) in (
         frozenset(("A","T")),frozenset(("C","G")))
     if pal:g["PALINDROMIC_REVIEW"]+=1
     if bb and ais and not pal:g["EXACT_BOTH_GWAS_NONPALINDROMIC"]+=1
     w.writerow({
      "region_id":rid,"QTL_dataset":did,"gene_id":gene,
      "source_variant_GRCh38":parts[5],"mapped_variant_GRCh37":mapped,
      "mapped_strand":strand,"QTL_p":qp,
      "QTL_beta_SOURCE_ORIENTATION_UNVERIFIED":qbeta,
      "QTL_se":qse,"QTL_maf":qmaf,"QTL_ac":qac,"QTL_an":qan,
      "BBJ_is_present":bool(bb),"BBJ_beta_ALT":bb["beta_alt"] if bb else "",
      "BBJ_p":bb["p"] if bb else "","BBJ_ALT_EAF":bb["af_alt"] if bb else "",
      "EAS_AIS_is_present":bool(ais),"EAS_AIS_beta_ALT":ais["beta_alt"] if ais else "",
      "EAS_AIS_p":ais["p"] if ais else "","EAS_AIS_ALT_EAF":ais["af_alt"] if ais else "",
      "palindromic_SNP_review":pal,
      "QTL_effect_allele_provenance":"UNVERIFIED_SOURCE_BETA_ALLELE",
      "science_gate":"SNP_OVERLAP_ONLY_NO_COLOC_NO_QTL_EFFECT_ALIGNMENT"
     })
   if c["source_rows"]!=int(s["regional_qtl_rows"]):
    raise ValueError("Raw source row count drift: "+did+" "+rid)
   summary.append({
    "region_id":rid,"QTL_dataset":did,"QTL_source_n":s["catalogue_n"],
    "source_rows":c["source_rows"],"target_gene_rows":c["source_rows"]-c["non_target_gene"],
    "target_genes_with_source_QTL_rows":len({g for r,ds,g in gene_stats if r==rid and ds==did}),
    "source_SNP_rows":c["source_rows"]-c["non_target_gene"]-c["non_SNP"],
    "unique_source_variants_liftover_cache_global":len(cache),
    "GRCH37_FASTA_REF_VERIFIED":c["GRCH37_FASTA_REF_VERIFIED"],
    "GRCH37_FASTA_REF_DISAGREES":c["GRCH37_FASTA_REF_DISAGREES"],
    "LIFTOVER_NOT_UNIQUE":c["LIFTOVER_NOT_UNIQUE"],
    "LIFTOVER_NON_INTEGER_OR_CHR_SHIFT":c["LIFTOVER_NON_INTEGER_OR_CHR_SHIFT"],
    "exact_any_GWAS_geneSNP_rows":c["EXACT_EITHER_GWAS"],
    "source_published_variant_id":"GRCh38",
    "GRCh37_independent_FASTA_verified":True,
    "native_QTL_beta_effect_allele_direction_verified":False,
    "independent_population_replication_established":False
   })
   print("QTL_GWAS_OVERLAP_FILE_DONE",rid,did,"target",summary[-1]["target_gene_rows"],
         "liftOK",c["GRCH37_FASTA_REF_VERIFIED"],"matchedEither",c["EXACT_EITHER_GWAS"],flush=True)
 records=[]
 for (rid,did,gene),x in sorted(gene_stats.items()):
  source_row=gene_lookup[gene]
  records.append({
   "region_id":rid,"QTL_dataset":did,"gene_id":gene,
   "gene_symbol":source_row["canonical_GENCODE_v19_gene_symbol"],
   "has_archived_GTEx_ABF":int(source_row["new_DIRECT_prior_grid_assays"])>0,
   "QTL_source_nominal_rows":x["QTL_raw_gene_variant_records"],
   "source_SNP_rows":x["QTL_SNP_rows"],
   "GRCh37_FASTA_VERIFIED_SNP_rows":x["GRCH37_FASTA_REF_VERIFIED"],
   "FASTA_REF_CONFLICT_rows":x["GRCH37_FASTA_REF_DISAGREES"],
   "QTL_valid_stats_rows":x["QTL_STATS_VALID"],
   "BBJ_exact_GRCh37_variant_rows":x["EXACT_REF_ALT_BBJ"],
   "EAS_AIS_exact_GRCh37_variant_rows":x["EXACT_REF_ALT_EAS_AIS"],
   "both_GWAS_exact_variant_rows":x["EXACT_BOTH_GWAS"],
   "both_GWAS_exact_nonpalindromic_rows":x["EXACT_BOTH_GWAS_NONPALINDROMIC"],
   "palindromic_GWAS_matched_rows":x["PALINDROMIC_REVIEW"],
   "eligible_for_coloc":False,
   "blocked_on":"QTL_BETA_EFFECT_ALLELE_UNVERIFIED;NO_MOLECULAR_LD;QC_AND_INDEPENDENCE"
  })
 save(a.out_dir/"IS_EAS5_GENE_TISSUE_GWAS_SNP_OVERLAP_QC.tsv",records)
 save(a.out_dir/"IS_EAS5_SOURCE_QTL_GWAS_SNP_OVERLAP_QC.tsv",summary)
 for s in summary:
  if s["target_gene_rows"]!=int(next(z["target_gene_qtl_rows"] for z in source if z["region_id"]==s["region_id"] and z["dataset_id"]==s["QTL_dataset"])):
   raise ValueError("Target source row count drift")
 totals=Counter()
 for row in records:
  for k in ("QTL_source_nominal_rows","GRCh37_FASTA_VERIFIED_SNP_rows",
   "BBJ_exact_GRCh37_variant_rows","EAS_AIS_exact_GRCh37_variant_rows",
   "both_GWAS_exact_variant_rows","both_GWAS_exact_nonpalindromic_rows"):
   totals[k]+=row[k]
 verified_gene={r["gene_id"] for r in records if r["both_GWAS_exact_nonpalindromic_rows"]>0}
 nominal_source_genes={r["gene_id"] for r in records}
 summaryobj={
  "status":"PUBLIC_NOMINAL_QTL_15_SLICES_GRCH37_FASTA_AND_GWAS_ALLELE_EXACT_OVERLAP_ONLY",
  "original_2225_candidate_genes_preserved":len(genes),
  "original_5_regions":list(REGIONS),
  "QTL_source_slices":len(source),
  "QTL_genes_present_in_raw_sources":len(nominal_source_genes),
  "genes_with_at_least_one_BBJ_AND_EAS_AIS_exact_nonpalindromic_SNP":len(verified_gene),
  "gene_tissue_records":len(records),"mapped_distinct_38_variant_cache_entries":len(cache),
  "quantities":dict(totals),"GWAS_read_counts":{"BBJ":bstats,"GWS_EAS_AIS":astats},
  "original_6_outcome_source_gwas_not_modified":True,
  "QTL_beta_effect_allele_independent_provenance_verified":False,
  "BBJ_and_GWS_independent_replication_established":False,
  "all_new_coloc_results_blocked":True,
  "QTL_raw_sha256":source_checks,
  "GWAS_and_build_input_SHA256":{str(p):sha(p) for p in (
    a.bbj,a.ais,a.regions,a.gene_universe,a.chain)}
 }
 (a.out_dir/"IS_EAS5_NOMINAL_QTL_SNP_OVERLAP_MANIFEST.json").write_text(json.dumps(summaryobj,indent=2)+"\n")
 print("IS_EAS5_GWAS_QTL_ALLELE_AND_FASTA_OVERLAP_AUDIT_COMPLETE",json.dumps({
  k:summaryobj[k] for k in ("QTL_source_slices","QTL_genes_present_in_raw_sources",
   "genes_with_at_least_one_BBJ_AND_EAS_AIS_exact_nonpalindromic_SNP","gene_tissue_records","quantities")}),flush=True)
 fa.close()
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for k in ("source-manifest","gene-universe","regions","bbj","ais",
           "chain","fasta","t1-raw","technical-raw","out-dir","vendor"):
  p.add_argument("--"+k,required=True,type=Path)
 a=p.parse_args()
 if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError("No output overwrite")
 sys.path.insert(0,str(a.vendor))
 run(a)
