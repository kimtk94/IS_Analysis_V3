#!/usr/bin/env python3
"""GTEx QTS000015 QTD000161 reprocessed cerebellar hemisphere cis universe audit.

No requester-pays access. Retrieve one compact permutation summary (n_variants)
and two distant tabix region queries from open EBI. The archived original BBJ
and GTEx summaries are not interchangeable with this reprocessing.
"""
import argparse,csv,gzip,json,hashlib,math,subprocess,time
from collections import Counter,defaultdict
from pathlib import Path

GENES={"ENSG00000138675":("FGF5","BBJ_IS_L001"),"ENSG00000197826":("C4orf22","BBJ_IS_L001"),
       "ENSG00000138172":("CALHM2","BBJ_IS_L002"),"ENSG00000107954":("NEURL","BBJ_IS_L002")}
TISSUE="GTEx_V8__Brain_Cerebellar_Hemisphere"
URL="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000161/QTD000161.all.tsv.gz"
COLS="gene_id chromosome position ref alt variant ma_samples maf pvalue beta se type ac an r2 molecular_trait_object_id molecular_trait_id median_tpm rsid".split()

def tsv(path):
 with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))

def get_permuted(p):
 with gzip.open(p,"rt") as f:
  rows=[x for x in csv.DictReader(f,delimiter="\t") if x["molecular_trait_object_id"] in GENES]
 if len(rows)!=4 or len({x["molecular_trait_object_id"] for x in rows})!=4:
  raise ValueError("Expected four unique gene-wide permutation rows")
 for x in rows:
  if x["molecular_trait_id"]!=x["molecular_trait_object_id"] or int(x["n_variants"])<1500:
   raise ValueError("Unexpected cis-variant count schema")
 return {x["molecular_trait_object_id"]:x for x in rows}

def parse_qtd_row(line):
 vals=line.rstrip("\n").split("\t")
 if len(vals)!=len(COLS):raise ValueError("Incorrect indexed summary column count "+str(len(vals)))
 r=dict(zip(COLS,vals))
 if not r["variant"].startswith("chr"+r["chromosome"]+"_"+r["position"]+"_"):
  raise ValueError("variant ID differs from indexed genome coordinates")
 if r["variant"]!="chr"+r["chromosome"]+"_"+r["position"]+"_"+r["ref"]+"_"+r["alt"]:
  raise ValueError("REF/ALT mismatch in indexed QTL file")
 return r

def fetch_region(region,out,timeout=150):
 """Bounded ONE tabix call per locus; atomic tabix filtered gene output."""
 if out.exists() and out.stat().st_size>0:return {"status":"EXISTING_PRESERVED","path":str(out)}
 tmp=out.with_suffix(".part")
 try:
  proc=subprocess.Popen(["timeout",str(timeout),"tabix",URL,region],
                       stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,bufsize=1,cwd=str(out.parent))
  stats=Counter()
  with tmp.open("w",newline="") as f:
   w=csv.DictWriter(f,fieldnames=COLS,delimiter="\t");w.writeheader()
   for line in proc.stdout:
    stats["all_region_rows"]+=1
    if stats["all_region_rows"]>2000000:
     proc.kill();raise ValueError("Region over 2M associations; refusing unchecked volume")
    r=parse_qtd_row(line)
    if r["gene_id"] in GENES:
     if (GENES[r["gene_id"]][1]=="BBJ_IS_L001" and region.startswith("4:")) or (GENES[r["gene_id"]][1]=="BBJ_IS_L002" and region.startswith("10:")):
      w.writerow(r);stats["priority_gene_rows"]+=1
  err=proc.stderr.read(1000);status=proc.wait()
  if status!=0:raise RuntimeError(f"tabix region {region} exit {status}: {err[:350]}")
  tmp.replace(out)
  return {"status":"FETCHED","path":str(out),**stats}
 finally:
  if tmp.exists():tmp.unlink()

def evaluate(permuted,archived,region_files,out):
 full=get_permuted(permuted)
 region_rows=[]
 for f in region_files:
  if not f.is_file():raise FileNotFoundError(f)
  region_rows+=tsv(f)
 raw_source=Path(archived)
 audit=[];summaries=[]
 for gene,(symbol,locus) in GENES.items():
  input_file=raw_source/f"{locus}__{TISSUE}__{gene}.tsv"
  archive=tsv(input_file)
  if len(archive)!= (1694 if locus=="BBJ_IS_L001" else 1583):
   raise ValueError("Original SNP archive cardinality changed")
  if len({x["variant_id_hg38"] for x in archive})!=len(archive):
   raise ValueError("Original GRCh38 SNP ID duplicated")
  catalogue=[x for x in region_rows if x["gene_id"]==gene]
  # Several dbSNP rsids can map to the same variant: keep distinct REF/ALT SNPs.
  variants=defaultdict(list)
  for r in catalogue:
   variants[(r["chromosome"],r["position"],r["ref"],r["alt"])].append(r)
  dup=sum(len(v)-1 for v in variants.values())
  source={tuple(x["variant_id_hg38"].split(":")):x for x in archive}
  matched=set(source)&set(variants)
  # catalogue is a re-analysis. Concordance is a diagnostic, NOT expected identity.
  same=[];freq=[];differences=[]
  for key in sorted(matched):
   q=variants[key][0];a=source[key]
   q_p=float(q["pvalue"]);a_p=float(a["eqtl_p"])
   if not (0<=q_p<=1 and 0<=a_p<=1):raise ValueError("Invalid nominal p-value")
   same.append((q_p,a_p))
   f1=float(q["maf"]);f2=float(a["eqtl_maf"])
   if math.isfinite(f1) and math.isfinite(f2):freq.append(abs(f1-f2))
   differences.append(dict(symbol=symbol,gene_id=gene,variant_hg38=":".join(key),
      archived_nominal_p=a_p,catalogue_nominal_p=q_p,
      archived_beta=float(a["eqtl_beta"]),catalogue_beta=float(q["beta"]),
      same_association_p_within_1e3=(abs(a_p-q_p)<.001),
      catalogue_molecular_n=float(q["an"])/2 if q["an"] not in ("NA",".","") else "",
      archived_molecular_n=a["eqtl_n"],
      scientific_gate="CATALOGUE_REANALYSIS_DO_NOT_REPLACE_ARCHIVED_QTL"))
  allq=[float(v[0]["pvalue"]) for v in variants.values()]
  refn=int(full[gene]["n_variants"])
  if len(variants)>refn:raise ValueError("Region subset exceeds full-cis tested SNP count")
  region_n=len(variants);orig_n=len(archive)
  # Position-bounded denominator can be a different genome-wide universe.
  ss=dict(symbol=symbol,gene_id=gene,locus=locus,tissue="Brain_Cerebellar_Hemisphere",
    archived_common_snp_count=orig_n,catalogue_full_cis_tested_snp_count=refn,
    catalogue_unique_snp_in_archived_coordinate_region=region_n,
    catalogue_region_rsID_duplicate_rows=dup,
    archived_overlap_with_catalogue_region=len(matched),
    archived_source_only_within_region=len(source)-len(matched),
    catalogue_region_only_over_archived=len(variants)-len(matched),
    archived_SNP_count_div_full_catalogue_cis_count=orig_n/refn,
    shared_SNP_count_div_full_catalogue_cis_count=len(matched)/refn,
    catalogue_region_fraction_of_full_cis=region_n/refn,
    fraction_archived_exact_snps_seen_in_catalogue=len(matched)/orig_n,
    catalogue_permutation_best_p=float(full[gene]["pvalue"]),
    catalogue_permutation_p_beta=float(full[gene]["p_beta"]),
    catalogue_permutation_n_variants=refn,
    catalogue_region_min_nominal_p=min(allq) if allq else "",
    archive_region_min_nominal_p=min(float(r["eqtl_p"]) for r in archive),
    region_matched_assoc_p_within_1e3_fraction=(sum(abs(x-y)<.001 for x,y in same)/len(same) if same else ""),
    median_abs_neglog10_p_difference=(__import__("statistics").median(
      abs(math.log10(max(x,1e-300))-math.log10(max(y,1e-300))) for x,y in same) if same else ""),
    archived_shared_p_lt_1e4=sum(x<1e-4 for x,y in same),
    catalogue_shared_p_lt_1e4=sum(y<1e-4 for x,y in same),
    shared_snp_beta_sign_concordance=(
      sum(a["archived_beta"]*a["catalogue_beta"]>0 for a in differences if a["gene_id"]==gene)/
      sum(a["archived_beta"]*a["catalogue_beta"]!=0 for a in differences if a["gene_id"]==gene)
      if any(a["archived_beta"]*a["catalogue_beta"]!=0 for a in differences if a["gene_id"]==gene) else ""),
    n_cis_QTL_ancestry_matched_LD_validated=False,
    claim_gate="EXTERNAL_GTEX_REPROCESSED_GTEX_NOT_ORIGINAL_COLOC_QTL_INPUT")
  summaries.append(ss);audit+=differences
 out.mkdir(parents=True,exist_ok=True)
 for name,rows in [("IS_PRIORITY4_REPROCESSED_CIS_UNIVERSE.tsv",summaries),
                    ("IS_PRIORITY4_ARCHIVED_VS_REPROCESSED_SHARED_SNP.tsv",audit)]:
  with (out/name).open("w",newline="") as f:
   if not rows:continue
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
 report=dict(source_dataset="QTS000015/QTD000161_GTEx_reprocessed_Brain_Cerebellar_Hemisphere",
    reprocessed_summary_not_original_gtEx_nominal=True,
    full_cis_denominators_provenance="eQTL_Catalogue_permuted_n_variants_each_gene",
    exact_original_archive_coloc_uses_all_tested_variants_verified=False,
    full_study_GTEx_QTL_LD_available=False,no_causal_gene_claim=True,original_candidate_genes=2225,
    rows=summaries,
    SHA256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [permuted,*region_files]})
 (out/"IS_PRIORITY4_REPROCESSED_CIS_READINESS_MANIFEST.json").write_text(json.dumps(report,indent=2)+"\n")
 return report

if __name__=="__main__":
 p=argparse.ArgumentParser()
 p.add_argument("--permuted",type=Path,required=True)
 p.add_argument("--source",type=Path,required=True)
 p.add_argument("--cache",type=Path,required=True)
 p.add_argument("--fetch",action="store_true")
 p.add_argument("--output",type=Path,required=True)
 a=p.parse_args()
 a.cache.mkdir(parents=True,exist_ok=True)
 regions={"BBJ_IS_L001":"4:79759933-80762996","BBJ_IS_L002":"10:103079395-104199932"}
 region_files=[]
 for i,(locus,region) in enumerate(regions.items()):
  dest=a.cache/(locus+".QTD000161.region.priority4.tsv")
  if a.fetch:
   print("FETCH",locus,region,fetch_region(region,dest),flush=True)
   if i==0:time.sleep(3)
  region_files.append(dest)
 print(json.dumps(evaluate(a.permuted,a.source,region_files,a.output),indent=2))
