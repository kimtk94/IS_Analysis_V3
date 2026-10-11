#!/usr/bin/env python3
"""Independent all-selected-SNP raw↔canonical GIGASTROKE EAS AIS provenance.

Checks every FASTA-anchored EAS AIS SNP used in 98 public eGene assays:
original EBI MD5, source raw vs derived canonical OA/EA/beta/SE/P/EAF,
GRCh37 FASTA REF assignment and historical ID other:effect convention.
Read-only scientific derivative. Not independent GWAS replication or coloc.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,math,sys
from collections import Counter
from pathlib import Path
MD5="904e12e33ec3c2eeefe4dc93bec4995c"
def md5(p):
 h=hashlib.md5()
 with p.open("rb") as f:
  for buf in iter(lambda:f.read(4*1024*1024),b""):h.update(buf)
 return h.hexdigest()
def expected(gz):
 variants=set()
 positions=set()
 rows=0
 with gzip.open(gz,"rt") as f:
  for r in csv.DictReader(f,delimiter="\t"):
   rows+=1
   if r["EAS_AIS_allele_exact_Fasta37"]!="True":continue
   key=r["mapped_genomic_variant_GRCh37"].upper()
   bits=key.split(":")
   if len(bits)!=4 or len(bits[2])!=1 or len(bits[3])!=1:
    raise ValueError("Unexpected non-SNV GIGASTROKE EAS normalized variant")
   variants.add(key)
   positions.add((bits[0],int(bits[1])))
 if rows<50000 or not variants:raise ValueError("Missing published full-cis matched record data")
 return variants,positions,rows

def pick(path,positions,is_raw):
 data={};stats=Counter()
 chromosomes={c for c,_ in positions}
 with gzip.open(path,"rt") as f:
  hdr=f.readline().rstrip("\r\n").split("\t")
  ix={c:i for i,c in enumerate(hdr)}
  ck="chromosome" if is_raw else "chr"
  pk="base_pair_location" if is_raw else "pos"
  cols=("effect_allele","other_allele","beta","standard_error" if is_raw else "se",
    "p_value" if is_raw else "p","effect_allele_frequency" if is_raw else "eaf")
  if not set((ck,pk,*cols)).issubset(ix):raise ValueError("GWAS source header drift")
  for line in f:
   if not line.strip():continue
   stats["total_source_records"]+=1
   c=line.rstrip("\r\n").split("\t")
   chrom=c[ix[ck]].removeprefix("chr")
   if chrom not in chromosomes:continue
   try:pos=int(c[ix[pk]])
   except ValueError:stats["invalid_position"]+=1;continue
   key=(chrom,pos)
   if key not in positions:continue
   vals=[c[ix[z]] for z in cols]
   ea,oa,b,se,pv,af=vals
   snpkey=(chrom,pos,ea.upper(),oa.upper())
   if snpkey in data:raise ValueError("Duplicated raw/canonical genomic-allele GWAS signature")
   canon=dict(chrom=chrom,pos=pos,effect_allele=ea.upper(),
     other_allele=oa.upper(),beta=b,se=se,p=pv,eaf=af)
   if not is_raw:
    canon["historical_variant_id"]=c[ix["variant_id"]]
    canon["build"]=c[ix["build"]]
    canon["dataset"]=c[ix["dataset"]]
    canon["ancestry"]=c[ix["ancestry"]]
   data[snpkey]=canon
   stats["rows_at_target_position"]+=1
 if stats["total_source_records"]<5000000:raise ValueError("Source file unexpectedly short")
 return data,stats
def eq_float(a,b):
 try:
  x=float(a);y=float(b)
  return math.isfinite(x) and math.isfinite(y) and math.isclose(x,y,rel_tol=1e-11,abs_tol=1e-12)
 except (TypeError,ValueError):return False
def audit(raw,can,targets,fa,normalize):
 if set(raw)!=set(can):
  missingraw=set(can)-set(raw);missingcan=set(raw)-set(can)
  raise ValueError(f"Original/derived source record set differs: rawmissing={len(missingraw)} canmissing={len(missingcan)}")
 qc=Counter();records=[]
 for key in sorted(raw):
  r=raw[key];c=can[key]
  for x in ("beta","se","p","eaf"):
   if not eq_float(r[x],c[x]):raise ValueError("GWAS original and derived numeric effect mismatch "+str(key)+" "+x)
  if c["build"]!="GRCh37" or c["dataset"]!="GCST90104545" or c["ancestry"]!="EAS":
   raise ValueError("Derived EAS AIS GWAS source cohort/build drift")
  chrom,pos,ea,oa=key
  # normalize_record expects canonical source with historical ID string,
  # *not* an invented biological REF:ALT tuple.
  source={"chr":chrom,"pos":str(pos),"build":"GRCh37",
   "effect_allele":ea,"other_allele":oa,"variant_id":c["historical_variant_id"],
   "beta":c["beta"],"eaf":c["eaf"]}
  ref=fa.get(chrom,pos)
  item=normalize(source,ref)
  if item.variant_id not in targets:
   qc["GWAS_OTHER_ALLELE_AT_SHARED_POSITION_BUT_NOT_MATCHED_EGENE_SNP"]+=1
   continue
  if item.historical_id_status not in ("FASTA_REF_ALT","OTHER_EFFECT_ORDER"):
   raise ValueError("Original GWAS variant ID order not independently verified")
  qc["RAW_CANONICAL_ALLELE_EFFECT_FREQ_NUMERIC_AND_FASTA_VERIFIED"]+=1
  qc[item.historical_id_status]+=1
  records.append({
   "source":"GIGASTROKE_EAS_AIS_GCST90104545",
   "FASTA_hg19_refalt_variant_id":item.variant_id,
   "historic_preGWAS_SSF_derived_variant_id":item.historical_id,
   "historical_id_order":item.historical_id_status,
   "beta_effect_ALT":item.beta_alt,
   "EAF_ALT":item.eaf_alt,
   "raw_vs_canonical_beta_SE_P_EAF_equal":True,
   "ref_allele_independent_GRCh37_FASTA":ref,
   "science_gate":"ASSOCIATION_INPUT_PROVENANCE_NOT_COLOCALIZATION"
  })
 if not records:raise ValueError("No raw to canonical matched EAS AIS GWAS variants")
 if len({r["FASTA_hg19_refalt_variant_id"] for r in records})!=len(records):
  raise ValueError("Two distinct source GWAS variants recoded to same FASTA REF_ALT ID")
 return records,dict(qc)
def main():
 ap=argparse.ArgumentParser()
 for name in ("match-gz","raw","canonical","fasta","out"):
  ap.add_argument("--"+name,required=True,type=Path)
 a=ap.parse_args()
 if a.out.exists() and any(a.out.iterdir()):raise FileExistsError("Refuse audit output overwrite")
 if md5(a.raw).lower()!=MD5:raise ValueError("Original GIGASTROKE EAS AIS MD5 failed")
 target,pos,source_records=expected(a.match_gz)
 r,rs=pick(a.raw,pos,True)
 c,cs=pick(a.canonical,pos,False)
 sys.path.insert(0,str(Path(__file__).resolve().parent))
 from audit_is_eas5_public_qtl_gwas_snp_overlap import Fasta37
 from gigastroke_fasta_variant_contract import normalize_record
 fa=Fasta37(a.fasta)
 try:rows,counts=audit(r,c,target,fa,normalize_record)
 finally:fa.close()
 a.out.mkdir(parents=True)
 with (a.out/"IS_EAS5_98_EGENE_GIGA_FASTA_ANCHORED_RAW_CANONICAL_SNP_PROVENANCE.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
  w.writeheader();w.writerows(rows)
 summary={
  "status":"SAMPLED_EGENE_SNP_UNIVERSE_ALL_CONTRIBUTING_GIGASTROKE_SOURCE_ROWS_RAW_CANONICAL_FASTA_QC_PASS",
  "all_original_study_records_scanned":rs["total_source_records"],
  "all_canonical_study_records_scanned":cs["total_source_records"],
  "unique_input_target_Fasta37_SNPs":len(target),
  "source_geneQTL_variant_rows_scanned":source_records,
  "GIGA_original_EBI_RAW_MD5":MD5,
  "verified_records":len(rows),
  "historical_allele_order_counts":counts,
  "no_independent_stroke_repl_or_new_coloc":True
 }
 (a.out/"IS_EAS5_98_GIGA_RAW_CANONICAL_FASTA_QC_MANIFEST.json").write_text(json.dumps(summary,indent=2)+"\n")
 print("IS_EAS5_98_GIGA_RAW_CANONICAL_FASTA_REFALT_ALL_MATCHED_SNP_PASS",json.dumps({
  k:summary[k] for k in ("all_original_study_records_scanned","unique_input_target_Fasta37_SNPs","verified_records","historical_allele_order_counts")}),flush=True)
if __name__=="__main__":main()
