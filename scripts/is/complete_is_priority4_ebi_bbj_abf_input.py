#!/usr/bin/env python3
"""Join unmodified canonical BBJ GWAS standard errors into four-gene QTD pilot.
Never derive beta standard errors from rounded GWAS P values.
"""
import argparse,csv,gzip,hashlib,json,math
from collections import defaultdict
from pathlib import Path

def read(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def execute(pre,bbj,out):
 if out.exists() and any(out.iterdir()):raise FileExistsError("Do not overwrite canonical pilot")
 r=read(pre)
 if len(r)!=11050:raise ValueError("Candidate SNP denominator drift")
 targets={x["variant_GRCh37_REF_ALT"] for x in r}
 match={}
 positions=defaultdict(set)
 for key in targets:
  a=key.split(":");positions[a[0]].add(int(a[1]))
 n_total=0;n_positions=0;n_valid=0
 with gzip.open(bbj,"rt") as f:
  line=f.readline().rstrip("\r\n").split("\t");ix={c:i for i,c in enumerate(line)}
  required={"chr","pos","ref","alt","effect_allele","other_allele","beta","se","p","eaf","variant_id"}
  if not required.issubset(ix):raise ValueError("Canonical BBJ schema drift")
  for line in f:
   n_total+=1
   s=line.rstrip("\r\n").split("\t");chrom=s[ix["chr"]]
   if chrom not in positions:continue
   try:pos=int(s[ix["pos"]])
   except ValueError:continue
   if pos not in positions[chrom]:continue
   n_positions+=1
   key=s[ix["variant_id"]]
   if key not in targets:continue
   ref=s[ix["ref"]];alt=s[ix["alt"]]
   if key!=f"{chrom}:{pos}:{ref}:{alt}":
    raise ValueError("BBJ canonical REF/ALT contract violated")
   ea=s[ix["effect_allele"]];oa=s[ix["other_allele"]]
   if (ea,oa) not in ((ref,alt),(alt,ref)):
    raise ValueError("BBJ effect allele differs from canonical REF/ALT")
   beta=float(s[ix["beta"]]);se=float(s[ix["se"]]);p=float(s[ix["p"]]);eaf=float(s[ix["eaf"]])
   if not all(math.isfinite(v) for v in (beta,se,p,eaf)) or se<=0 or not 0<=p<=1 or not 0<eaf<1:
    raise ValueError("Invalid GWAS numeric source")
   beta_alt=beta if ea==alt else -beta
   eaf_alt=eaf if ea==alt else 1-eaf
   if key in match:raise ValueError("Multiple BBJ canonical GWAS variants with same REF/ALT")
   match[key]=(beta_alt,se,p,eaf_alt)
   n_valid+=1
 if n_total<13000000:raise ValueError("Full BBJ source truncated")
 if match.keys()!=targets:raise ValueError(f"Missing BBJ exactly matched variants {len(targets-set(match))}")
 for x in r:
  v=match[x["variant_GRCh37_REF_ALT"]]
  if not math.isclose(float(x["BBJ_beta_ALT"]),v[0],abs_tol=1e-10,rel_tol=1e-9):
   raise ValueError("BBJ beta ALT changed since FASTA-verified join")
  if not math.isclose(float(x["BBJ_p"]),v[2],abs_tol=1e-10,rel_tol=1e-9):
   raise ValueError("BBJ p value changed since FASTA-verified join")
  if not math.isclose(float(x["BBJ_EAF_ALT"]),v[3],abs_tol=1e-10,rel_tol=1e-9):
   raise ValueError("BBJ EAF ALT changed since FASTA-verified join")
  x["BBJ_se"]=v[1]
 out.mkdir(parents=True)
 dest=out/"IS_PRIORITY4_EBI_BBJ_ABF_QC_VALIDATED_INPUT.tsv"
 with dest.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(r[0]),delimiter="\t");w.writeheader();w.writerows(r)
 stats=dict(n_original_ebi_full_cis_geneSNP_records=20689,
   n_qtl_bbj_allele_exact_nonpal_geneSNP_records=len(r),
   unique_BBJ_exact_REF_ALT=len(targets),
   BBJ_full_source_rows_scanned=n_total,BBJ_target_positions_seen=n_positions,
   BBJ_unique_target_snp_found=n_valid,
   BBJ_N_assumed=174686,BBJ_case_fraction_assumed=22664/174686,
   QTL_EFFECTIVE_N_NOT_ATTESTED=True,
   QTL_genotype_AN_div_2_ONLY_LIMIT_UPPER_SAMPLE_PROXY=True,
   status="DESCRIPTIVE_SENSITIVITY_ONLY_NO_COHORT_LD_OR_ORIGINAL_GTEX_CONFIRMED")
 (out/"IS_PRIORITY4_EBI_BBJ_ABF_PILOT_INPUT_MANIFEST.json").write_text(json.dumps(stats,indent=2)+"\n")
 return stats
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for x in ("pre","bbj","out"):p.add_argument("--"+x,type=Path,required=True)
 a=p.parse_args()
 print(json.dumps(execute(a.pre,a.bbj,a.out),indent=2))
