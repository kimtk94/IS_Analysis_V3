#!/usr/bin/env python3
"""CI-only scientific prerequisite gates. No inferred patient joins or invented source evidence."""
import csv, hashlib, json, os, sys
from pathlib import Path
from collections import Counter
key=sys.argv[1]
spec={
"is":("GWAS SNP harmonization prerequisites",["chr","pos","effect_allele","other_allele","beta","se","p"],"gwas"),
"glaucoma":("DEX donor-paired design and count provenance",["donor","sample","condition","pool"],"dex"),
"hcc":("PDC aliquot-case-TMT one-to-one consistency",["aliquot","case","channel"],"aliquot"),
"radiomics":("image-RNA-patient link uniqueness",["patient_id","image","rna"],"radiomics"),
"ocular":("ENA 46-run identity and MD5 declarations",["run_accession","sample_accession","experiment_accession","fastq_md5"],"ena"),
"muscle":("20 seed loci distinctness and coordinates",["rsid"],"muscle"),
"oral":("oral donor-sample-endpoint consistency",["donor","sample","endpoint"],"oral"),
"skin":("MR summary-statistics allele and numeric sanity",["snp","effect_allele","other_allele","beta","se","pval"],"skin")
}
aliases={"chr":["chromosome"],"pos":["position"],"effect_allele":["ea"],"other_allele":["oa"],"p":["pval","p_value"],"pval":["p","p_value"],"donor":["donor_id"],"sample":["sample_id"],"pool":["batch"],"condition":["treatment"],"aliquot":["aliquot_id"],"case":["case_id"],"channel":["tmt_channel"],"patient_id":["case_id"],"image":["image_path"],"rna":["rna_sample"],"run_accession":["run"],"sample_accession":["sample"],"experiment_accession":["experiment"],"fastq_md5":["md5"],"rsid":["snp","variant_id"],"snp":["rsid"],"endpoint":["outcome"]}
label,required,token=spec[key]
out={"date":"2026-11-21","research":key,"task":label,"technical_status":"HOLD","scientific_status":"HOLD","source_scope":"tracked small files only; not primary source","evidence":[],"blockers":[]}
candidates=[]
for root,dirs,files in os.walk("."):
 dirs[:]=[d for d in dirs if d not in (".git","node_modules","raw","data","results",".venv")]
 for name in files:
  p=Path(root)/name
  if p.suffix.lower() not in (".csv",".tsv") or token not in str(p).lower() and not (key=="muscle" and name=="stage1_seed_loci.tsv"):continue
  if 0<p.stat().st_size<=2_000_000:candidates.append(p)
for p in sorted(candidates)[:50]:
 try:
  with p.open(encoding="utf-8-sig",newline="") as f:
   rd=csv.DictReader(f,delimiter="\t" if p.suffix.lower()==".tsv" else ",")
   names={n.strip().lower():n for n in (rd.fieldnames or [])}
   cols={c:next((names[n] for n in [c]+aliases.get(c,[]) if n in names),None) for c in required}
   if any(v is None for v in cols.values()):continue
   rows=list(rd)
  missing=sum(any(not str(r.get(v) or "").strip() for v in cols.values()) for r in rows)
  dup=0; invalid=0; seen=set(); identities={}
  for r in rows:
   if key=="hcc": identity=(r[cols["aliquot"]],r[cols["channel"]]); target=r[cols["case"]]
   elif key=="ocular": identity=(r[cols["run_accession"]],); target=(r[cols["sample_accession"]],r[cols["experiment_accession"]])
   elif key=="oral": identity=(r[cols["donor"]],r[cols["sample"]]);target=r[cols["endpoint"]]
   elif key=="radiomics":identity=(r[cols["patient_id"]],r[cols["image"]]);target=r[cols["rna"]]
   elif key=="glaucoma":identity=(r[cols["sample"]],);target=(r[cols["donor"]],r[cols["condition"]],r[cols["pool"]])
   elif key=="muscle":identity=(r[cols["rsid"]],);target=None
   elif key=="is":identity=(r[cols["chr"]],r[cols["pos"]],r[cols["effect_allele"]],r[cols["other_allele"]]);target=None
   else:identity=(r[cols["snp"]],);target=None
   if identity in seen:dup+=1
   seen.add(identity)
   if identity in identities and identities[identity]!=target:invalid+=1
   identities[identity]=target
   if key in ("is","skin"):
    ea=str(r[cols["effect_allele"]]).upper();oa=str(r[cols["other_allele"]]).upper()
    if ea not in "ACGT" or oa not in "ACGT" or ea==oa:invalid+=1
    try:
     b=float(r[cols["beta"]]);s=float(r[cols["se"]]);pv=float(r[cols["pval" if key=="skin" else "p"]])
     if not all(map(__import__("math").isfinite,(b,s,pv))) or s<=0 or not 0<=pv<=1:invalid+=1
    except (ValueError,TypeError):invalid+=1
   if key=="ocular" and len(str(r[cols["fastq_md5"]]).strip()) not in (32,65):invalid+=1
  ev={"path":str(p),"sha256":hashlib.sha256(p.read_bytes()).hexdigest(),"rows":len(rows),"missing":missing,"duplicate_keys":dup,"invalid_or_conflicting":invalid}
  if key=="ocular":ev["equals_46_runs"]=len(rows)==46
  if key=="muscle":ev["equals_20_loci"]=len(rows)==20
  out["evidence"].append(ev)
 except (OSError,UnicodeError,csv.Error) as e:out["blockers"].append("INPUT_READ_ERROR:"+type(e).__name__)
if not out["evidence"]:out["blockers"].append("NO_MATCHING_VERIFIED_PRIMARY_INPUT")
else:
 if any(e["missing"] or e["duplicate_keys"] or e["invalid_or_conflicting"] for e in out["evidence"]):out["blockers"].append("TRACKED_INPUT_QC_ANOMALIES")
 if key=="ocular" and not any(e["equals_46_runs"] for e in out["evidence"]):out["blockers"].append("RUN_COUNT_NOT_46")
 if key=="muscle" and not any(e["equals_20_loci"] for e in out["evidence"]):out["blockers"].append("SEED_COUNT_NOT_20")
 if not out["blockers"]:out["technical_status"]="PASS"
out["blockers"].append("PRIMARY_SOURCE_CHECKSUM_AND_INDEPENDENT_SCIENTIFIC_VALIDATION_NOT_ESTABLISHED")
out["next_ci_action"]={
"is":"full EUR/EAS source harmonization with genome build and ancestry LD",
"glaucoma":"obtain donor-linked integer count matrix and leave-pool-out CI",
"hcc":"official PDC aliquot-TMT-to-TCGA case readback",
"radiomics":"IBSI/NIfTI spacing and exact patient ID MRI-RNA-clinical join",
"ocular":"download 46 FASTQ in CI and independently verify full content MD5",
"muscle":"ancestry-matched EAS LD plus GTEx eQTL coloc",
"oral":"donor-level expression and clinical endpoint crosswalk",
"skin":"instrument clumping, harmonization and pleiotropy sensitivity"
}[key]
Path("results").mkdir(exist_ok=True)
Path("results",key+".json").write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
print(json.dumps({"research":key,"technical":out["technical_status"],"scientific":out["scientific_status"],"evidence_files":len(out["evidence"]),"blockers":out["blockers"]}))
