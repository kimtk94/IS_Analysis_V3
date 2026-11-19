#!/usr/bin/env python3
"""CI-only bounded row-level prerequisite audit. Never substitutes for primary-source verification."""
import csv,hashlib,json,sys
from pathlib import Path
from collections import Counter
K=sys.argv[1]
S={
"is":("gwas",["chr","pos","effect_allele","other_allele","beta","se","p"],"EUR/EAS full genome build and allele harmonization"),
"glaucoma":("dex",["donor","sample","condition","pool"],"paired raw integer count matrix and leave-pool-out"),
"hcc":("aliquot",["aliquot","case","channel"],"official PDC aliquot to TMT channel and RNA case mapping"),
"radiomics":("radiomics",["patient_id","image","rna"],"NIfTI header spacing and image-RNA-clinical patient mapping"),
"ocular":("ena",["run_accession","sample_accession","experiment_accession","fastq_md5"],"all 46 FASTQ content MD5 and strain/timepoint"),
"muscle":("muscle",["rsid"],"20 seed loci and EAS ancestry-matched LD"),
"oral":("oral",["donor","sample","endpoint"],"donor-expression-clinical endpoint one-to-one join"),
"skin":("skin",["snp","effect_allele","other_allele","beta","se","pval"],"complete GWAS effect harmonization and instrument validity")}
token,cols,gate=S[K]
aliases={"chr":["chromosome"],"pos":["position"],"rsid":["snp","variant_id"],"patient_id":["case_id"],"image":["image_path"],"rna":["rna_sample"],"donor":["donor_id"],"sample":["sample_id"],"condition":["treatment"],"pool":["batch"],"aliquot":["aliquot_id"],"case":["case_id","patient_id"],"channel":["tmt_channel"],"run_accession":["run"],"sample_accession":["sample"],"experiment_accession":["experiment"],"fastq_md5":["md5"],"endpoint":["outcome"],"snp":["rsid"],"effect_allele":["ea"],"other_allele":["oa"],"pval":["p_value","p"]}
out={"date":"2026-11-20","research":K,"technical_status":"HOLD","scientific_status":"HOLD","gate":gate,"files":[],"reason":"NO_VERIFIED_PRIMARY_DATA"}
paths=[]
for p in Path(".").rglob("*"):
 if not p.is_file() or p.suffix.lower() not in (".tsv",".csv"):continue
 if any(z in p.parts for z in (".git","node_modules","data","raw","fastq","results")):continue
 if token not in str(p).lower() and not (K=="muscle" and p.name=="stage1_seed_loci.tsv"):continue
 if 0<p.stat().st_size<=2000000:paths.append(p)
for p in sorted(paths)[:100]:
 try:
  with p.open(encoding="utf-8-sig",errors="replace",newline="") as f:
   reader=csv.DictReader(f,delimiter="\t" if p.suffix.lower()==".tsv" else ",")
   names={x.lower().strip():x for x in (reader.fieldnames or [])}
   match={c:next((names[x] for x in [c]+aliases.get(c,[]) if x in names),None) for c in cols}
   if any(v is None for v in match.values()):continue
   n=0; missing=0; duplicate=0; seen=set(); bad_alleles=0
   for row in reader:
    n+=1
    if any(not str(row.get(v,"")).strip() for v in match.values()):missing+=1
    if K in ("hcc","ocular","oral"):
     key=tuple(row.get(match[c],"") for c in (["aliquot","channel"] if K=="hcc" else ["run_accession"] if K=="ocular" else ["donor","sample"]))
     if key in seen:duplicate+=1
     seen.add(key)
    if K in ("is","skin"):
     for c in ("effect_allele","other_allele"):
      if c in match and str(row.get(match[c],"")).upper() not in ("A","C","G","T"):bad_alleles+=1
   out["files"].append({"path":str(p),"sha256":hashlib.sha256(p.read_bytes()).hexdigest(),"rows":n,"missing_required_rows":missing,"duplicate_keys":duplicate,"non_acgt_alleles":bad_alleles})
 except (OSError,UnicodeError,csv.Error) as e:
  out.setdefault("read_errors",[]).append({"path":str(p),"error":type(e).__name__})
if out["files"]:
 out["technical_status"]="PASS"
 out["reason"]="BOUNDED_TRACKED_FILE_AUDIT_ONLY_NOT_PRIMARY_SOURCE_OR_STATISTICAL_PASS"
if K=="muscle" and out["files"]:
 out["seed_loci_count_check"]=[{"path":f["path"],"observed":f["rows"],"equals_20":f["rows"]==20} for f in out["files"]]
Path("results").mkdir(exist_ok=True)
Path("results",K+".json").write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
print(json.dumps(out,ensure_ascii=False))
