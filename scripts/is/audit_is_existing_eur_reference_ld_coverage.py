#!/usr/bin/env python3
"""Read-only candidate EUR reference region SNP overlap for 30 archived IS assays."""
import csv,json,hashlib,argparse,gzip
from pathlib import Path
from collections import defaultdict

def rows(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def read_pvar(p):
 opener=gzip.open if str(p).endswith(".gz") else open
 with opener(p,"rt") as f:
  for line in f:
   if line.startswith("##"):continue
   if line.startswith("#"):columns=line.lstrip("#").rstrip("\n").split("\t");continue
   x=dict(zip(columns,line.rstrip("\n").split("\t")))
   yield x
def run(assays,source,eur,out):
 q=[x for x in rows(assays) if x["status"]=="PASS"]
 if len(q)!=30:raise ValueError("Require 30 source assays")
 panels=[]
 for f in sorted(eur.rglob("*.pgen")):
  base=f.with_suffix("")
  pvars=[Path(str(base)+".pvar"),Path(str(base)+".pvar.zst"),Path(str(base)+".pvar.gz")]
  pv=next((p for p in pvars if p.exists() and p.suffix!=".zst"),None)
  if pv is None:continue
  psam=Path(str(base)+".psam")
  ids=set();sites=set();n=0
  for x in read_pvar(pv):
   n+=1
   if x.get("ID"):ids.add(x["ID"])
   if all(z in x for z in ("#CHROM","POS","REF","ALT")):
    sites.add((x["#CHROM"].replace("chr",""),x["POS"],x["REF"],x["ALT"]))
  panels.append((str(f),str(pv),len(ids),sites,ids,n,psam.exists()))
 results=[]
 for a in q:
  file=source/a["filename"]
  if not file.exists():raise ValueError("Missing original source "+str(file))
  raw=rows(file)
  target={(x["chromosome"],x["position"],x["ref"],x["alt"]) for x in raw}
  hits=[]
  for f,pv,nid,sites,ids,n,haspsam in panels:
   overlap=len(target & sites)
   if overlap:
    hits.append(dict(panel=f,matched_exact_grch38_sites=overlap,panel_n_variants=n,has_psam=haspsam))
  hits.sort(key=lambda x:-x["matched_exact_grch38_sites"])
  results.append(dict(gene_id=a["gene_base"],locus=a["locus"],dataset_key=a["dataset_key"],
   source_SNPs=len(raw),best_available_EUR_pgen_exact_GRCh38_matching_SNPs=(hits[0]["matched_exact_grch38_sites"] if hits else 0),
   best_EUR_panel=(hits[0]["panel"] if hits else ""),
   best_panel_n_variants=(hits[0]["panel_n_variants"] if hits else 0),
   independent_GTEx_donor_LD_available=False,
   usable_as_GTex_in_study_LD=False,
   confidence="EXPLORATORY_POSITION_REF_ALT_COVERAGE_ONLY"))
 out.mkdir(parents=True,exist_ok=True)
 with (out/"IS_30_GTEX_EUR_EXISTING_REFERENCE_PANEL_COVERAGE.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter="\t");w.writeheader();w.writerows(results)
 j=dict(audited_archived_assays=len(results),available_EUR_PGEN_with_readable_PVAR=len(panels),
  assays_with_any_EUR_site_overlap=sum(x["best_available_EUR_pgen_exact_GRCh38_matching_SNPs"]>0 for x in results),
  assays_with_90pct_coverage=sum(x["best_available_EUR_pgen_exact_GRCh38_matching_SNPs"]>=.9*x["source_SNPs"] for x in results),
  verified_actual_GTex_cohort_LD_panels=0,
  caveat="External 1000G EUR reference is not GTEx study cohort LD; some panels may be GRCh37 and raw comparison requires build QC. Empty overlap does not imply unavailable EUR genotype elsewhere.",
  assays_sha256=hashlib.sha256(assays.read_bytes()).hexdigest())
 (out/"IS_30_EUR_REFERENCE_SEARCH_SUMMARY.json").write_text(json.dumps(j,indent=2)+"\n")
 return j
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for key in ("assays","source","eur","out"):p.add_argument("--"+key,type=Path,required=True)
 a=p.parse_args()
 print(json.dumps(run(a.assays,a.source,a.eur,a.out),indent=2))
