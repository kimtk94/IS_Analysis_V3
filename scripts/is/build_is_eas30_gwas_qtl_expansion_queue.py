#!/usr/bin/env python3
"""Prioritize EAS 30 IS regions by observed cross-outcome GWAS and QTL gaps.

A cross-outcome regional P association is NOT replication or gene causality;
GWS regions can overlap biologically and involve overlapping study cohorts.
Output is a workflow queue only, preserving all 2225 positional genes.
"""
import argparse,csv,json,hashlib
from collections import Counter,defaultdict
from pathlib import Path
def read(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def hashfile(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for buf in iter(lambda:f.read(1048576),b""):h.update(buf)
 return h.hexdigest()
def build(gwas,regions,genes,pilot=None):
 if len(gwas)!=180 or len(regions)!=30 or len(genes)!=2225:
  raise ValueError("Expected 180 six-outcome GWAS + 30 regional + 2225 genes")
 if any(len(x["region_id"])<10 for x in gwas):raise ValueError("Region missing")
 g=defaultdict(dict)
 for r in gwas:
  rid=r["region_id"];trait=r["source_dataset"]
  if trait in g[rid]:raise ValueError("Duplicate region outcome")
  g[rid][trait]=r
 pilotmap=defaultdict(list)
 if pilot is not None:
  if len(pilot) not in (6,15):raise ValueError("Requires 6 or 15 verified public nominal QTL slices")
  for r in pilot:
   if r["status"]!="PUBLIC_NOMINAL_QTL_REGION_DOWNLOADED":
    continue
   pilotmap[r["region_id"]].append(r)
 output=[]
 for row in regions:
  rid=row["EAS_region_id"]
  if set(g[rid])!={"BBJ_IS","GWS_EAS_AS","GWS_EAS_AIS","GWS_EAS_CES","GWS_EAS_LAS","GWS_EAS_SVS"}:
   raise ValueError("Incomplete six outcome source coverage")
  p=lambda trait:float(g[rid][trait]["best_p"]) if g[rid][trait]["best_p"] else None
  bbj,ais,as_,ces,las,svs=[p(t) for t in (
    "BBJ_IS","GWS_EAS_AIS","GWS_EAS_AS","GWS_EAS_CES","GWS_EAS_LAS","GWS_EAS_SVS")]
  if bbj is None or ais is None:raise ValueError("Missing EAS AIS / BBJ source coverage")
  if bbj<5e-8 and ais<5e-8:
   tier="T1_BBJ_AND_EAS_AIS_REGION_ASSOCIATED_NOT_INDEPENDENT_REPLICATION"
  elif ces is not None and ces<5e-8:
   tier="T2_CES_SUBTYPE_ASSOCIATED_EAS_AIS_NOT_GWS"
  elif as_ is not None and as_<5e-8:
   tier="T3_AS_ASSOCIATED_EAS_AIS_NOT_GWS"
  else:
   tier="T4_SUGGESTIVE_REGION_NO_5E8_SIGNAL_IN_PRIMARY_OUTCOMES"
  events=pilotmap[rid]
  output.append({
   "region_id":rid,"primary_GWAS_EAS_BBJ_AIS_tier":tier,
   "region_is_provisional_not_independent":True,
   "n_positional_gene_ids":int(row["n_positional_genes"]),
   "n_archived_QTL_tested_gene_ids":int(row["n_genes_archived_GTEx_ABF_tested"]),
   "n_EAS_positional_genes_without_archived_coloc":int(row["n_genes_not_tested"]),
   "BBJ_500kb_min_p":bbj,"GWS_EAS_AIS_500kb_min_p":ais,
   "GWS_EAS_AS_500kb_min_p":as_,"GWS_EAS_CES_500kb_min_p":ces,
   "GWS_EAS_LAS_500kb_min_p":las,"GWS_EAS_SVS_500kb_min_p":svs,
   "BBJ_lead_exact_source_present":g[rid]["BBJ_IS"]["group_source_lead_present"],
   "GWS_EAS_AIS_lead_exact_source_present":g[rid]["GWS_EAS_AIS"]["group_source_lead_present"],
   "public_QTL_pilot_datasets_downloaded":len(events),
   "public_QTL_pilot_datasets": ";".join(sorted(r["dataset_id"] for r in events)),
   "new_public_QTL_nominal_row_coverage_not_coloc":bool(events),
   "new_cohort_matched_QTL_LD_verified":False,
   "cross_source_genetic_replication_tested":False,
   "causal_gene_claim":"NOT_ESTABLISHED",
   "triage_reason":"GWAS_TIER_AND_UNTESTED_QTL_COVERAGE_ONLY"
  })
 if len(output)!=30 or sum(x["n_positional_gene_ids"] for x in output)!=869 or sum(x["n_EAS_positional_genes_without_archived_coloc"] for x in output)!=826 or sum(x["n_archived_QTL_tested_gene_ids"] for x in output)!=43:
  raise ValueError("EAS source/gene coverage counts changed")
 tiers=Counter(x["primary_GWAS_EAS_BBJ_AIS_tier"][:2] for x in output)
 if dict(tiers)!={"T1":4,"T2":1,"T3":2,"T4":23}:
  raise ValueError("Observed source-derived region support tier count changed")
 if sum(x["n_EAS_positional_genes_without_archived_coloc"] for x in output if x["primary_GWAS_EAS_BBJ_AIS_tier"].startswith("T1"))!=164:
  raise ValueError("164 untested genes in four primary-associated EAS+BBJ regions changed")
 output.sort(key=lambda x:(x["primary_GWAS_EAS_BBJ_AIS_tier"],-x["n_EAS_positional_genes_without_archived_coloc"],x["region_id"]))
 return output,tiers
def main():
 p=argparse.ArgumentParser()
 for k in ("gwas-context","region-gap","gene-universe","out-dir"):
  p.add_argument("--"+k,required=True,type=Path)
 p.add_argument("--public-QTL-pilot",type=Path)
 a=p.parse_args()
 ins={"GWS_REGION":a.gwas_context,"REGION_GAP":a.region_gap,"GENES":a.gene_universe}
 if a.public_QTL_pilot:ins["QTL_PILOT"]=a.public_QTL_pilot
 for x in ins.values():
  if not x.is_file():raise FileNotFoundError(x)
 if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError("No overwrite")
 data,tiers=build(read(a.gwas_context),read(a.region_gap),read(a.gene_universe),
  read(a.public_QTL_pilot) if a.public_QTL_pilot else None)
 a.out_dir.mkdir(parents=True,exist_ok=True)
 with (a.out_dir/"IS_EAS30_OUTCOME_QTL_EXPANSION_PRIORITY.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter="\t")
  w.writeheader();w.writerows(data)
 manifest={
  "status":"DESCRIPTIVE_SIX_OUTCOME_REGION_ASSOCIATION_AND_QTL_COVERAGE_QUEUE",
  "all_original_EAS_regions":30,"all_original_positional_genes_preserved":2225,
  "n_EAS_region_genes":869,"EAS_QTL_coloc_historically_tested_genes":43,
  "EAS_no_archived_coloc_genes":826,
  "GWAS_tier_counts":dict(tiers),
  "T1_regions":[x["region_id"] for x in data if x["primary_GWAS_EAS_BBJ_AIS_tier"].startswith("T1")],
  "T1_untested_QTL_genes":164,
  "GWAS_sources_are_not_independent_replications":True,
  "nominal_QTL_pilot_is_not_coloc":True,
  "input_SHA256":{k:hashfile(v) for k,v in ins.items()}
 }
 (a.out_dir/"IS_EAS30_GWAS_QTL_RESEARCH_QUEUE_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
 print("IS_EAS30_OUTCOME_QTL_PRIORITIES_PASS",json.dumps(manifest),flush=True)
if __name__=="__main__":main()
