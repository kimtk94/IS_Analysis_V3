#!/usr/bin/env python3
"""Build 2,225-gene / 2,425 gene-region IS evidence overlays from existing files.

No new databases, no source modifications, no inferential evidence promotions.

Principles:
- Stable Ensembl ID is the gene key; never merge by alias alone.
- The 1 legacy-only NEURL1 anchor is kept separate, not inserted as 2,226th
  unique positional gene.
- Historic coloc is linked only to EAS/Japanese gene-region rows that explicitly
  match the BBJ locus AND gene stable ID. No legacy coloc leakage to EUR.
- GTEx coloc H4 is a descriptive archived posterior, *not* a causal gene score.
- Missing QTL/feature status is NOT_TESTED/NOT_ASSESSABLE, never negative.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter,defaultdict
from pathlib import Path

REL={
 "expanded":"stage5_functional/broad_discovery_v2_ancestry/IS_ANCESTRY_GENE_EVIDENCE_V3.tsv",
 "windows":"stage5_functional/broad_discovery_v2_ancestry/IS_ALL_GENE_WINDOW_UNIVERSE.tsv",
 "universe":"stage5_functional/broad_discovery_v2_ancestry/IS_ALL_GENE_UNIVERSE_SUMMARY.json",
 "coloc":"stage5_functional/phase9c_convergence/COLOC_ABF_MASTER_ANNOTATED_V2.tsv",
 "replay":"audits/legacy_coloc_646_batch_cache30_20261010_v1/IS_LEGACY_646_SNP_REPLAY_STATUS.tsv",
 "priority8":"audits/is_priority8_source_verified_20261010_v2/IS_PRIORITY8_SOURCE_VERIFIED.tsv",
 "literature":"stage5_functional/phase9d_literature_benchmark/GENE_EVIDENCE_LITERATURE_MASTER.tsv",
 "cell_gene":"stage5_functional/phase9f_e_human_r3_7_donor_import_20261010/HUMAN_TARGET_GENE_SUMMARY.tsv",
 "cell_feature":"stage5_functional/phase9f_e_human_r3_7_donor_import_20261010/HUMAN_DONOR_FEATURE_STATUS.tsv",
 "cell_pairs":"stage5_functional/phase9f_e_human_r3_7_donor_import_20261010/HUMAN_DONOR_PAIRED_COMPARISONS.tsv",
}
def read(path):
 with path.open(newline="",encoding="utf-8") as f:
  return list(csv.DictReader(f,delimiter="\t"))
def digest(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b""):h.update(chunk)
 return h.hexdigest()
def fnum(x):
 try: return float(x)
 except (ValueError,TypeError): return None
def tswrite(path,rows,cols):
 with path.open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=cols,delimiter="\t",extrasaction="raise")
  w.writeheader();w.writerows(rows)
def exact_row_key(row):
 return row["group_id"],row["gene_id_stable"]
def stable(val):
 return val.split(".")[0]

def build(sources):
 expanded=sources["expanded"]
 windows=sources["windows"]
 legacy=sources["coloc"]
 replay=sources["replay"]
 verified_replay=sources.get("verified_replay")
 priority=sources["priority8"]
 lit=sources["literature"]
 cell=sources["cell_gene"]
 feature=sources["cell_feature"]
 pairs=sources["cell_pairs"]
 u=sources["universe"]
 positional=[r for r in expanded if r["gene_id_stable"]]
 legacy_only=[r for r in expanded if not r["gene_id_stable"]]
 if len(positional)!=2425 or len(legacy_only)!=1 or len(expanded)!=2426:
  raise ValueError("Expanded positional universe changed")
 if len(windows)!=2425 or len({exact_row_key(r) for r in positional})!=2425:
  raise ValueError("Duplicate/missing group and gene ID in positional source")
 if len({(r["group_id"],stable(r["gene_id"])) for r in windows})!=2425:
  raise ValueError("Duplicate GENCODE window mapping")
 win={(r["group_id"],stable(r["gene_id"])):r for r in windows}
 if set(win)!={exact_row_key(r) for r in positional}:
  raise ValueError("Expanded gene ID versus gene-window joins changed")
 if len({r["gene_id_stable"] for r in positional})!=2225 or u["unique_gene_ids"]!=2225:
  raise ValueError("Expected 2225 distinct positional Ensembl stable IDs")
 if len({r["group_id"] for r in positional})!=80 or u["regions"]!=80:
  raise ValueError("Region universe changed")

 # GWAS-locus QTL provenance includes phenotype and tissue; keep all records
 # but attach only to explicit EAS locus and matching stable ID.
 if len(legacy)!=646 or len(replay)!=646:
  raise ValueError("Expected 646 source ABF and replay entries")
 locus_by_id=defaultdict(list)
 qtl_by_key={}
 for r in legacy:
  key=r["locus"],r["dataset_key"],r["gene_base"]
  if key in qtl_by_key:raise ValueError("Duplicate archived QTL assay key")
  qtl_by_key[key]=r
  locus_by_id[(r["locus"],r["gene_base"])].append(r)
 replay_by_key={}
 for r in replay:
  key=r["locus"],r["dataset_key"],r["gene_base"]
  if key in replay_by_key:raise ValueError("Duplicate replay assay key")
  replay_by_key[key]=r
 if set(qtl_by_key)!=set(replay_by_key):
  raise ValueError("Archived QTL assay and replay sources disagree")
 verified_by_key={}
 if verified_replay is not None:
  if len(verified_replay)!=646:raise ValueError("Verified replays must have 646 records")
  for item in verified_replay:
   key=item["locus"],item["dataset_key"],item["gene_base"]
   if key in verified_by_key:raise ValueError("Duplicate new replay assay")
   if item["status"]=="PASS" and (fnum(item["delta_max"]) is None or fnum(item["delta_max"])>1e-8):
    raise ValueError("Invalid PASS posterior delta")
   verified_by_key[key]=item
  if set(verified_by_key)!=set(qtl_by_key):raise ValueError("New replay assay keys changed")

 # Historic gene symbols are aliases; convert all relevant priority/literature/
 # healthy-cell tags through the exact archived gene_ID (locus+dataset+symbol).
 def resolve_symbol_and_locus(symbol,locus):
  ids={r["gene_base"] for r in legacy if r["locus"]==locus and r["gene_symbol"]==symbol}
  return sorted(ids)
 priority_by_id=defaultdict(list)
 priority_unmapped=[]
 for r in priority:
  if r["status"]!="ORIGINAL_FIVE_HYPOTHESES_SNP_REPLAY_PASS_NO_CAUSAL_CLAIM":
   raise ValueError("Unexpected priority source QC state")
  ids={q["gene_base"] for q in legacy if q["locus"]==r["locus"] and
       q["dataset_key"]==r["dataset_key"] and q["gene_symbol"]==r["gene"]}
  if len(ids)!=1:
   priority_unmapped.append({"category":"PRIORITY8","label":r["gene"],
       "legacy_locus":r["locus"],"status":"AMBIGUOUS_GENE_ID","candidates":";".join(sorted(ids))})
  else:priority_by_id[(r["locus"],next(iter(ids)))].append(r)
 if len(priority)!=8:raise ValueError("Expected eight priority assays")
 literature_by_id=defaultdict(list)
 for r in lit:
  ids=resolve_symbol_and_locus(r["gene"],r["locus"])
  if len(ids)!=1:
   priority_unmapped.append({"category":"LITERATURE","label":r["gene"],"legacy_locus":r["locus"],
                              "status":"AMBIGUOUS_GENE_ID","candidates":";".join(ids)})
  else:
   literature_by_id[(r["locus"],ids[0])].append(r)
 aliases=defaultdict(set)
 for r in legacy:aliases[r["gene_base"]].add(r["gene_symbol"])
 # Human reference annotations can attach by an exact current stable-annotation
 # symbol or by unambiguous archived GTEx stable-ID mapping.
 symbols_by_id=defaultdict(set)
 for r in positional: symbols_by_id[r["gene_id_stable"]].add(r["gene_symbol"])
 candidate_ids=set(symbols_by_id)
 symbol_mapping={}
 for r in cell:
  symbol=r["gene"]
  exact={gid for gid,names in symbols_by_id.items() if symbol in names}
  legacy_alias={gid for gid in candidate_ids if symbol in aliases[gid]}
  union=exact|legacy_alias
  if len(union)==1:
   symbol_mapping[symbol]=(next(iter(union)),
                            "GENCODE_SYMBOL_MATCH" if exact else "ARCHIVED_GTEX_GENE_ID_ALIAS_MATCH")
  else:
   priority_unmapped.append({"category":"HUMAN_REFERENCE","label":symbol,"legacy_locus":"",
                             "status":"NO_UNIQUE_GENE_ID","candidates":";".join(sorted(union))})
 cell_by_id={}
 for r in cell:
  if r["gene"] in symbol_mapping:
   gid,method=symbol_mapping[r["gene"]]
   if gid in cell_by_id:raise ValueError("Ambiguous two healthy reference genes map to one stable ID")
   cell_by_id[gid]=(r,method)
 feature_by_symbol={r["gene"]:r for r in feature}
 paired_by_symbol=defaultdict(list)
 for r in pairs: paired_by_symbol[r["gene"]].append(r)

 contextual=[]
 for r in positional:
  k=exact_row_key(r)
  w=win[k]
  if stable(w["gene_id"])!=r["gene_id_stable"] or w["gene_symbol"]!=r["gene_symbol"]:
   raise ValueError("Window gene-symbol/id mismatch")
  loc=r["legacy_original_bbj_locus"]
  archived=(locus_by_id.get((loc,r["gene_id_stable"]),[])
       if r["ancestry"]=="EAS_AND_JAPANESE" and loc and
          r["legacy_coloc_status"]=="LOCUS_AND_GENE_ID_MATCHED" else [])
  statuses=Counter(replay_by_key[(q["locus"],q["dataset_key"],q["gene_base"])]["status"]
                   for q in archived)
  if archived and statuses.total()!=len(archived):raise ValueError("Replay sample count unexpected")
  recomputed=Counter(verified_by_key[(q["locus"],q["dataset_key"],q["gene_base"])]["status"]
                 for q in archived) if verified_by_key else Counter()
  if r["eqtl_status"]=="REASSESS_LEGACY" and not archived:
   raise ValueError(f"Legacy eQTL label lacks matching locus+Ensembl ID: {k}")
  max_h4=max((fnum(q["PP.H4"]) for q in archived if fnum(q["PP.H4"]) is not None),default=None)
  priority_rows=priority_by_id.get((loc,r["gene_id_stable"]),[]) if archived else []
  literature_rows=literature_by_id.get((loc,r["gene_id_stable"]),[]) if archived else []
  contextual.append({
   **r,
   "genomic_chr_GRCh37":w["chr"],
   "gene_start_GRCh37":w["gene_start"],"gene_end_GRCh37":w["gene_end"],
   "genomic_annotation_version":w["annotation_version"],
   "archive_QTL_scope":"EAS_LEGACY_EXACT_LOCUS_GENE_ID" if archived else "NOT_APPLICABLE_OR_NOT_TESTED",
   "archived_ABF_assays":len(archived),"archived_max_PP_H4":max_h4 if max_h4 is not None else "",
   "archived_replay_PASS":statuses.get("PASS",0),
   "archived_replay_MISSING_INPUT":statuses.get("MISSING_INPUT",0),
   "new_independent_SNP_ABF_replay_PASS":recomputed.get("PASS",0),
   "new_independent_SNP_ABF_replay_REVIEW":len(archived)-recomputed.get("PASS",0) if verified_by_key else "",
   "priority8_SNP_replay":len(priority_rows),
   "literature_assays_linked":len(literature_rows),
   "legacy_symbols":";".join(sorted(aliases.get(r["gene_id_stable"],set()))),
   "interpretation":"DESCRIPTIVE_EVIDENCE_NOT_CAUSAL_GENE",
  })

 per=defaultdict(list)
 for row in contextual:per[row["gene_id_stable"]].append(row)
 full=[]
 for gid,rs in sorted(per.items()):
  ancestries=Counter(x["ancestry"] for x in rs)
  primary=sorted(rs,key=lambda x: (x["ancestry"]!="EAS_AND_JAPANESE",x["gene_symbol"]))[0]
  symbols=sorted({x["gene_symbol"] for x in rs})
  grps=sorted({x["group_id"] for x in rs})
  coding=sorted({x["biotype"] for x in rs})
  pvalues=[fnum(x["region_best_p"]) for x in rs if fnum(x["region_best_p"]) is not None]
  distances=[int(x["gene_lead_distance_bp"]) for x in rs if x["gene_lead_distance_bp"]]
  qassays=sum(x["archived_ABF_assays"] for x in rs)
  qpass=sum(x["archived_replay_PASS"] for x in rs)
  qmissing=sum(x["archived_replay_MISSING_INPUT"] for x in rs)
  verified_pass=sum(x["new_independent_SNP_ABF_replay_PASS"] for x in rs)
  pri=sum(x["priority8_SNP_replay"] for x in rs)
  lcount=sum(x["literature_assays_linked"] for x in rs)
  h4s=[x["archived_max_PP_H4"] for x in rs if x["archived_max_PP_H4"]!=""]
  cellval,method=cell_by_id.get(gid,(None,""))
  if cellval:
   symbol=cellval["gene"]
   feat=feature_by_symbol.get(symbol,{}).get("feature_status","NOT_ASSESSED")
   healthy_status=feat
   top_cell=cellval["top_celltype"] if cellval["present"]=="TRUE" else ""
   donor_pairs=paired_by_symbol[symbol]
   donors_pass=sum(x["qc_status"]=="DESCRIPTIVE_PAIRED_QC_PASS" for x in donor_pairs)
  else:
   healthy_status="NOT_ASSESSED_IN_SELECTED_NINE_GENE_HUMAN_REFERENCE"
   top_cell="";donors_pass=0
  gws=sum(x["region_gws"]=="1" for x in rs)
  level=("SNP_ABF_REPLAY_CONTEXT_READY" if (verified_pass if verified_by_key else qpass)>0 else
         "LEGACY_ABF_INPUT_PENDING" if qassays>0 else
         "GWAS_REGION_POSITIONAL_FOLLOWUP" if gws>0 else
         "POSITIONAL_OR_SUGGESTIVE_FOLLOWUP")
  full.append({
   "gene_id_stable":gid,"canonical_GENCODE_v19_gene_symbol":primary["gene_symbol"],
   "all_positional_symbols":";".join(symbols),
   "archived_GTEx_gene_symbols":";".join(sorted(aliases.get(gid,set()))),
   "biotypes":";".join(coding),
   "n_positional_regions":len(rs),
   "region_ids":";".join(grps),
   "n_EAS_Japanese_regions":ancestries.get("EAS_AND_JAPANESE",0),
   "n_EUR_regions":ancestries.get("EUR",0),
   "n_GWS_regions":gws,
   "minimum_region_lead_p":min(pvalues) if pvalues else "",
   "minimum_gene_to_lead_distance_bp":min(distances) if distances else "",
   "archived_EAS_GTEx_ABF_assays":qassays,
   "archived_EAS_GTEx_max_PP_H4":max(h4s) if h4s else "",
   "archived_EAS_GTEx_SNP_replay_PASS":qpass,
   "archived_EAS_GTEx_SNP_replay_MISSING_INPUT":qmissing,
   "new_SNP_ABF_replay_PASS":verified_pass if verified_by_key else "",
   "SNP_replayed_priority8_evidence":pri,
   "linked_literature_gene_locus_records":lcount,
   "human_reference_feature_status":healthy_status,
   "human_reference_mapping_method":method,
   "human_reference_top_celltype":top_cell,
   "human_reference_donor_paired_QC_PASS":donors_pass,
   "sQTL_status":"NOT_ASSESSED",
   "pQTL_status":"NOT_ASSESSED",
   "direct_causal_gene_status":"NOT_ESTABLISHED",
   "readiness_for_exploratory_followup":level,
   "scientific_interpretation":"TRIAGE_ONLY_NOT_CAUSAL_PROOF",
  })
 if len(full)!=2225 or sum(int(x["n_positional_regions"]) for x in full)!=2425:
  raise ValueError("Final gene-level partition broke the 2225/2425 scope")
 if sum(int(x["archived_EAS_GTEx_ABF_assays"]) for x in full)>646:
  raise ValueError("Accidentally double counted historic QTL assays")

 # NEURL1 legacy anchor without stable ID exists in historical version: capture
 # exact GTEx stable gene_ID as separate alias linkage without expanding scope.
 anchor=[]
 for r in legacy_only:
  sy=r["gene_symbol"]
  matches=sorted({q["gene_base"] for q in legacy if q["gene_symbol"]==sy})
  status=("UNIQUE_LEGACY_GTEX_ID_NOT_POSITIONAL" if len(matches)==1 and
          matches[0] not in per else
          "ALIASED_TO_EXISTING_POSITIONAL_STABLE_ID" if len(matches)==1 else
          "LEGACY_SYMBOL_UNRESOLVED")
  anchor.append({"legacy_symbol":sy,"legacy_source":r["source"],
                 "legacy_candidate_stable_ids":";".join(matches),
                 "stable_id_already_in_2225":"YES" if any(g in per for g in matches) else "NO",
                 "status":status,"included_in_2225_gene_count":"NO_ADDITIONAL_GENE"})
  # Attaching NEURL1 symbol directly to NEURL stable ID is permitted only via GTEx ID.
  if sy=="NEURL1" and not (len(matches)==1 and
                            matches[0]=="ENSG00000107954" and matches[0] in per):
   raise ValueError("NEURL1 archived-GTEx-ID-to-GENCODE mapping changed")

 diagnostics={
  "status":"LOCALLY_RECONCILED_EXISTING_DATA_ONLY",
  "stable_gene_ids":len(full),"gene_region_pairs":len(contextual),
  "ancestry_pairs":dict(Counter(r["ancestry"] for r in contextual)),
  "priority8_archived_assays":len(priority),
  "priority8_attached_assays":sum(x["priority8_SNP_replay"] for x in contextual),
  "historic_ABF_records":len(legacy),
  "historic_replay_status":dict(Counter(r["status"] for r in replay)),
  "new_SNP_ABF_replay_status":dict(Counter(r["status"] for r in verified_replay)) if verified_by_key else "NOT_PROVIDED",
  "linked_historical_ABF_assays_to_EAS_positional_gene_regions":sum(x["archived_ABF_assays"] for x in contextual),
  "linked_replay_pass_to_EAS_positional_gene_regions":sum(x["archived_replay_PASS"] for x in contextual),
  "new_replay_pass_linked_to_EAS_genes":sum(x["new_independent_SNP_ABF_replay_PASS"] for x in contextual),
  "legacy_only_anchors":anchor,
  "unmapped_support_context":priority_unmapped,
  "gene_followup_tiers":dict(Counter(x["readiness_for_exploratory_followup"] for x in full)),
  "human_gene_id_mappings":{symbol:{"gene_id":gid,"method":how} for symbol,(gid,how) in symbol_mapping.items()},
  "human_reference_not_disease_comparison":True,
  "legacy_coloc_not_applied_to_EUR":True,
  "cohort_matched_LD_verified":False,
  "no_new_DB":True,
  "ALDH2_alcohol_finemap_status":"BLOCKED",
  "ADH1B_alcohol_finemap_status":"EXPLORATORY",
  "causal_gene_established":0,
  "no_hard_candidate_filtering":True,
 }
 return contextual,full,anchor,diagnostics

def main():
 p=argparse.ArgumentParser()
 p.add_argument("--results-root",type=Path,default=Path("/srv/is-analysis/results/is"))
 p.add_argument("--out-dir",required=True,type=Path)
 p.add_argument("--verified-replay",type=Path,help="Full source-file-based SNP-level coloc rerun TSV; overrides readiness only, never canonical sources")
 a=p.parse_args()
 paths={name:a.results_root/rel for name,rel in REL.items()}
 for k,v in paths.items():
  if not v.exists():raise FileNotFoundError(str(v))
 sources={k:json.loads(v.read_text()) if v.suffix==".json" else read(v) for k,v in paths.items()}
 if a.verified_replay:
  sources["verified_replay"]=read(a.verified_replay)
  paths["verified_replay"]=a.verified_replay
 contextual,genes,anchors,manifest=build(sources)
 if a.out_dir.exists() and any(a.out_dir.iterdir()):
  raise FileExistsError("Will not overwrite existing results")
 a.out_dir.mkdir(parents=True,exist_ok=True)
 tswrite(a.out_dir/"IS_2425_GENE_REGION_EVIDENCE.tsv",contextual,list(contextual[0]))
 tswrite(a.out_dir/"IS_2225_GENE_EVIDENCE_MATRIX.tsv",genes,list(genes[0]))
 tswrite(a.out_dir/"IS_LEGACY_ONLY_ANCHORS.tsv",anchors,list(anchors[0]))
 manifest["sources_sha256"]={k:digest(v) for k,v in paths.items()}
 manifest["sources"]={k:str(v) for k,v in paths.items()}
 (a.out_dir/"IS_GENE_EVIDENCE_MATRIX_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
 print("IS_EVIDENCE_MATRIX_COMPLETE","GENES",len(genes),"REGION_PAIRS",len(contextual))
 print("REPLAY",manifest["historic_replay_status"])
 print("NEW_REPLAY",manifest["new_SNP_ABF_replay_status"])
 print("LINKED_ABF",manifest["linked_historical_ABF_assays_to_EAS_positional_gene_regions"])
 print("PRIORITY8_LINKED",manifest["priority8_attached_assays"])
 print("HUMAN",manifest["human_gene_id_mappings"])
 print("TOP_CONTEXT")
 for g in genes:
  if g["SNP_replayed_priority8_evidence"]:
   print(g["canonical_GENCODE_v19_gene_symbol"],g["gene_id_stable"],
         "GTEx_ALIAS",g["archived_GTEx_gene_symbols"],
         "ASSAYS",g["archived_EAS_GTEx_ABF_assays"],
         "REPLAY",g["archived_EAS_GTEx_SNP_replay_PASS"])

if __name__=="__main__":main()
