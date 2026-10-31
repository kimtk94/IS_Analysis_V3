#!/usr/bin/env python3
"""Gate B alternative cell labels: donor-paired feasibility without individual data.

Only aggregate public metadata; no raw read download and no genotype/phenotype.
Cell label is post-exposure and cannot remove causal treatment/library confounding.
"""
import argparse,csv,gzip,json,math,tempfile
from collections import defaultdict,Counter
from pathlib import Path
from ocular_p1_integrity import BASE,download
from ocular_p4_fullcohort import pool_name

TYPES=("BeamA_any","BeamA_050","BeamA_075","scpred_TM_Fibro","BeamA_050_AND_TM_Fibro",
       "BeamA_075_AND_TM_Fibro","JCT_any","JCT_075")
THRESHOLDS=(10,30,50)
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    donors=defaultdict(Counter);totals=defaultdict(Counter);donor_pool=defaultdict(set)
    observed=0
    with tempfile.TemporaryDirectory(prefix="ocular-b-alternates-") as td:
        filename=Path(td)/"metadata.gz"
        download(BASE+"series/GSE313nnn/GSE313736/suppl/GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",filename,12*1024*1024)
        with gzip.open(filename,"rt",encoding="utf8") as h:
            reader=csv.DictReader(h,delimiter="\t")
            required={"donor_id_demuxlet","droplet_type","library_id","DEX Treatment",
                      "celltypist_chen2025_fine","celltypist_chen2025_fine_score",
                      "scpred_vanzyl2022"}
            if not required.issubset(reader.fieldnames or []):raise ValueError("CELL_LABEL_SOURCE_SCHEMA_MISSING")
            for row in reader:
                observed+=1
                if row["droplet_type"]!="SNG":continue
                arm=row["DEX Treatment"].strip().lower()
                donor=row["donor_id_demuxlet"].strip()
                if arm not in ("dex","no-dex") or not donor or donor.lower() in ("na","doublet","unassigned"):
                    continue
                pool=pool_name(row["library_id"])
                if not pool:continue
                donor_pool[donor].add(pool)
                lab=row["celltypist_chen2025_fine"].strip()
                scpred=row["scpred_vanzyl2022"].strip()=="TM_Fibro"
                try:score=float(row["celltypist_chen2025_fine_score"])
                except (TypeError,ValueError):score=float("nan")
                isbeam=lab=="BeamA";isjct=lab=="JCT"
                accepted={"BeamA_any":isbeam,"BeamA_050":isbeam and score>=.5,
                  "BeamA_075":isbeam and score>=.75,
                  "scpred_TM_Fibro":scpred,
                  "BeamA_050_AND_TM_Fibro":isbeam and score>=.5 and scpred,
                  "BeamA_075_AND_TM_Fibro":isbeam and score>=.75 and scpred,
                  "JCT_any":isjct,"JCT_075":isjct and score>=.75}
                for name,val in accepted.items():
                    if val:
                        donors[(donor,arm)][name]+=1
                        totals[arm][name]+=1
        if observed!=169398:raise ValueError("FINAL_METADATA_ROWS_CHANGED")
        if any(len(sets)!=1 for sets in donor_pool.values()):
            raise ValueError("DONOR_ALLOCATED_ACROSS_MULTIPLE_POOLS")
        donor_ids={d for d,arm in donors}
        result=[]
        pools={}
        by30={}
        for subset in TYPES:
            stats={t:set(d for d in donor_ids
                  if donors[d,"dex"][subset]>=t and donors[d,"no-dex"][subset]>=t)
                  for t in THRESHOLDS}
            by30[subset]=stats[30]
            entry={"selection":subset,
              "all_sng_metadata_cells_vehicle":totals["no-dex"][subset],
              "all_sng_metadata_cells_dex":totals["dex"][subset],
              "paired_donors_min10_per_arm":len(stats[10]),
              "paired_donors_min30_per_arm":len(stats[30]),
              "paired_donors_min50_per_arm":len(stats[50]),
              "selected_metadata_cells_30_donors_vehicle":sum(donors[d,"no-dex"][subset] for d in stats[30]),
              "selected_metadata_cells_30_donors_dex":sum(donors[d,"dex"][subset] for d in stats[30]),
              "pool_count_min30":len({next(iter(donor_pool[d])) for d in stats[30]}),
              "excluded_due_to_post_exposure_annotation":False,
              "alternative_predefined_label_not_independent_replication":True}
            result.append(entry)
            pools[subset]=dict(Counter(next(iter(donor_pool[d])) for d in stats[30]))
        if len(by30["BeamA_075"])!=57 or len(by30["BeamA_050"])!=57:
            raise ValueError("BASELINE_BEAMA_57_CHANGED")
        baseline=by30["BeamA_075"]
        for item in result:
            item["selected_donors_overlap_primary_57_n"]=len(by30[item["selection"]]&baseline)
            item["new_or_different_donors_n"]=len(by30[item["selection"]]-baseline)
        with (a.out/"GATE_B_alternative_cell_annotation_donor_eligibility.tsv").open("w",encoding="utf8",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(result[0]),delimiter="\t")
            writer.writeheader();writer.writerows(result)
        summary={"status":"PASS_METADATA_ALTERNATE_ANNOTATION_ELIGIBILITY_AUDIT",
          "original_GEO_cell_rows":observed,"selection_classes":list(TYPES),
          "paired_independent_donors_min30_by_selection":{r["selection"]:r["paired_donors_min30_per_arm"] for r in result},
          "primary_BeamA_075_donor_count":57,
          "per_source_pool_eligible_donors_by_selection":pools,
          "response_labels_post_treatment":True,
          "no_donor_ids_or_per_donor_counts_uploaded":True,
          "gene_expression_significance_computed_here":False,
          "interpretation":"Only metadata feasibility. Conditions processed in separate 10x libraries; different annotation labels are not randomized biological strata. JCT low donor N limits inference."}
        (a.out/"GATE_B_alternative_cell_annotation_readiness.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
        print("ALTERNATIVE_ANNOTATION",json.dumps(summary["paired_independent_donors_min30_by_selection"],ensure_ascii=False),flush=True)
if __name__=="__main__":main()
