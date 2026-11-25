#!/usr/bin/env python3
"""Fixed publisher 67 targets in same-study GSE301517 CPM. Descriptive only."""
import argparse,collections,csv,gzip,json,math,statistics
from pathlib import Path
from ocular_frozen_241_external import CPM_PIN,STRAINS,get_source,sha
def req(ok,msg):
    if not ok:raise ValueError(msg)
def read_source(path,ids):
    req(sha(path)==CPM_PIN,"Source hash")
    found={};n=0
    with gzip.open(path,"rt",encoding="utf8",newline="") as f:
        t=csv.DictReader(f,delimiter="\t")
        required={""}|{f"{s}_{arm}_{i}" for s,k in STRAINS.items() for arm in ("Control","Dex") for i in range(1,k+1)}
        req(required<=set(t.fieldnames or []) and len(t.fieldnames)==len(set(t.fieldnames)),"CPM header")
        for row in t:
            n+=1;req(n<=85298,"Source row count cap")
            g=row[""].split(".")[0]
            if g not in ids:continue
            req(g not in found,"Duplicate target ID")
            means={}
            for s,k in STRAINS.items():
                for arm in ("Control","Dex"):
                    vals=[float(row[f"{s}_{arm}_{i}"]) for i in range(1,k+1)]
                    req(all(math.isfinite(v) and v>=0 for v in vals),"Nonfinite or negative CPM")
                    means[s,arm]=statistics.mean(vals)
            eff={s:math.log2((means[s,"Dex"]+1)/(means[s,"Control"]+1)) for s in STRAINS}
            found[g]={"strain_log2_ratios":eff,"mean_log2_ratio":statistics.mean(eff.values()),"min_arm_mean_cpm":min(statistics.mean(means[s,arm] for s in STRAINS) for arm in ("Control","Dex"))}
    req(n==85298,"Unexpected row count")
    return found,n
def run(a):
    raw=json.loads(a.lineage.read_text());e=json.loads(a.evidence.read_text())
    req(e["provenance"]["run_id"]==45578473621 and e["provenance"]["artifact_id"]==15004819277,"Lineage provenance")
    req({k:v for k,v in raw.items() if k!="records"}==e["result"] and raw["status"]=="FIXED_67_SELECTION_LINEAGE_AUDIT","Aggregate parity")
    d={v["gene_id"]:v for v in raw["records"]}
    req(len(d)==len(raw["records"])==67 and not any(v["in_frozen241"] for v in d.values()),"Fixed distinct targets")
    ref={k:v["donor"]["log2fc"] for k,v in d.items() if v["donor"] is not None and v["donor"]["log2fc"]!=0}
    both=set(raw["original_models"]["both_fdr05_same_sign"]["target_ids"]);strong=set(raw["original_models"]["both_fdr05_same_sign_abs_ge1"]["target_ids"])
    req(len(both)==37 and len(strong)==4 and strong<=both<=set(ref),"Fixed subgroups")
    url,attempts=get_source(a.cpm)
    if not url:return {"status":"HOLD_7DAY_SOURCE","attempts":attempts}
    obs,n=read_source(a.cpm,set(d))
    rows=[]
    for g,v in sorted(d.items()):
        r={"gene_id":g,"donor":v["donor"],"pool":v["pool"],"internal_status":v["original_internal_status"],"observed_7d":obs.get(g)}
        if g in ref and g in obs:
            sign=1 if ref[g]>0 else -1
            x=obs[g];signed=[sign*z for z in x["strain_log2_ratios"].values()]
            r["directional"]={"signed_mean":sign*x["mean_log2_ratio"],"same_mean":sign*x["mean_log2_ratio"]>0,"all_three_same":all(z>0 for z in signed),"all_three_abs_ge025":all(z>=.25 for z in signed),"leave_one_strain_out_all_same":all(statistics.mean(signed[:i]+signed[i+1:])>0 for i in range(3))}
        rows.append(r)
    out={"status":"FIXED_PNAS67_7DAY_DESCRIPTIVE","source":url,"source_sha256":CPM_PIN,"source_rows":n,"all_targets":67,"mapped":len(obs),"missing_ids":sorted(set(d)-set(obs)),"direction_reference_available":len(ref),"groups":{},"internal_strong_4":[r for r in rows if r["gene_id"] in strong],"limits":["GSE301517 belongs to the same PNAS study GSE301525 that defined the 67 promoter-contact targets; not independent validation of the PNAS study. Donor/pool direction uses the same GSE313736 discovery cohort. Subgroups fixed from already observed internal model results, not prospective confirmation.","Report all 67; no imputed direction for untested genes, no missing=zero. CPM normalized values; strain technical repeats averaged, equal strain weights; log2((DEX mean+1)/(control mean+1)). No gene-level inferential P, pooled study test, causal or clinical prediction claim.","Separate exploratory target audit, not a change to the frozen241/ECM16 panel. No claim of expression specificity, effect equality or time trend."]}
    for name,keys in (("all67",set(d)),("donor_direction",set(ref)),("both_significant_concordant37",both),("internal_strong4",strong)):
        group=[v for v in rows if v["gene_id"] in keys];strata=[]
        for cutoff in (0,1,5):
            m=[v for v in group if v["observed_7d"] is not None and v["observed_7d"]["min_arm_mean_cpm"]>=cutoff]
            z=[v["directional"] for v in m if "directional" in v]
            strata.append({"min_arm_mean_cpm":cutoff,"mapped_expression_eligible":len(m),"direction_evaluable":len(z),"mean_direction_concordant":sum(v["same_mean"] for v in z),"all_three_direction_concordant":sum(v["all_three_same"] for v in z),"all_three_abs_ge025":sum(v["all_three_abs_ge025"] for v in z),"all_leave_one_strain_out_same":sum(v["leave_one_strain_out_all_same"] for v in z),"median_signed_mean_log2_ratio":statistics.median(v["signed_mean"] for v in z) if z else None})
        out["groups"][name]={"total":len(group),"strata":strata}
    out["records"]=rows
    return out
def main():
    p=argparse.ArgumentParser()
    for k in ("lineage","evidence","cpm","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    try:x=run(a)
    except (OSError,ValueError,KeyError,TypeError) as ex:x={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(ex)[:250]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps({k:v for k,v in x.items() if k!="records"},indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
