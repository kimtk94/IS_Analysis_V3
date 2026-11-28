#!/usr/bin/env python3
"""Fixed PNAS 67-target audit of original selection estimates; no panel refit."""
import argparse,collections,hashlib,json,math,zipfile
from pathlib import Path
from publication_sprint import POOL_SHA,read_gene_table,same_sign,sha
from ocular_frozen_241_external import read_panel,PIN
from ocular_p5j_official_TM_overlap import parse_official,TM_ALL
def req(ok,msg):
    if not ok:raise ValueError(msg)
def run(a):
    x=json.loads(a.evidence.read_text())
    d=x["result"]
    req(x["provenance"]["run_id"]==45578116782 and x["provenance"]["artifact_id"]==15004188841 and d["status"]=="PUBLISHED_COUNTS_RECONCILED_EXACT_ID_OVERLAP","PNAS provenance")
    ids=d["source_gene_ids"]
    req(len(ids)==len(set(ids))==67 and d["groups"]["all_all"]["panel_overlap_count"]==0,"Frozen PNAS targets")
    req(sha(a.pool)==POOL_SHA and sha(a.panel)==PIN,"Original source archive hash")
    panel=read_panel(a.panel)
    with zipfile.ZipFile(a.pool) as z:
        donor=read_gene_table(z,"p5_all_gene_results.tsv.gz","logFC","FDR",16050)
        pool=read_gene_table(z,"p5d_poolunit_all_gene_results.tsv.gz","pool_logFC","pool_FDR",21081)
    out={"status":"FIXED_67_SELECTION_LINEAGE_AUDIT","target_count":67,"pool_archive_sha256":POOL_SHA,"panel_sha256":PIN,"source_67_run":45578116782,"original_models":{},"selection_48h":{"status":"HOLD_SOURCE_UNVERIFIED"},"limits":["All 67 exact stable Ensembl IDs, not P-selected. Original donor and pool are the same GSE313736 cohort, not independent studies. Missing gene not zero.","Original selection-stage 48h source included only if exact prior byte hash matches. Different tested universes and models; no enrichment P or confirmation of glaucoma mechanism. No changes to fixed 241/16 membership.","Descriptive original selection lineage; no external validation or causal inference."]}
    shared=set(donor)&set(pool);both={g for g in shared if donor[g][1]<.05 and pool[g][1]<.05 and same_sign(donor[g][0],pool[g][0])};strong={g for g in both if min(abs(donor[g][0]),abs(pool[g][0]))>=1}
    for name,group in (("donor_tested",set(donor)),("pool_tested",set(pool)),("both_tested",shared),("donor_fdr05",{g for g,v in donor.items() if v[1]<.05}),("pool_fdr05",{g for g,v in pool.items() if v[1]<.05}),("both_fdr05_same_sign",both),("both_fdr05_same_sign_abs_ge1",strong),("frozen241",set(panel))):
        out["original_models"][name]={"universe_count":len(group),"pnas_targets":len(set(ids)&group),"target_ids":sorted(set(ids)&group)}
    records=[]
    for g in ids:
        r={"gene_id":g,"in_frozen241":g in panel}
        for name,src in (("donor",donor),("pool",pool)):
            r[name]={"symbol":src[g][2],"log2fc":src[g][0],"fdr":src[g][1]} if g in src else None
        r["original_internal_status"]="NOT_IN_BOTH_MODELS" if g not in shared else "NOT_BOTH_SIGNIFICANT_AND_CONCORDANT" if g not in both else "CONCORDANT_BELOW_STRONG_THRESHOLD" if g not in strong else "INTERNAL_STRONG_ELIGIBLE"
        records.append(r)
    out["internal_status_counts"]=dict(collections.Counter(r["original_internal_status"] for r in records))
    # Optional same-byte selection-stage table only. Do not replace with a new version.
    if a.external and a.external.exists() and a.external_hash:
        observed=sha(a.external)
        if observed==a.external_hash:
            try:
                external=parse_official(a.external)[TM_ALL];three=set(external)&shared
                req(len(three)==12251,"THREE_MODEL_GENE_COUNT")
                eligible={g for g in three&both if external[g]["fdr"] is not None and external[g]["fdr"]<.05 and same_sign(donor[g][0],external[g]["log2fc"])}
                req(eligible&strong==set(panel),"FROZEN_PANEL_PARITY")
                out["selection_48h"]={"status":"PINNED_SELECTION_SOURCE_RECONCILED","sha256":observed,"three_model_tested":len(three),"three_model_concordant_fdr05":len(eligible),"targets_in_three_model":len(set(ids)&three),"targets_concordant_fdr05":len(set(ids)&eligible),"targets_strong":len(set(ids)&eligible&strong)}
                for r in records:r["selection_48h"]=external.get(r["gene_id"])
            except (ValueError,KeyError,OSError,zipfile.BadZipFile) as ex:out["selection_48h"]={"status":"HOLD_SCHEMA_OR_PANEL_PARITY","sha256":observed,"reason":str(ex)[:150]}
        else:out["selection_48h"]={"status":"HOLD_SOURCE_HASH_CHANGED","observed_sha256":observed,"expected_sha256":a.external_hash}
    out["records"]=records
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--evidence",type=Path,required=True);p.add_argument("--pool",type=Path,required=True);p.add_argument("--panel",type=Path,required=True);p.add_argument("--external",type=Path);p.add_argument("--external-hash");p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    try:x=run(a)
    except (OSError,ValueError,KeyError,TypeError,zipfile.BadZipFile) as ex:x={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(ex)[:250]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps({k:v for k,v in x.items() if k!="records"},indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
