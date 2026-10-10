#!/usr/bin/env python3
"""Reproduce regional Koyanagi alcohol GWAS / EAS LD PLINK2 sensitivity grid.

Default: plan-only. --execute writes results only under noncanonical research
results and LD reference directory, never into production pipelines.
"""
import argparse,json,subprocess
from pathlib import Path
from extract_alcohol_eas_region_reference_for_clump import SENTINELS
from audit_alcohol_eas_regional_clump_sensitivity import SRC,REF
def plan(root,ref):
    tasks=[]
    for maf in (None,.01,.05):
        for chrom in SENTINELS:
            name=f"chr{chrom}.EAS504.alcohol_r2_0p1" if maf is None else f"chr{chrom}.EAS504.alcohol_maf_{maf:.2f}_r2_0p1"
            args=["plink2","--pfile",str(ref/f"chr{chrom}.EAS504.regions")]
            if maf is not None:args+=["--maf",str(maf)]
            args+=["--clump",str(root/f"chr{chrom}.alcohol_regional.clump_input.tsv"),
                "--clump-p1","5e-8","--clump-p2","1e-4",
                "--clump-r2","0.1","--clump-kb","1000",
                "--out",str(root/name)]
            tasks.append({"chromosome":chrom,"maf":maf,"r2":.1,
                          "task":name,"command":args})
    for chrom in ("4","12"):
        for r2 in (.2,.5):
            name=f"chr{chrom}.EAS504.alcohol_maf_0.01_r2_{r2}"
            args=["plink2","--pfile",str(ref/f"chr{chrom}.EAS504.regions"),
                "--maf","0.01","--clump",str(root/f"chr{chrom}.alcohol_regional.clump_input.tsv"),
                "--clump-p1","5e-8","--clump-p2","1e-4",
                "--clump-r2",str(r2),"--clump-kb","1000",
                "--out",str(root/name)]
            tasks.append({"chromosome":chrom,"maf":.01,"r2":r2,
                          "task":name,"command":args})
    for chrom in SENTINELS:
        prefix=str(ref/f"chr{chrom}.EAS504.regions")
        tasks.append({"chromosome":chrom,"task":f"chr{chrom}.EAS504.freq",
                     "command":["plink2","--pfile",prefix,
                         "--freq","--out",prefix]})
    return tasks
def run(root,ref,execute):
    tasks=plan(root,ref)
    if not execute:
        print(json.dumps({"mode":"PLAN_ONLY","task_count":len(tasks),
          "chromosomes":list(SENTINELS),"sources_EAS_504":True,
          "MR_NOT_RUN":True,"commands":tasks},indent=2))
        return {"task_count":len(tasks)}
    if not(root/"IS_ALCOHOL_REGIONAL_GWAS_REFERENCE_COVERAGE.json").is_file():
        raise FileNotFoundError("Source GWAS allele-match audit absent")
    for chrom in SENTINELS:
        reference=json.loads((ref/f"chr{chrom}.region_source_audit.json").read_text())
        if not reference.get("validated") or reference.get("plink_sample_count")!=504:
            raise ValueError("Invalid regional 1KG EAS source")
        if not(root/f"chr{chrom}.alcohol_regional.clump_input.tsv").exists():
            raise FileNotFoundError("Source allele-harmonized clump associations missing")
    for task in tasks:
        print("RUN",task["task"],flush=True)
        p=subprocess.run(task["command"],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                         text=True,timeout=150)
        if p.returncode:
            raise RuntimeError(f"PLINK2 {task['task']} failed: {p.stdout[-900:]}")
        if "clump" in task["task"] and not Path(task["command"][-1]+".clumps").is_file():
            raise RuntimeError("Expected clump summary missing")
    from audit_alcohol_eas_regional_clump_sensitivity import run as audit
    result=audit(root,ref)
    if result["causal_eligible_MR_IVs"]!=0:
        raise ValueError("Invalid promotion of mere clumps to causal IVs")
    print("REPRODUCIBLE_REGIONAL_EAS_CLUMP_GRID_PASS",len(tasks))
    return result
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=SRC)
    p.add_argument("--reference",type=Path,default=REF)
    p.add_argument("--execute",action="store_true")
    a=p.parse_args();run(a.out,a.reference,a.execute)
