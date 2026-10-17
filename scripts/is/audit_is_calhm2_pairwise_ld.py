#!/usr/bin/env python3
"""Audited pairwise LD for CALHM2 QTL-tag SNPs in public 1000G EAS/JPT.

Not fine-mapping; no inference about causal number of signals.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import pandas as pd

LEADS=("10:105459834:T:C","10:105599770:T:G")
POPS=("EAS504","JPT104","otherEAS400")

def load_qtl(path, cutoff=1e-4):
    with Path(path).open(newline="") as f:
        records=list(csv.DictReader(f,delimiter="\t"))
    seen={}
    for r in records:
        if r["locus"]!="BBJ_IS_L002" or r["symbol"]!="CALHM2":
            continue
        if r["BBJ_exact_Fasta37"]!="True" or r["GIGA_EAS_AIS_exact_Fasta37"]!="True" or r["palindromic_review"]!="False":
            continue
        if float(r["QTD000161_nominal_p"])>=cutoff:
            continue
        k=r["SNP_GRCh37"]
        if k in seen and seen[k]!=r:
            raise ValueError("Duplicate QTL SNP with conflicting row")
        seen[k]=r
    if len(seen)!=3:
        raise ValueError(f"Expected exactly 3 CALHM2 SNPs; observed {len(seen)}")
    return seen

def oriented_dosage(raw, snps):
    df=pd.read_csv(raw,sep="\t")
    columns={col.rsplit("_",1)[0]:col for col in df.columns[6:]}
    if any(k not in columns for k in snps):
        raise ValueError("Missing target SNP genotype; no silent imputation")
    samples=df["IID"].astype(str)
    if samples.duplicated().any():
        raise ValueError("Duplicate individual IDs")
    data={}
    for k in snps:
        col=columns[k]
        counted=col.rsplit("_",1)[1]
        chrom,pos,ref,alt=k.split(":")
        values=df[col].to_numpy(float)
        if counted==ref:
            values=2-values
        elif counted!=alt:
            raise ValueError(f"Counted allele not REF or ALT: {k}")
        if not np.isfinite(values).all():
            raise ValueError(f"Missing genotypes: {k}")
        data[k]=values
    return samples.to_numpy(),data

def build(evidence,genotypes,panel,out):
    out=Path(out)
    if out.exists():
        raise FileExistsError("Output path exists; use a new immutable version")
    qtl=load_qtl(evidence)
    keys=sorted(qtl,key=lambda k:int(k.split(":")[1]))
    all_keys=list(LEADS)+keys
    if len(set(all_keys))!=len(all_keys):
        raise ValueError("Lead overlaps target key; review before running")
    samples,dosages=oriented_dosage(genotypes,all_keys)
    pop=pd.read_csv(panel,sep="\t")
    if pop["sample"].duplicated().any():
        raise ValueError("Population metadata duplicates")
    jpt=set(pop.loc[pop["pop"]=="JPT","sample"].astype(str))
    if len(samples)!=504 or len(jpt)!=104 or len(set(samples)&jpt)!=104:
        raise ValueError("Expected 504 EAS with 104 JPT individuals")
    select={"EAS504":np.ones(len(samples),dtype=bool),
            "JPT104":np.isin(samples,list(jpt))}
    select["otherEAS400"]=~select["JPT104"]
    results=[]
    for population in POPS:
        mask=select[population]
        if sum(mask)!={"EAS504":504,"JPT104":104,"otherEAS400":400}[population]:
            raise ValueError("Population N drift")
        for i,k1 in enumerate(all_keys):
            for k2 in all_keys[i+1:]:
                a=dosages[k1][mask];b=dosages[k2][mask]
                if np.std(a)==0 or np.std(b)==0:
                    raise ValueError("Monomorphic requested variant")
                r=float(np.corrcoef(a,b)[0,1])
                results.append(dict(population=population,n=int(sum(mask)),
                    variant_1=k1,variant_2=k2,signed_r=r,r2=r*r,
                    high_LD_r2_ge_0p8=bool(r*r>=.8),
                    SNP_pair_not_independent_causal_evidence=True))
    out.mkdir(parents=True)
    with (out/"CALHM2_THREE_SNP_PAIRWISE_LD.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter="\t")
        w.writeheader();w.writerows(results)
    manifest={"state":"REFERENCE_PAIRWISE_LD_ONLY",
      "source_evidence":str(evidence),"genotype_source":str(genotypes),
      "population_panel":str(panel),"qtl_snp_keys":keys,
      "populations":{p:int(sum(m)) for p,m in select.items()},
      "source_genotypes_unmodified":True,"not_credible_set":True,
      "not_causal_colocalization":True}
    (out/"CALHM2_PAIRWISE_LD_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
    return results

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--evidence",type=Path,required=True)
    p.add_argument("--genotypes",type=Path,required=True)
    p.add_argument("--panel",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(build(a.evidence,a.genotypes,a.panel,a.out),indent=2))
