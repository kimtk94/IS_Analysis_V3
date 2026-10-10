#!/usr/bin/env python3
"""Fixed-seed ancestry panel size control: 1KG EAS104 subsampling.

Use identical source SNPs/dosage as JPT104. Sample 104 subjects without
replacement from 504 EAS individuals for independent repeat controls.
This isolates reference-N uncertainty from selected JPT ancestry, but does
not recreate GWAS-study LD nor provide causal inference. Repeats are
correlated and descriptive, not an inferential permutation p test.
"""
import argparse,csv,json
from pathlib import Path
import numpy as np
from prepare_alcohol_jpt104_ld_sensitivity import BASE,SAMPLES,LOCI,read_rows
def run(base=BASE,sample_file=SAMPLES,reps=12,seed=20261010):
    if reps<2 or reps>30:raise ValueError("Replicates restricted to 2-30")
    sample_set=set(sample_file.read_text().splitlines())
    if len(sample_set)!=504:raise ValueError("EAS source individuals n504 unverified")
    gen=np.random.default_rng(seed)
    index_sets=None;result={}
    for locus in LOCI:
        d=base/locus
        v=read_rows(d/"variants.tsv")
        names=[x["ID"] for x in v]
        arr={}
        with (d/"reference_eas504.traw").open() as f:
            rr=csv.reader(f,delimiter="\t")
            header=next(rr)
            if header[:6]!=["CHR","SNP","(C)M","POS","COUNTED","ALT"]:
                raise ValueError("PLINK traw wrong schema")
            ids=[x.split("_")[0] for x in header[6:]]
            if len(ids)!=504 or len(set(ids))!=504 or set(ids)!=sample_set:
                raise ValueError("Wrong 504 source samples")
            if index_sets is None:
                # Fixed shared draws across all loci. No post-hoc selection on GWAS s.
                index_sets=[np.sort(gen.choice(504,104,replace=False))
                            for _ in range(reps)]
                provenance_ids=ids
            elif ids!=provenance_ids:
                raise ValueError("Cross-chrom sample order changed")
            for row in rr:
                snp=row[1];counted,other=row[4:6]
                if snp in arr:raise ValueError("Duplicate source genotyped SNP")
                ref,alt=snp.split(":")[2:]
                if {counted,other}!={ref,alt}:raise ValueError("Wrong reference allele set")
                doses=np.array(row[6:],dtype=float)
                if len(doses)!=504 or not np.isfinite(doses).all():
                    raise ValueError("Missing genotype in reference sample")
                if counted!=alt:doses=2-doses
                arr[snp]=doses
        if set(arr)!=set(names):raise ValueError("Source SNP counts not same as original matrix")
        x=np.vstack([arr[vid] for vid in names])
        if np.std(x,axis=1,ddof=1).min()<=0:raise ValueError("Full EAS monomorphic")
        n=len(names);entries=[]
        out=d/"EAS104_random_reference_controls"
        out.mkdir(exist_ok=True)
        for it,subset in enumerate(index_sets,1):
            g=x[:,subset]
            stdev=g.std(axis=1,ddof=1)
            bad=np.flatnonzero(stdev<1e-10)
            key=f"reference_draw_{it:02d}"
            if len(bad):
                # record rather than impute/pretend r2=0 for missing-LD
                entries.append({"id":key,"sample_count":104,
                   "n_SNPs_requested":n,"n_SNPs_monomorphic":len(bad),
                   "status":"MONOMORPHIC_ABORT","file":None})
                continue
            g=(g-g.mean(axis=1,keepdims=True))/stdev[:,None]
            R=(g@g.T)/103
            if np.max(abs(np.diag(R)-1))>1e-6:raise ValueError("Broken draw correlation")
            file=out/f"{key}.ld.f64.rowmajor"
            file.write_bytes(R.astype("<f8").tobytes())
            entries.append({"id":key,"sample_count":104,
                   "n_SNPs_requested":n,"n_SNPs_monomorphic":0,
                   "status":"COMPLETE","file":file.name})
        plan={"locus":locus,"source_1000G_EAS_n":504,"draw_n":104,
          "draws_requested":reps,"reproducibility_seed":seed,
          "same_sample_draw_indices_across_ADH1B_ALDH2":True,
          "draw_sampling_without_replacement":True,
          "random_draws_are_not_independent_GWAS_datasets":True,
          "no_causal_claim":True,"source_SNPs":n,"draws":entries}
        (out/"random_reference_manifest.json").write_text(json.dumps(plan,indent=2))
        print("EAS104_DOWNSAMPLE",locus,"draws_complete",sum(a["status"]=="COMPLETE" for a in entries),
          "monomorphic_abort",sum(a["status"]!="COMPLETE" for a in entries),flush=True)
        result[locus]=plan
    return result
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=BASE)
    p.add_argument("--replicates",type=int,default=12)
    p.add_argument("--seed",type=int,default=20261010)
    a=p.parse_args();run(a.root,reps=a.replicates,seed=a.seed)
