#!/usr/bin/env python3
"""Reference-only EUR503 LD blocks, rank and shrinkage diagnostics.

Blocks are contiguous genomic tiles, NOT evidence of independent causal
signals. All sample correlations rely on 1000G EUR, not GTEx in-cohort LD.
"""
import argparse,csv,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd

EXPECTED={"BBJ_IS_L001":2834,"BBJ_IS_L002":2828}
BLOCK_BP=250000
SHRINKAGE=(0.0,0.01,0.05,0.10,0.20)
def analyze(source,out):
    out.mkdir(parents=True,exist_ok=True)
    rows=[];global_summary=[]
    for loc,expected in EXPECTED.items():
        fp=source/f"{loc}.EUR503.raw"
        df=pd.read_csv(fp,sep="\t")
        if len(df)!=503:raise ValueError("EUR sample count drift")
        ids=[];dosage=df.iloc[:,6:].to_numpy(dtype=float,copy=True)
        for j,col in enumerate(df.columns[6:]):
            key,observed=col.rsplit("_",1)
            parts=key.split(":")
            if len(parts)!=4:raise ValueError("SNP key invalid")
            ref,alt=parts[-2:]
            if observed==ref:dosage[:,j]=2-dosage[:,j]
            elif observed!=alt:raise ValueError("Allele orientation not verified")
            ids.append((key,int(parts[1])))
        if len(ids)!=expected or len(set(x[0] for x in ids))!=expected:raise ValueError("Input SNP set drift")
        if np.any(~np.isfinite(dosage)):raise ValueError("Unexpected genotype missingness")
        if np.any((dosage<0)|(dosage>2)):raise ValueError("Invalid genotype")
        order=np.argsort([x[1] for x in ids])
        ids=[ids[k] for k in order];dosage=dosage[:,order]
        variant=np.var(dosage,axis=0)>0
        tiles={}
        for j,(_,pos) in enumerate(ids):
            tiles.setdefault((pos//BLOCK_BP),[]).append(j)
        good=0;largest=0
        for tile,indices in sorted(tiles.items()):
            inds=np.array([i for i in indices if variant[i]])
            if len(inds)<2:continue
            good+=1;largest=max(largest,len(inds))
            x=dosage[:,inds]
            x=(x-x.mean(axis=0))/x.std(axis=0,ddof=1)
            # Rank via sample dimension: number of SNPs can be > n_samples
            gram=(x@x.T)/(len(x)-1)
            eig=np.linalg.eigvalsh(gram)
            rank=int(np.sum(eig>max(eig.max(),1.0)*1e-9))
            # Gram positive eigenvalues equal nonzero SNP-LD eigenvalues.
            pos_eig=eig[eig>max(eig.max(),1.0)*1e-9]
            top=float(pos_eig.max());low=float(pos_eig.min())
            for alpha in SHRINKAGE:
                # (1-alpha) * R + alpha * I, so for nvars > rank,
                # lambda_min=alpha; otherwise use sample eigen minimum.
                nvars=len(inds)
                lower=(0.0 if nvars>rank else low)
                min_eig=alpha+(1-alpha)*lower
                max_eig=alpha+(1-alpha)*top
                cond=max_eig/min_eig if min_eig>0 else None
                rows.append(dict(locus=loc,block_index=tile,
                 chromosome=ids[inds[0]][0].split(":")[0],
                 start_b37=ids[inds[0]][1],end_b37=ids[inds[-1]][1],
                 n_reference_samples=503,n_polymorphic_snps=nvars,
                 effective_rank=rank,rank_deficient=nvars>rank,
                 sample_nonzero_eigen_min=low,sample_eigen_max=top,
                 shrinkage_alpha=alpha,estimated_lambda_min=min_eig,
                 estimated_lambda_max=max_eig,estimated_condition_number=cond if cond is not None else "",
                 interpretation="EUR_REFERENCE_ONLY_NOT_QTL_COHORT_OR_SIGNAL_BLOCK"))
        global_summary.append(dict(locus=loc,variant_keys=len(ids),
          monomorphic=int(np.sum(~variant)),contiguous_250kb_blocks=len(tiles),
          blocks_with_2plus_polymorphic=good,largest_block_polymorphic_snps=largest,
          max_locus_full_ld_rank=502,GTEx_donor_LD_verified=False,
          interpretation="NUMERICAL_QC_ONLY_NO_FORMAL_SUSIE_COLOC"))
    with (out/"IS_PRIORITY4_EUR503_250KB_BLOCK_SHRINKAGE_QC.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
    with (out/"IS_PRIORITY4_EUR503_BLOCK_QC_SUMMARY.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(global_summary[0]),delimiter="\t");w.writeheader();w.writerows(global_summary)
    summary=dict(status="REFERENCE_LD_SHRINKAGE_ONLY",blocks=global_summary,
      shrinkage_grid=list(SHRINKAGE),block_width_bp=BLOCK_BP,
      no_cohort_matched_QTL_LD=True,no_SuSiE_coloc=True,no_causal_claim=True,
      original_candidate_genes=2225,
      raw_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in source.glob("*.raw")})
    (out/"IS_PRIORITY4_EUR503_BLOCK_SHRINKAGE_MANIFEST.json").write_text(json.dumps(summary,indent=2)+"\n")
    return summary

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(analyze(a.input,a.out),indent=2))
