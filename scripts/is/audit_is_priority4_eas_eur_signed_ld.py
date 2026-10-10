#!/usr/bin/env python3
"""SNP-level 1000G EAS504 vs EUR503 allele-aligned reference signed-LD diagnostic.

No GTEx cohort LD, no formal coloc, no inference of gene causality.
All PLINK --export A counted alleles are converted to GRCh37 ALT dosages.
"""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

GENES={
 "BBJ_IS_L001":[("FGF5","ENSG00000138675"),("C4orf22","ENSG00000197826")],
 "BBJ_IS_L002":[("CALHM2","ENSG00000138172"),("NEURL","ENSG00000107954")]
}
TISSUE="GTEx_V8__Brain_Cerebellar_Hemisphere"

def alt_dosage(values, counted, ref, alt):
    if counted==alt:return values
    if counted==ref:return 2-values
    raise ValueError(f"Counted allele {counted} neither REF {ref} nor ALT {alt}")

def signed_r(a,b,min_pairs=100):
    mask=np.isfinite(a)&np.isfinite(b)
    n=int(mask.sum())
    if n<min_pairs:return float("nan"),n
    x=a[mask];y=b[mask]
    dx=x-x.mean();dy=y-y.mean()
    den=float(np.sqrt((dx@dx)*(dy@dy)))
    if den<=0:return float("nan"),n
    return float((dx@dy)/den),n

def load_raw(path,expected):
    df=pd.read_csv(path,sep="\t")
    if len(df)!=expected:raise ValueError(f"{path} sample count {len(df)} != {expected}")
    matrix=df.iloc[:,6:].to_numpy(dtype=float,copy=True)
    ids=[]
    for j,col in enumerate(df.columns[6:]):
        if "_" not in col:raise ValueError("Missing counted allele suffix")
        snp,counted=col.rsplit("_",1)
        alleles=snp.split(":")
        if len(alleles)!=4:raise ValueError(f"Not a chrom:pos:ref:alt key {snp}")
        matrix[:,j]=alt_dosage(matrix[:,j],counted,alleles[2],alleles[3])
        ids.append(snp)
    if len(ids)!=len(set(ids)):raise ValueError("Duplicated SNP in PLINK raw")
    if np.any(np.isfinite(matrix)&((matrix<0)|(matrix>2))):raise ValueError("Invalid hardcall dosage")
    return {v:matrix[:,i] for i,v in enumerate(ids)}

def tsv_write(path,records):
    if not records:raise ValueError("No rows to write")
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t")
        w.writeheader();w.writerows(records)

def simple_float(v):
    try:return float(v)
    except (TypeError,ValueError):return float("nan")

def run(rawdir,source,out):
    out.mkdir(parents=True,exist_ok=True)
    details=[]; summaryrows=[];coverage=[]
    for locus,genes in GENES.items():
        short=locus[-4:] # L001
        eas=load_raw(rawdir/f"{short}.EAS504.raw",504)
        eur=load_raw(rawdir/f"{short}.EUR503.raw",503)
        ref=list(eas)
        if set(eas)!=set(eur):
            raise ValueError(f"Ancestry SNP set mismatch {locus}: {len(eas)} {len(eur)}")
        if len(ref)!=({"BBJ_IS_L001":1694,"BBJ_IS_L002":1583}[locus]):
            raise ValueError("Unexpected locus candidate SNP set")
        original={}
        for symbol,gene in genes:
            sp=source/f"{locus}__{TISSUE}__{gene}.tsv"
            with sp.open(newline="") as f:rows=list(csv.DictReader(f,delimiter="\t"))
            if len(rows)!=len(ref):raise ValueError("Original assay variant count changed")
            by_id={r["variant_id"]:r for r in rows}
            if set(by_id)!=set(eas):raise ValueError("Original molecular SNP IDs mismatch reference")
            topg=min(rows,key=lambda r:float(r["gwas_p"]))
            topq=min(rows,key=lambda r:float(r["eqtl_p"]))
            original[symbol]={"topg":topg,"topq":topq,"rows":by_id}
        topg=original[genes[0][0]]["topg"]["variant_id"]
        for symbol,_ in genes[1:]:
            if original[symbol]["topg"]["variant_id"]!=topg:
                raise ValueError("GWAS lead differs within same locus")
        frequencies={}
        for var in ref:
            ae=eas[var];au=eur[var]
            f_e=float(np.nanmean(ae)/2);f_u=float(np.nanmean(au)/2)
            if not 0<=f_e<=1 or not 0<=f_u<=1:raise ValueError("Nonphysical ALT frequency")
            frequencies[var]={"EAS_ALT_EAF":f_e,"EUR_ALT_EAF":f_u,
                "EAS_MAF":min(f_e,1-f_e),"EUR_MAF":min(f_u,1-f_u),
                "EAS_missingness":float(np.isnan(ae).mean()),
                "EUR_missingness":float(np.isnan(au).mean())}
        for var in ref:
            re,ne=signed_r(eas[topg],eas[var])
            ru,nu=signed_r(eur[topg],eur[var])
            x=frequencies[var]
            details.append(dict(locus=locus,variant_b37=var,lead_GWAS_SNP_b37=topg,
                EAS_ALT_EAF=x["EAS_ALT_EAF"],EUR_ALT_EAF=x["EUR_ALT_EAF"],
                EAS_MAF=x["EAS_MAF"],EUR_MAF=x["EUR_MAF"],
                EAS_missingness=x["EAS_missingness"],EUR_missingness=x["EUR_missingness"],
                EAS_signed_r_with_GWAS_lead=re,EUR_signed_r_with_GWAS_lead=ru,
                EAS_r2_with_GWAS_lead=re*re,EUR_r2_with_GWAS_lead=ru*ru,
                EAS_pairwise_N=ne,EUR_pairwise_N=nu,
                abs_delta_r=(abs(re-ru) if np.isfinite(re) and np.isfinite(ru) else ""),
                abs_delta_r2=(abs(re*re-ru*ru) if np.isfinite(re) and np.isfinite(ru) else ""),
                interpretation="REFERENCE_LD_ONLY_NOT_STUDY_MATCHED"))
        for symbol,gene in genes:
            gw=original[symbol]["topg"];q=original[symbol]["topq"];v=q["variant_id"];rg=topg
            if rg not in eas or v not in eas:raise ValueError("Lead SNP absent")
            re,ne=signed_r(eas[rg],eas[v]);ru,nu=signed_r(eur[rg],eur[v])
            gr=original[symbol]["rows"]
            maf_eas_g=np.array([abs(frequencies[k]["EAS_MAF"]-float(gr[k]["gwas_maf"])) for k in ref])
            maf_eur_q=np.array([abs(frequencies[k]["EUR_MAF"]-float(gr[k]["eqtl_maf"])) for k in ref])
            sumrow=dict(symbol=symbol,gene_id=gene,locus=locus,n_locus_snps=len(ref),
                gwas_lead_b37=rg,gwas_lead_gwas_p=float(gw["gwas_p"]),
                qtl_lead_b37=v,qtl_lead_p=float(q["eqtl_p"]),
                EAS_reference_signed_r=re,EUR_reference_signed_r=ru,
                EAS_reference_r2=re*re,EUR_reference_r2=ru*ru,
                EAS_reference_n=ne,EUR_reference_n=nu,
                abs_delta_signed_r=(abs(re-ru) if np.isfinite(re) and np.isfinite(ru) else ""),
                abs_delta_r2=(abs(re*re-ru*ru) if np.isfinite(re) and np.isfinite(ru) else ""),
                EAS_vs_BBJ_GWAS_median_abs_MAF_difference=float(np.median(maf_eas_g)),
                EUR_vs_GTEx_QTL_median_abs_MAF_difference=float(np.median(maf_eur_q)),
                EAS_vs_BBJ_GWAS_frac_MAF_diff_gt_0p1=float(np.mean(maf_eas_g>.1)),
                EUR_vs_GTEx_QTL_frac_MAF_diff_gt_0p1=float(np.mean(maf_eur_q>.1)),
                claim_status="CROSS_ANCESTRY_REFERENCE_LD_TAGGING_DIAGNOSTIC_NOT_COLOC")
            summaryrows.append(sumrow)
        both=[x for x in details if x["locus"]==locus and x["abs_delta_r"]!=""]
        coverage.append(dict(locus=locus,n_snps=len(ref),n_with_informative_both_ancestry_r=len(both),
            informative_fraction=len(both)/len(ref),
            max_EAS_missingness=max(frequencies[v]["EAS_missingness"] for v in ref),
            max_EUR_missingness=max(frequencies[v]["EUR_missingness"] for v in ref),
            median_abs_EAS_EUR_ALT_EAF_diff=float(np.median([abs(frequencies[v]["EAS_ALT_EAF"]-frequencies[v]["EUR_ALT_EAF"]) for v in ref])),
            frac_informative_snps_with_abs_delta_r2_gt_0p2=float(np.mean([x["abs_delta_r2"]>.2 for x in both])) if both else None))
    tsv_write(out/"IS_PRIORITY4_EUR_EAS_GWAS_LEAD_TAGGING_BY_SNP.tsv",details)
    tsv_write(out/"IS_PRIORITY4_EUR_EAS_GENE_MAIN_SIGNED_LD.tsv",summaryrows)
    tsv_write(out/"IS_PRIORITY4_EUR_EAS_LOCI_QC.tsv",coverage)
    manifest=dict(status="EAS504_EUR503_REFERENCE_SIGNED_LD_ONLY",gene_count=len(summaryrows),
                  loci=coverage,top_genes=summaryrows,original_candidate_universe=2225,
                  no_gtEx_cohort_LD=True,no_coloc_susie_executed=True,
                  source_ancestry_counts={"EAS":504,"EUR":503},
                  source_raw_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(rawdir.glob("L00*.raw"))},
                  warning="Signed LD is calculated on original GRCh37 ALT allele direction. EAS504 and EUR503 are independent 1000G ancestry strata; neither panel equals BBJ GWAS or GTEx QTL in-study LD. Locus lead selection uses minimum p within available archived coloc SNP intersection; not necessarily genome-wide or full tested-set primary lead.")
    (out/"IS_PRIORITY4_EUR_EAS_SIGNED_LD_SUMMARY.json").write_text(json.dumps(manifest,indent=2)+"\n")
    return manifest

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for key in ("rawdir","source","out"):p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    result=run(a.rawdir,a.source,a.out)
    for x in result["top_genes"]:
        print(x["symbol"],"EAS signed_r",round(x["EAS_reference_signed_r"],4),"EUR signed_r",round(x["EUR_reference_signed_r"],4),"QTL P",x["qtl_lead_p"])
    print("QC",result["loci"])
