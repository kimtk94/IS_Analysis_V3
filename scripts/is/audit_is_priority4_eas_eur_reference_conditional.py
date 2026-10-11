#!/usr/bin/env python3
"""Matched EUR503 vs EAS504 1000G ALT-dose LD and BBJ lead-residual sensitivity.

This is reference-LD proxy screening only. BBJ GWAS is East Asian, and a
heuristic lead-residual z is NOT a formally calibrated conditional GWAS test.
No credible sets, SuSiE, coloc, independence or causal effects are inferred.
"""
import argparse,csv,hashlib,json,math
from pathlib import Path
import numpy as np
import pandas as pd

LOCI={"BBJ_IS_L001":("4:81182554:T:C",2834),
      "BBJ_IS_L002":("10:105459834:T:C",2828)}
CANDIDATE="10:105694487:G:T"
ALPHA=(0.0,0.01,0.05,0.10,0.20)
Z_SCREEN=5.45

def read_source(src):
    with src.open(newline="") as f: rr=list(csv.DictReader(f,delimiter="\t"))
    if len(rr)!=11050:raise ValueError("Frozen QTL × BBJ gene-SNP source changed")
    locus={}
    for r in rr:
        k=(r["locus"],r["variant_GRCh37_REF_ALT"])
        if k in locus:
            x=locus[k]
            for field in ("BBJ_beta_ALT","BBJ_se","BBJ_p","BBJ_EAF_ALT"):
                if r[field]!=x[field]:raise ValueError("BBJ source inconsistent between gene assays")
        locus[k]=r
    return locus

def extract_alt(path, expected_n, expected_count):
    frame=pd.read_csv(path,sep="\t")
    if len(frame)!=expected_n:raise ValueError("1000G population sample denominator mismatch")
    x=frame.iloc[:,6:].to_numpy(dtype=float,copy=True)
    identifiers=[]
    for i,col in enumerate(frame.columns[6:]):
        k,counted=col.rsplit("_",1)
        p=k.split(":")
        if len(p)!=4 or not p[1].isdigit():raise ValueError("Invalid VCF/PGEN genomic key")
        ref,alt=p[2:]
        if counted==ref:x[:,i]=2-x[:,i]
        elif counted!=alt:raise ValueError("Counted allele not genomic REF or ALT")
        identifiers.append(k)
    if len(identifiers)!=expected_count or len(set(identifiers))!=len(identifiers):
        raise ValueError("Matched 1000G reference SNP set drift")
    if not np.isfinite(x).all():raise ValueError("Unexpected missing 1000G genotypes")
    if np.any((x<0)|(x>2)):raise ValueError("Invalid phased hardcall ALT dose")
    af=x.mean(axis=0)/2.0
    centered=x-x.mean(axis=0)
    length=np.sqrt(np.sum(centered*centered,axis=0))
    return identifiers,centered,length,af

def lead_correlations(identifiers,x,length,lead):
    j=identifiers.index(lead)
    denom=length*length[j]
    r=np.full(len(identifiers),np.nan)
    ok=denom>0
    r[ok]=(x[:,j]@x[:,ok])/denom[ok]
    r[ok]=np.clip(r[ok],-1.0,1.0)
    if abs(r[j]-1.0)>1e-10:raise ValueError("Lead correlation not unity")
    return r

def write_table(path,rows):
    if not rows:raise ValueError("Empty output")
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader();w.writerows(rows)

def run(source,eas_dir,eur_dir,out):
    if out.exists() and any(out.iterdir()):raise FileExistsError("Refuse overwrite of research results")
    src=read_source(source)
    table=[];summary=[];candidates=[];scenarios=[]
    for locus,(lead,expected) in LOCI.items():
        files={
           "EAS504":(eas_dir/f"{locus}.EAS504.fullcis.raw",504),
           "EUR503":(eur_dir/f"{locus}.EUR503.raw",503)}
        cohorts={}
        for name,(path,n) in files.items():
            keys,x,l,af=extract_alt(path,n,expected)
            cohorts[name]=(keys,lead_correlations(keys,x,l,lead),af,l)
        a,ra,afa,la=cohorts["EAS504"]
        b,rb,afb,lb=cohorts["EUR503"]
        if set(a)!=set(b):raise ValueError("Ancestry reference keys differ")
        if lead not in a or (locus,lead) not in src:raise ValueError("Missing fixed BBJ lead association")
        j_a={k:i for i,k in enumerate(a)}
        j_b={k:i for i,k in enumerate(b)}
        lead_z=float(src[(locus,lead)]["BBJ_beta_ALT"])/float(src[(locus,lead)]["BBJ_se"])
        if not math.isfinite(lead_z):raise ValueError("Invalid BBJ lead z")
        valid_comparable=0;monomorphic_EAS=0;monomorphic_EUR=0
        sign_discord=0;material_r2=0
        screen_counts={(anc,alpha):0 for anc in ("EAS504","EUR503") for alpha in ALPHA}
        eligible_counts={(anc,alpha):0 for anc in ("EAS504","EUR503") for alpha in ALPHA}
        eligible={k:src[(locus,k)] for k in a if (locus,k) in src}
        if len(eligible)!=expected:raise ValueError("ABF SNP/source matching incomplete")
        for snp in sorted(eligible):
            i=j_a[snp];j=j_b[snp]
            r_e=float(ra[i]);r_u=float(rb[j])
            poly_e=math.isfinite(r_e);poly_u=math.isfinite(r_u)
            if not poly_e:monomorphic_EAS+=1
            if not poly_u:monomorphic_EUR+=1
            both=poly_e and poly_u
            z=float(eligible[snp]["BBJ_beta_ALT"])/float(eligible[snp]["BBJ_se"])
            if not math.isfinite(z):raise ValueError("Invalid BBJ SNP z")
            if both:
                valid_comparable+=1
                if abs(r_e)>=.1 and abs(r_u)>=.1 and r_e*r_u<0:sign_discord+=1
                if abs(r_e*r_e-r_u*r_u)>.2:material_r2+=1
            table.append(dict(locus=locus,variant_GRCh37=snp,
                BBJ_marginal_z=z,lead_variant_GRCh37=lead,
                EAS504_signed_r=r_e if poly_e else "",
                EUR503_signed_r=r_u if poly_u else "",
                EAS504_r2=r_e*r_e if poly_e else "",
                EUR503_r2=r_u*r_u if poly_u else "",
                EAS504_ALT_EAF=float(afa[i]),EUR503_ALT_EAF=float(afb[j]),
                signed_r_opposite_if_abs_r_both_ge_0p1=both and abs(r_e)>=.1 and abs(r_u)>=.1 and r_e*r_u<0,
                delta_r2_EAS_minus_EUR=(r_e*r_e-r_u*r_u) if both else "",
                both_informative=both,
                ancestry_reference_only_not_BBJ_or_GTEx_study_LD=True))
            if not both or snp==lead:
                continue
            for alpha in ALPHA:
                for anc,rr in (("EAS504",r_e),("EUR503",r_u)):
                    q=(1-alpha)*rr
                    denom=1-q*q
                    if denom<=1e-9:continue
                    zz=(z-q*lead_z)/math.sqrt(denom)
                    eligible_counts[(anc,alpha)]+=1
                    screen_counts[(anc,alpha)]+=abs(zz)>=Z_SCREEN
                    if snp==CANDIDATE:
                        candidates.append(dict(locus=locus,variant_GRCh37=snp,
                            ancestry=anc,n_reference=504 if anc=="EAS504" else 503,
                            lead_GRCh37=lead,BBJ_lead_z=lead_z,
                            BBJ_candidate_marginal_z=z,
                            lead_candidate_signed_r=rr,lead_candidate_r2=rr*rr,
                            reference_ALT_EAF= float(afa[i] if anc=="EAS504" else afb[j]),
                            alpha=alpha,approx_lead_residual_z=zz,
                            absZ_gt_5p45=abs(zz)>=Z_SCREEN,
                            provenance="REFERENCE_LD_HEURISTIC_NOT_CONDITIONAL_GWAS"))
        for alpha in ALPHA:
            ea=screen_counts["EAS504",alpha];eu=screen_counts["EUR503",alpha]
            ne=eligible_counts["EAS504",alpha];nu=eligible_counts["EUR503",alpha]
            scenarios.append(dict(locus=locus,alpha=alpha,
                EAS504_evaluable=ne,EUR503_evaluable=nu,
                EAS504_approx_residual_absZ_gt_5p45=ea,
                EUR503_approx_residual_absZ_gt_5p45=eu,
                threshold=Z_SCREEN,
                not_formal_conditional_GWAS=True))
        summary.append(dict(locus=locus,lead_b37=lead,
            n_EAS504=504,n_EUR503=503,n_exact_shared_SNPs=len(a),
            n_allele_match_to_BBJ_and_QTD=len(eligible),
            monomorphic_EAS=monomorphic_EAS,monomorphic_EUR=monomorphic_EUR,
            both_ancestry_lead_r_informative=valid_comparable,
            direction_opposite_abs_r_ge_0p1_both=sign_discord,
            delta_r2_abs_gt_0p2=material_r2,
            LD_cohort_specific_to_GTEx_or_BBJ=False,
            formal_GWAS_conditional_or_SuSiE_performed=False))
    if len(candidates)!=10:raise ValueError("Expected 5 alpha x 2 ancestries for one candidate")
    out.mkdir(parents=True)
    write_table(out/"IS_PRIORITY4_EAS_EUR_COMMON_SNP_LEAD_LD.tsv",table)
    write_table(out/"IS_PRIORITY4_EAS_EUR_LEAD_LD_SUMMARY.tsv",summary)
    write_table(out/"IS_CHR10_105694487_EAS_EUR_CANDIDATE_SENSITIVITY.tsv",candidates)
    write_table(out/"IS_PRIORITY4_EAS_EUR_APPROX_CONDITIONAL_SCREEN.tsv",scenarios)
    manifest=dict(schema="IS_PRIORITY4_MATCHED_ANCESTRY_REFERENCE_CONDITIONAL_V1",
        n_loci=2,n_EAS=504,n_EUR=503,alpha_grid=list(ALPHA),
        candidate=CANDIDATE,reference_only=True,
        status="ANCESTRY_LD_DIAGNOSTIC_NOT_BBJ_OR_GTEX_COHORT_MATCHED",
        approximate_z_not_calibrated_conditioning=True,
        all_2225_genes_preserved=True,original_646_ABF_not_modified=True,
        newer_120_ABF_not_modified=True,
        summaries=summary,
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in (source,*[eas_dir/f"{l}.EAS504.fullcis.raw" for l in LOCI],
                             *[eur_dir/f"{l}.EUR503.raw" for l in LOCI])})
    (out/"IS_PRIORITY4_EAS_EUR_REFERENCE_CONDITIONAL_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
    return manifest

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for arg in ("source","eas-dir","eur-dir","out"):
        p.add_argument("--"+arg,type=Path,required=True)
    a=p.parse_args()
    m=run(a.source,a.eas_dir,a.eur_dir,a.out)
    print(json.dumps({k:v for k,v in m.items() if k!="source_sha256"},indent=2))
