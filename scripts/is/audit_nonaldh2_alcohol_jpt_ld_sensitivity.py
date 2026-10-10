#!/usr/bin/env python3
"""Japanese JPT (1000G n104) vs EAS n504 LD sensitivity on seven alcohol SNPs.

Uses verified already downloaded 1KG phase3 v5b genotype dosages;
does not redownload or infer any clinical causality. Monomorphic rarer
variants in small JPT are flagged, never replaced with artificial r2=0.
"""
import argparse,csv,json,math
from pathlib import Path
import numpy as np
LD=Path("/srv/is-analysis/data/is/ld_reference/nonaldh2_alcohol_eas_20261010")
PANEL=Path("/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/integrated_call_samples_v3.20130502.ALL.panel")
REGIONS={"2":["2:27730940:T:C"],
         "4":["4:39413780:A:G","4:100239319:T:C"],
         "9":["9:38395928:T:C","9:75461066:T:C"],
         "12":["12:106750302:A:G","12:112241766:G:A"]}
def run(root,panel):
    with panel.open() as f:
        all_samples=list(csv.DictReader(f,delimiter="\t"))
    jap={r["sample"] for r in all_samples if r["pop"]=="JPT" and r["super_pop"]=="EAS"}
    if len(jap)!=104:raise ValueError("Expected reference JPT n104")
    g={}
    sample_ref=None
    for c,vars_ in REGIONS.items():
        with (root/f"chr{c}.7snp_eas.dosages.tsv").open() as f:
            rows=list(csv.DictReader(f,delimiter="\t"))
        ids=[x["EAS_sample_ID"] for x in rows]
        if len(ids)!=504 or (sample_ref is not None and ids!=sample_ref):
            raise ValueError("EAS sample mismatch across VCF chromosomes")
        sample_ref=ids
        for v in vars_:
            g[v]=np.asarray([float(r[v]) for r in rows])
    selected=np.asarray([a in jap for a in sample_ref])
    if selected.sum()!=104:raise ValueError("JPT subset not represented completely")
    ids=[x for xx in REGIONS.values() for x in xx]
    freq=[];out=[]
    for v in ids:
        e=g[v];j=e[selected]
        if not(np.isfinite(e).all() and np.isfinite(j).all()):
            raise ValueError("Unexpected missing genotype")
        maf_e=float(np.mean(e)/2);maf_j=float(np.mean(j)/2)
        freq.append({"SNP":v,"EAS_ALT_EAF":maf_e,"JPT_ALT_EAF":maf_j,
            "JPT_minor_allele_count":int(min(sum(j),2*len(j)-sum(j))),
            "JPT_genotype_variance_positive":int(np.std(j)>0),
            "JPT_n":104,"EAS_n":504})
    for i,va in enumerate(ids):
        for vb in ids[i+1:]:
            easx,easy=g[va],g[vb];jpx,jpy=easx[selected],easy[selected]
            r_e=float(np.corrcoef(easx,easy)[0,1])
            v1=float(np.std(jpx));v2=float(np.std(jpy))
            if v1<=1e-14 or v2<=1e-14:
                jpt_r2="";status="UNDEFINED_JPT_MONOMORPHIC_OR_ZERO_VARIANCE"
            else:
                r_j=float(np.corrcoef(jpx,jpy)[0,1]);jpt_r2=r_j*r_j
                status="JPT_ESTIMATED_104_SAMPLES_RARE_VARIANT_UNCERTAINTY"
            out.append({"SNP1":va,"SNP2":vb,
              "chr1":va.split(":")[0],"chr2":vb.split(":")[0],
              "same_chr":int(va.split(":")[0]==vb.split(":")[0]),
              "EAS_504_ALT_dosage_r2":r_e*r_e,
              "JPT_104_ALT_dosage_r2":jpt_r2,
              "JPT_r2_status":status,
              "JPT_reference_n":104,"full_EAS_n":504,
              "source_genotype_ref_alt_aligned":True,
              "LD_independence_does_not_prove_pleiotropy_absent":True})
    with (root/"IS_RS671_7SNP_EAS504_JPT104_LD_SENSITIVITY.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter="\t")
        w.writeheader();w.writerows(out)
    with (root/"IS_RS671_7SNP_EAS_JPT_ALLELE_FREQUENCIES.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(freq[0]),delimiter="\t")
        w.writeheader();w.writerows(freq)
    measured=[x for x in out if isinstance(x["JPT_104_ALT_dosage_r2"],float)]
    report={"original_ref":"1000G_phase3_EAS504_JPT104",
       "sample_count_EAS":504,"sample_count_JPT":104,
       "real_variant_count":len(ids),
       "pair_count":len(out),
       "JPT_pairs_estimable":len(measured),
       "JPT_pairs_undefined_zero_genotype_variance":len(out)-len(measured),
       "max_measurable_JPT_r2":max(float(x["JPT_104_ALT_dosage_r2"]) for x in measured) if measured else None,
       "max_full_EAS_r2":max(float(x["EAS_504_ALT_dosage_r2"]) for x in out),
       "JPT_frequencies":{x["SNP"]:x["JPT_ALT_EAF"] for x in freq},
       "JPT_r2_sample_limited_no_causal_MR_gate":True,
       "source_JPT_pop_IDs_from_official_phase3_panel":True,
       "status":"JPT_REFERENCE_SENSITIVITY_ONLY_UNCERTAIN_RARE_MAF"}
    (root/"IS_RS671_7SNP_JPT_SENSITIVITY_SUMMARY.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    return report
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=LD)
    p.add_argument("--panel",type=Path,default=PANEL)
    a=p.parse_args();run(a.root,a.panel)
