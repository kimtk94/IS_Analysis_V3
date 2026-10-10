#!/usr/bin/env python3
"""Source-verified Japanese JPT104 LD sensitivity for alcohol ADH1B/ALDH2.

Read the SAME genotype export as the established EAS504 LD; match all sample
identifiers to official 1000G Phase3 panel; verify independently reconstructed
EAS signed dosage LD agrees with existing EAS matrix to 1e-9 before writing
JPT output. Do not infer that small JPT n104 solves reference LD mismatch.
"""
import argparse,csv,json,hashlib
from pathlib import Path
import numpy as np
BASE=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_conditional_susie_v1")
PANEL=Path("/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/integrated_call_samples_v3.20130502.ALL.panel")
SAMPLES=Path("/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/EAS.samples.txt")
LOCI=("ADH1B","ALDH2")
def read_rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def build(base=BASE,panel=PANEL,samples=SAMPLES,loci=LOCI):
    pa=read_rows(panel)
    jpt={r["sample"] for r in pa if r["pop"]=="JPT" and r["super_pop"]=="EAS"}
    expected=[x.strip() for x in samples.read_text().splitlines() if x.strip()]
    if len(jpt)!=104 or len(expected)!=504 or len(set(expected))!=504 or not jpt.issubset(set(expected)):
        raise ValueError("Official 1000G EAS/JPT reference panel mismatch")
    results={}
    for locus in loci:
        folder=base/locus
        state=json.loads((folder/"input_qc.json").read_text())
        if state["status"]!="ALCOHOL_REGIONAL_REAL_GWAS_EAS_LD_INPUT_PASS" or state["reference_EAS_n"]!=504:
            raise ValueError("Unverified EAS source genotype for "+locus)
        variants=read_rows(folder/"variants.tsv")
        byid={v["ID"]:v for v in variants}
        if len(byid)!=len(variants) or not 50<=len(variants)<=1600:raise ValueError("Source SNP count or ID duplication")
        data=[]
        observed=[]
        with (folder/"reference_eas504.traw").open() as f:
            r=csv.reader(f,delimiter="\t")
            header=next(r)
            if header[:6]!=["CHR","SNP","(C)M","POS","COUNTED","ALT"]:
                raise ValueError("Unexpected export Av schema")
            ids=[i.split("_")[0] for i in header[6:]]
            if len(ids)!=504 or len(set(ids))!=504 or set(ids)!=set(expected):
                raise ValueError("PLINK source individual EAS identities not complete")
            mask=np.array([a in jpt for a in ids])
            if int(mask.sum())!=104:raise ValueError("JPT subset incomplete")
            for row in r:
                if len(row)!=6+504:raise ValueError("Unexpected genotype row width")
                chrom,var,cm,pos,counted,other=row[:6]
                q=byid.get(var)
                if q is None:raise ValueError("Unexpected variant in source dosage")
                ref,alt=var.split(":")[2:]
                if {counted,other}!={ref,alt} or q["reference_ALT"]!=alt:
                    raise ValueError("Reference ALT dosage orientation mismatch for "+var)
                x=np.array([float(y) for y in row[6:]],dtype=np.float64)
                if not np.isfinite(x).all() or not np.isin(x,[0,1,2]).all():
                    raise ValueError("Unexpected missing/non-diploid 1KG genotype "+var)
                if counted!=alt:x=2-x
                if abs(float(x.mean()/2)-float(q["reference_EAS_ALT_EAF"]))>1e-9:
                    raise ValueError("Source EAS frequency from traw does not match preparer")
                data.append(x)
                observed.append(var)
        if len(data)!=len(byid) or set(observed)!=set(byid):
            raise ValueError("SNP coverages do not match original source")
        # SOURCE LD is in variants.tsv order, not the PLINK traw order.
        order={vid:i for i,vid in enumerate(observed)}
        A=np.stack([data[order[v["ID"]]] for v in variants],axis=0)
        def corr_genotypes(g):
            if any((g.std(axis=1,ddof=1)<=1e-9)):
                raise ValueError("Reference subset contains monomorphic SNP")
            x=(g-g.mean(axis=1,keepdims=True))/g.std(axis=1,ddof=1,keepdims=True)
            return (x@x.T)/(g.shape[1]-1)
        eas=corr_genotypes(A)
        p=len(variants)
        saved=np.fromfile(folder/"ld.f64.rowmajor",dtype="<f8").reshape((p,p))
        delta=float(np.max(np.abs(saved-eas)))
        if delta>1e-9:raise ValueError("Existing EAS LD cannot be reconstructed: matrix/individual/sample allele order mismatch")
        jap=A[:,mask]
        maf=np.mean(jap,axis=1)/2
        monomorphic=[variants[i]["ID"] for i,m in enumerate(maf) if m<=0 or m>=1]
        keep=[i for i,m in enumerate(maf) if m>0 and m<1]
        rj=corr_genotypes(jap[keep])
        # Store only polymorphic JPT variants; never force undefined LD to r2=0.
        output=folder/"JPT104_sensitivity"
        output.mkdir(exist_ok=True)
        selected=[variants[i] for i in keep]
        with (output/"variants.tsv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(selected[0])+["JPT_ALT_freq","JPT_minor_allele_count"],
              delimiter="\t")
            w.writeheader()
            for i,v in zip(keep,selected):
                w.writerow({**v,"JPT_ALT_freq":float(maf[i]),
                  "JPT_minor_allele_count":int(min(jap[i].sum(),208-jap[i].sum()))})
        (output/"ld.f64.rowmajor").write_bytes(rj.astype("<f8").tobytes())
        info={"locus":locus,"status":"JPT104_LD_SOURCE_RECONSTRUCTED_EAS_BASELINE_CROSSCHECK_PASS",
           "1000G_reference":"GRCh37_20130502_Phase3_EAS504_JPT104",
           "source_original_EAS_n":504,"JPT_n":104,"source_snps":p,
           "JPT_polymorphic_snps":len(keep),"JPT_monomorphic_snps":len(monomorphic),
           "JPT_monomorphic_variant_ids":monomorphic,
           "max_reconstructed_existing_EAS_signed_LD_absolute_difference":delta,
           "EAS_versus_JPT_af_mean_difference":float(np.mean(np.abs(A.mean(axis=1)/2-maf))),
           "max_JPT_vs_original_Japanese_GWAS_ALT_af_difference":
              float(max(abs(float(v["Japanese_original_ALT_EAF"])-float(maf[i])) for i,v in enumerate(variants))),
           "JPT_reference_matrix_rank_at_most":103,
           "EAS_reference_matrix_rank_at_most":503,
           "not_study_matched_GWAS_LD":True,
           "exclusion_restriction_untested":True,
           "no_causal_finemap_from_JPT_alone":True}
        (output/"input_qc.json").write_text(json.dumps(info,indent=2))
        print(json.dumps(info,indent=2),flush=True)
        results[locus]=info
    return results
if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--root",type=Path,default=BASE)
    a.add_argument("--panel",type=Path,default=PANEL)
    a.add_argument("--samples",type=Path,default=SAMPLES)
    opt=a.parse_args();build(opt.root,opt.panel,opt.samples)
