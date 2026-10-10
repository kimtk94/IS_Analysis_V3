#!/usr/bin/env python3
"""Audit G0022 OASIS gene-cell cis peaks relative to ALDH2 rs671 in EAS/JPT.

Input is public 1000G EAS504 PLINK-exported genotype .traw, not OASIS LD.
No full-cis completeness, beta harmonization, colocalization, MR or independent
causal gene evidence is inferred. JPT104 is a subset of EAS504 (not independent).
"""
import argparse,csv,hashlib,json,math
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np

ANCHOR="12:112241766:G:A"
PROXY="12:112168009:G:A"
EXPECTED_PIP_CS={
 "12:112241766:G:A":("rs671",1.0),
 "12:112168009:G:A":("rs11066015",.987528196),
 "12:112468206:C:T":("rs11066132",.97423724),
 "12:112736118:A:G":("rs77768175",.869614381),
}

def panel_mapping(path, samples):
    x={}
    with path.open() as f:
        for i,row in enumerate(csv.DictReader(f,delimiter="\t")):
            sample=row.get("sample")
            pop=row.get("pop")
            if sample in x:
                raise ValueError("Duplicate 1KG panel sample")
            x[sample]=pop
    got=[str(sample).split("_",1)[0] for sample in samples]
    if len(got)!=504 or len(set(got))!=504:
        raise ValueError("Not 504 unique EAS panel people")
    missing=[sample for sample in got if sample not in x]
    if missing:
        raise ValueError("Samples missing original reference panel")
    jpt=np.asarray([i for i,s in enumerate(got) if x[s]=="JPT"],dtype=int)
    if len(jpt)!=104:
        raise ValueError("JPT sample count differs from verified source n104")
    pops=Counter(x[s] for s in got)
    if set(pops)!={"JPT","CHB","CHS","CDX","KHV"}:
        raise ValueError("Reference population composition changed")
    return jpt,dict(pops)

def read_ld_genotypes(traw,panel, requested):
    with traw.open(newline="") as f:
        r=csv.reader(f,delimiter="\t")
        head=next(r)
        if head[:6]!=["CHR","SNP","(C)M","POS","COUNTED","ALT"]:
            raise ValueError("Unexpected PLINK dosage schema")
        jpt,pops=panel_mapping(panel,head[6:])
        rows={}
        for row in r:
            if len(row)!=len(head):
                raise ValueError("Truncated PLINK source row")
            ident=row[1]
            if ident not in requested and ident not in EXPECTED_PIP_CS:
                continue
            if ident in rows:
                raise ValueError("Duplicate genotype SNP")
            chrom,pos,ref,alt=ident.split(":")
            if row[0]!=chrom or row[3]!=pos:
                raise ValueError("CHR POS genotype key disagreement")
            if (row[4],row[5]) not in ((ref,alt),(alt,ref)):
                raise ValueError("PLINK COUNTED/ALT inconsistent with variant key")
            dos=np.array([np.nan if v in ("NA",".","") else float(v) for v in row[6:]],dtype=float)
            if np.any(np.isfinite(dos)&((dos<0)|(dos>2))):
                raise ValueError("Genotype dosage outside [0,2]")
            # PLINK Av counts COUNTED allele; make all effects ALT-positive
            alt_dos=2-dos if row[4]==ref else dos
            rows[ident]=alt_dos
    if ANCHOR not in rows or PROXY not in rows:
        raise ValueError("Required four-locus LD markers unavailable")
    if not set(EXPECTED_PIP_CS).issubset(rows):
        raise ValueError("Original four GWAS-CS reference genotypes absent")
    return rows,jpt,pops

def pairwise_ld(anchor,variant,sample_idx=None,min_n=80):
    if sample_idx is not None:
        a=anchor[sample_idx];v=variant[sample_idx]
    else:
        a=anchor;v=variant
    mask=np.isfinite(a)&np.isfinite(v)
    n=int(mask.sum())
    if n<min_n:
        return {"n_pair":n,"r":None,"r2":None}
    a=a[mask];v=v[mask]
    sa=np.std(a,ddof=1);sv=np.std(v,ddof=1)
    if sa<1e-8 or sv<1e-8:
        return {"n_pair":n,"r":None,"r2":None}
    corr=float(np.corrcoef(a,v)[0,1])
    return {"n_pair":n,"r":corr,"r2":corr*corr}

def read_qtl_pairs(path):
    by=defaultdict(list)
    with path.open(newline="") as f:
        for x in csv.DictReader(f,delimiter="\t"):
            if x["source_ref_alt_id_match"]!="EXACT_AFTER_PRIMARY_UCSC_CHAIN":
                raise ValueError("Unverified source allele join")
            if x["valid_for_full_locus_coloc"]!="False":
                raise ValueError("Previously closed coloc gate changed")
            q=float(x["oasis_p"]);g=float(x["gwas_p"])
            if not (0<q<=1 and 0<g<=1):
                raise ValueError("Malformed source QTL/GWAS p")
            by[(x["gene"],x["cell"])].append(x)
    if set(by)!={("ALDH2","Mono-L1"),("ALDH2","B_Activated-L2"),
                 ("BRAP","Mono-L1"),("BRAP","B_Activated-L2"),
                 ("RPH3A","Mono-L1")}:
        raise ValueError("Unexpected gene-cell source state")
    return by

def analyze_group(gene,cell,items,genos,jpt):
    with_ref=[]
    for r in items:
        key=r["gwas_variant_grch37"]
        if key not in genos:
            continue
        eas=pairwise_ld(genos[ANCHOR],genos[key])
        jp=pairwise_ld(genos[ANCHOR],genos[key],jpt)
        if eas["r2"] is None or jp["r2"] is None:
            continue
        with_ref.append((r,eas,jp))
    if not with_ref:
        raise ValueError("No QC-passing overlapping LD positions")
    top_qtl=min(items,key=lambda r:(float(r["oasis_p"]),r["gwas_variant_grch37"]))
    top_gwas=min(items,key=lambda r:(float(r["gwas_p"]),r["gwas_variant_grch37"]))
    def lead_r2(row):
        key=row["gwas_variant_grch37"]
        if key not in genos:
            return {"in_reference":False, "EAS_r2":None,"JPT_r2":None}
        e=pairwise_ld(genos[ANCHOR],genos[key])
        j=pairwise_ld(genos[ANCHOR],genos[key],jpt)
        return {"in_reference":True, "EAS_r2":e["r2"],"JPT_r2":j["r2"]}
    for r in (top_qtl,top_gwas):
        if not float(r["oasis_p"])<=1 or not float(r["gwas_p"])<=1:
            raise ValueError("invalid top stats")
    high=[(r,e,j) for r,e,j in with_ref if e["r2"]>=.8]
    high_qtl_sig=[(r,e,j) for r,e,j in high if float(r["oasis_p"])<.05]
    all_qtl_sig=sum(float(r["oasis_p"])<.05 for r in items)
    return dict(
       gene=gene,cell=cell,n_exact_GWAS_OASIS_SNVS=len(items),
       n_genotypes_ref_ld_QC=len(with_ref),n_genotypes_ref_ld_absent_or_QC_failed=len(items)-len(with_ref),
       n_EAS_rs671_r2_ge_08=len(high),n_EAS_rs671_r2_ge_08_QTL_nominal_p_lt_05=len(high_qtl_sig),
       n_original_OASIS_QTL_nominal_p_lt_05=all_qtl_sig,
       strongest_oasis_QTL_SNP=top_qtl["gwas_variant_grch37"],
       strongest_oasis_QTL_rsid=top_qtl["oasis_rsid"],
       strongest_oasis_QTL_p=float(top_qtl["oasis_p"]),
       strongest_oasis_QTL_GWAS_p=float(top_qtl["gwas_p"]),
       strongest_oasis_QTL_ref_LD=lead_r2(top_qtl),
       strongest_GWAS_SNP_in_OASIS_intersect=top_gwas["gwas_variant_grch37"],
       strongest_GWAS_p_in_OASIS_intersect=float(top_gwas["gwas_p"]),
       strongest_GWAS_SNP_QTL_p=float(top_gwas["oasis_p"]),
       strongest_GWAS_SNP_ref_LD=lead_r2(top_gwas),
       median_EAS_r2_to_rs671=float(np.median([e["r2"] for r,e,j in with_ref])),
       median_JPT_r2_to_rs671=float(np.median([j["r2"] for r,e,j in with_ref])),
       highLD_QTL_min_p=min((float(r["oasis_p"]) for r,e,j in high),default=None),
       claim="PANEL_TAGGING_DESCRIPTIVE_ONLY")
def run(args):
    pairs=read_qtl_pairs(args.matched)
    requested={r["gwas_variant_grch37"] for v in pairs.values() for r in v}
    ld,jpt,pops=read_ld_genotypes(args.traw,args.panel,requested)
    anchor=ld[ANCHOR]
    # Recapitulate independent anchored LD values as a hard regression check.
    known=[]
    for marker,(rsid,expected) in EXPECTED_PIP_CS.items():
        a=pairwise_ld(anchor,ld[marker])
        b=pairwise_ld(anchor,ld[marker],jpt)
        if a["r2"] is None or abs(a["r2"]-expected)>.0001:
            raise ValueError("rs671->GWAS CS EAS reference LD no longer agrees with original")
        known.append(dict(rsid=rsid,variant_id=marker,
                          eas504_r2=a["r2"],jpt104_r2=b["r2"],
                          eas_pair_n=a["n_pair"],jpt_pair_n=b["n_pair"]))
    results=[analyze_group(g,c,items,ld,jpt) for (g,c),items in sorted(pairs.items())]
    all_ld_rows=[]
    for (gene,cell),items in sorted(pairs.items()):
        for item in items:
            vid=item["gwas_variant_grch37"]
            dos=ld.get(vid)
            if dos is None:
                continue
            eas=pairwise_ld(anchor,dos)
            jp=pairwise_ld(anchor,dos,jpt)
            if eas["r2"] is None or jp["r2"] is None:
                continue
            all_ld_rows.append({
                "gene":gene,"cell":cell,
                "GWAS_SNP_grch37":vid,
                "oasis_rsid":item["oasis_rsid"],
                "gwas_p":item["gwas_p"],
                "qtl_p":item["oasis_p"],
                "EAS504_signed_r_to_rs671":eas["r"],
                "EAS504_r2_to_rs671":eas["r2"],
                "JPT104_signed_r_to_rs671":jp["r"],
                "JPT104_r2_to_rs671":jp["r2"],
                "EAS_n_pair":eas["n_pair"],"JPT_n_pair":jp["n_pair"],
                "source_completeness":"BOKEH_ONLY",
                "coloc_claim":"NOT_VALID"
            })
    result={
       "audit":"IS_G0022_OASIS_RS671_EAS504_JPT104_REFERENCE_TAGGING_V1",
       "GWAS_and_QTL_exact_overlap_rows":sum(x["n_exact_GWAS_OASIS_SNVS"] for x in results),
       "gene_cell_layers":len(results),
       "LD_reference":"1000_Genomes_phase3_EAS_504; JPT_104_nested_subset_not_independent",
       "LD_genotypes_snp_count_loaded":len(ld),
       "JPT_reference_n":len(jpt),"reference_panel_pops":pops,
       "source_traw_sha256":hashlib.sha256(args.traw.read_bytes()).hexdigest(),
       "source_1kg_panel_sha256":hashlib.sha256(args.panel.read_bytes()).hexdigest(),
       "source_GWAS_OASIS_intersect_sha256":hashlib.sha256(args.matched.read_bytes()).hexdigest(),
       "rs671_anchor_ld_checkpoints":known,
       "gene_cell_results":results,
       "per_gene_cell_SNP_ref_tag_rows":len(all_ld_rows),
       "all_four_gwas_cs_covered_by_OASIS_browser":False,
       "published_OASIS_donor_genotype_LD_verified":False,
       "all_tested_OASIS_cis_SNP_denominator_verified":False,
       "effect_allele_direction_source_BIM_verified":False,
       "colocalizations_performed":0,"causal_gene_established":False,
       "scientific_warning":"Population-reference LD proximity of a marginal eQTL is NOT GWAS/QTL colocalization; EAS and JPT comparisons are nested, and nearby rs671 proxy is not independent evidence. Ranking genes by qtl lead p and LD is exploratory descriptive only."
    }
    args.outdir.mkdir(parents=True,exist_ok=True)
    with (args.outdir/"G0022_OASIS_EAS_JPT_RS671_REFERENCE_TAGGING_BY_SNP.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(all_ld_rows[0]),delimiter="\t")
        w.writeheader();w.writerows(all_ld_rows)
    p=args.outdir/"G0022_OASIS_EAS_JPT_RS671_REFERENCE_TAGGING_SUMMARY.json"
    p.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    table=[]
    for r in results:
        row={k:v for k,v in r.items() if not isinstance(v,dict)}
        row.update({"top_qtl_EAS_r2":r["strongest_oasis_QTL_ref_LD"]["EAS_r2"],
                    "top_qtl_JPT_r2":r["strongest_oasis_QTL_ref_LD"]["JPT_r2"],
                    "top_gwas_EAS_r2":r["strongest_GWAS_SNP_ref_LD"]["EAS_r2"],
                    "top_gwas_JPT_r2":r["strongest_GWAS_SNP_ref_LD"]["JPT_r2"]})
        table.append(row)
    with (args.outdir/"G0022_OASIS_GENE_CELL_RS671_LD_REFERENCE_SCORECARD.tsv").open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=list(table[0]),delimiter="\t")
        wr.writeheader();wr.writerows(table)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--traw",type=Path,required=True)
    p.add_argument("--panel",type=Path,required=True)
    p.add_argument("--matched",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    j=run(p.parse_args())
    for x in j["gene_cell_results"]:
        print(x["gene"],x["cell"],"overlap",x["n_exact_GWAS_OASIS_SNVS"],
           "EAS_ref",x["n_genotypes_ref_ld_QC"],"QTL_lead",x["strongest_oasis_QTL_rsid"],
           x["strongest_oasis_QTL_p"],"lead_EAS_r2",x["strongest_oasis_QTL_ref_LD"]["EAS_r2"],
           "lead_JPT_r2",x["strongest_oasis_QTL_ref_LD"]["JPT_r2"])
if __name__=="__main__":
    main()
