#!/usr/bin/env python3
"""Audit six independent study families at G0022 4-SNP CS and eight genes.

Public eQTL Catalogue GRCh38 source records, reciprocal UCSC liftover and
explicit allele equivalence to GIGASTROKE EAS AIS GRCh37. No causal claims.
"""
import argparse,csv,gzip,json,hashlib,math
from collections import Counter,defaultdict
from pathlib import Path
from audit_g0022_gtex_v8_dapg import LiftOver,stable,complement
from prepare_g0022_eqtl_catalogue_aorta_coloc import parse_nominal,gwas as load_gwas_effects
from prepare_g0022_non_gtex_qtl_sources import OUT as SOURCE_LIST
GENES=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/IS_ALL_GENE_WINDOW_UNIVERSE.tsv")
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
RAW=Path("/srv/is-analysis/data/is/qtl/eqtl_catalogue/G0022_v1/non_gtex")
CHAIN=Path("/srv/is-analysis/data/is/reference/hg19ToHg38.over.chain.gz")
REVERSE=Path("/srv/is-analysis/data/reference/hg38ToHg19.over.chain")
OUTPUT=SOURCE_LIST
SYMS={"ACAD10","ALDH2","NAA25","HECTD4","BRAP","PTPN11","IFT81","ATP2A2"}
def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def freq(x):
    try:
        v=float(x)
        return v if math.isfinite(v) else None
    except (ValueError,TypeError):return None
def mapping(ref,lo,re):
    ret={}
    for r in ref:
        old=r["variant_id"];chrom,position,a,b=old.split(":")
        if chrom!="12":raise ValueError("Non-chr12 GWAS input")
        pos=int(position)
        forward=lo.convert_coordinate("chr12",pos-1)
        if len(forward)!=1:
            continue
        name,new0,strand,_=forward[0]
        if name!="chr12" or strand not in ("+","-"):
            continue
        back=re.convert_coordinate("chr12",int(new0))
        if len(back)!=1 or back[0][0]!="chr12" or int(back[0][1])!=pos-1:
            continue
        if strand=="-":a,b=complement(a),complement(b)
        key=(int(new0)+1,a,b)
        if key in ret and ret[key]["variant_id"]!=old:
            raise ValueError("Many-to-one SNP allele collision")
        ret[key]=r
    return ret
def build(root,raw,out,candidates,genes,chain,reverse):
    if LiftOver is None:raise RuntimeError("pyliftover unavailable")
    sources=rows(candidates)
    if len(sources)!=6:raise ValueError("Expected exactly six source cohorts")
    canonical=json.loads((root/"G0022_FULL_AIS_SUSIE_SUMMARY.json").read_text())
    cs={x for z in canonical["credible_set_variant_ids"].values() for x in z}
    if len(cs)!=4:raise ValueError("Expected four CS alleles")
    genotype=rows(root/"variants.tsv")
    gene=rows(genes)
    names={stable(x["gene_id"]):x["gene_symbol"] for x in gene
           if x["group_id"]=="IS_XDATA_G0022" and x["gene_symbol"] in SYMS}
    if set(names.values())!=SYMS:raise RuntimeError("Gene positional universe altered")
    l=LiftOver(str(chain));re=LiftOver(str(reverse))
    crosswalk=mapping(genotype,l,re)
    gwas_effects=load_gwas_effects(Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz"))
    by_id={x["variant_id"]:x for x in genotype}
    mapping_by_cs={}
    for key,record in crosswalk.items():
        if record["variant_id"] in cs:mapping_by_cs[record["variant_id"]]=key
    if set(mapping_by_cs)!=cs:raise ValueError("Unstable assembly conversion for GWAS CS")
    audit=[];by_gene=[]
    for s in sources:
        ds=s["dataset_id"];path=raw/f"{ds}_chr12_111650000_112370000.tsv.gz"
        if not path.exists():raise FileNotFoundError(path)
        count=Counter();positions=defaultdict(list)
        targets=defaultdict(dict)
        with gzip.open(path,"rt") as f:
            for line in f:
                if not line.strip():continue
                r=parse_nominal(line);count["source_records"]+=1
                if r["chromosome"]!="12":raise ValueError("Different chromosome in selected region")
                if r["type"]!="SNP" or len(r["ref"])!=1 or len(r["alt"])!=1:
                    count["non_SNP"]+=1;continue
                if r["variant"] != "chr12_{}_{}_{}".format(r["position"],r["ref"],r["alt"]):
                    raise RuntimeError("eQTL Catalogue SNP tuple mismatch")
                key=(int(r["position"]),r["ref"],r["alt"])
                for variant,cskey in mapping_by_cs.items():
                    if key[0]==cskey[0]:
                        positions[variant].append(r)
                symbol=names.get(stable(r["gene_id"]))
                if symbol is None:continue
                count["target_gene_source_records"]+=1
                gwas=crosswalk.get(key)
                if gwas is None:
                    count["no_exact_GWAS_match"]+=1
                    continue
                qc=(r["beta"],r["se"],r["pvalue"],r["maf"])
                if any(freq(val) is None for val in qc) or freq(r["se"])<=0:
                    count["invalid_effect"]+=1;continue
                if {gwas["ref"],gwas["alt"]} in ({"A","T"},{"G","C"}):
                    count["palindrome"]+=1;continue
                effect=gwas_effects.get(gwas["variant_id"])
                if effect is None:raise RuntimeError("Full-locus SNP missing source GWAS effect")
                old=targets[symbol].get(gwas["variant_id"])
                if old:
                    if old["qtl_beta_alt"]!=r["beta"] or old["qtl_se"]!=r["se"]:
                        raise ValueError("Duplicate conflicting SNP-gene effect")
                    count["duplicate_rsid_record"]+=1;continue
                targets[symbol][gwas["variant_id"]]={
                  "study_id":s["study_id"],"dataset_id":ds,"study_label":s["study_label"],
                  "tissue":s["sample_group"],"gene_symbol":symbol,
                  "gene_id":stable(r["gene_id"]),"gwas_snp_grch37":gwas["variant_id"],
                  "qtl_snp_grch38":r["variant"],"gwas_beta_alt":effect["alt_effect_beta"],
                  "gwas_se":effect["se"],"gwas_p":effect["p"],
                  "qtl_beta_alt":r["beta"],"qtl_se":r["se"],"qtl_p":r["pvalue"],
                  "qtl_maf":r["maf"],"in_gwas_CredibleSet":int(gwas["variant_id"] in cs),
                  "scientific_gate":"NOT_VALIDATED_AIS_QTL_COLOCALIZATION"}
                count["gene_snp_matched"]+=1
        coverage={}
        for variant,key in mapping_by_cs.items():
            hit=positions.get(variant,[])
            exact=[r for r in hit if (int(r["position"]),r["ref"],r["alt"])==key]
            swapped=[r for r in hit if int(r["position"])==key[0] and
                     (r["ref"],r["alt"])==(key[2],key[1])]
            status=("POSITION_ABSENT" if not hit else
                  "EXACT_REF_ALT_PRESENT" if exact else
                  "ALLELE_PAIR_SWAPPED_PRESENT" if swapped else
                  "POSITION_PRESENT_ALLELE_DIFFERENT")
            audit.append({
                "dataset_id":ds,"study_label":s["study_label"],"tissue":s["sample_group"],
                "sample_size_manifest":s["sample_size"],
                "GWAS_CS_GRCh37":variant,"QTL_grch38_position":key[0],
                "QTL_expected_ref":key[1],"QTL_expected_alt":key[2],
                "qtl_records_at_position":len(hit),
                "qtl_records_matching_ref_alt":len(exact),
                "allele_QC_status":status,
                "interpretation":"MEASUREMENT_COVERAGE_ONLY_NOT_BIOLOGICAL_NEGATIVE"})
            coverage[variant]=status
        print(json.dumps({"dataset":ds,"study":s["study_label"],"records":count["source_records"],
            "matched_gene_snps":count["gene_snp_matched"],"CS":coverage}),flush=True)
        for sym in sorted(SYMS):
            vals=list(targets.get(sym,{}).values())
            vals.sort(key=lambda x:(int(x["gwas_snp_grch37"].split(":")[1]),x["gwas_snp_grch37"]))
            if vals:
                path_out=out/ds
                path_out.mkdir(parents=True,exist_ok=True)
                with (path_out/f"G0022_{ds}_{sym}_AIS_QTL_INPUT.tsv").open("w",newline="") as f:
                    w=csv.DictWriter(f,fieldnames=list(vals[0]),delimiter="\t");w.writeheader();w.writerows(vals)
            by_gene.append({
              "dataset_id":ds,"study_label":s["study_label"],"tissue":s["sample_group"],
              "gene":sym,"n_matched_allele_snp":len(vals),
              "min_gwas_p":min((float(x["gwas_p"]) for x in vals),default=""),
              "min_qtl_p":min((float(x["qtl_p"]) for x in vals),default=""),
              "CS_variant_matched":sum(x["in_gwas_CredibleSet"] for x in vals),
              "qtl_status":"EXPLORATORY_NOT_VALIDATED"})
    out.mkdir(parents=True,exist_ok=True)
    for filename,items in (
      ("G0022_NON_GTEX_CS_SOURCE_COVERAGE.tsv",audit),
      ("G0022_NON_GTEX_GENE_SNP_COVERAGE.tsv",by_gene)):
        with (out/filename).open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(items[0]),delimiter="\t");w.writeheader();w.writerows(items)
    by_status=Counter(a["allele_QC_status"] for a in audit)
    summary={
        "source_dataset_count":len(sources),"cs_variant_dataset_tests":len(audit),
        "study_ancestry_verified_as_EAS":0,
        "source_test_status":dict(by_status),
        "gene_snp_matched_total":sum(x["n_matched_allele_snp"] for x in by_gene),
        "gene_tissue_combinations":len(by_gene),
        "dataset_counts":Counter(x["dataset_id"] for x in by_gene),
        "causal_gene_validated":0,
        "scientific_gate":"SOURCE_VARIANT_COVERAGE_ONLY_DO_NOT_CALL_POSITIVE_COLOCALIZATION"}
    (out/"G0022_NON_GTEX_QTL_CS_COVERAGE_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    return summary
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--raw",type=Path,default=RAW)
    p.add_argument("--out",type=Path,default=OUTPUT)
    p.add_argument("--candidates",type=Path,default=SOURCE_LIST/"G0022_NON_GTEX_QTL_SOURCE_CANDIDATES.tsv")
    p.add_argument("--genes",type=Path,default=GENES)
    p.add_argument("--chain",type=Path,default=CHAIN)
    p.add_argument("--reverse",type=Path,default=REVERSE)
    a=p.parse_args()
    build(a.root,a.raw,a.out,a.candidates,a.genes,a.chain,a.reverse)
