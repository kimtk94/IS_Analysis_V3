#!/usr/bin/env python3
"""G0022 GTEx gene-level eQTL Catalogue: multi-tissue assembly+allele audit.

Unmodified nominal QTL GRCh38 regional extracts -> exact GRCh37 REF/ALT matches
with EAS AIS GWAS. Coverage-first discovery: NEVER call "no QTL" from missing SNPs.
"""
import argparse,csv,gzip,hashlib,json,math,sys
from collections import defaultdict,Counter
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_g0022_gtex_v8_dapg import variant_liftover,stable,LiftOver
from prepare_g0022_eqtl_catalogue_aorta_coloc import (
    parse_nominal,gwas,candidate_genes,GENES,CHAIN,GWAS,ANNOTATION)
RAW=Path("/srv/is-analysis/data/is/qtl/eqtl_catalogue/G0022_v1")
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
OUT=ROOT/"eqtl_catalogue_multitissue_v1"
REGION="12:111650000-112370000"
SOURCES={
    "QTD000131":("artery_aorta",387),
    "QTD000136":("artery_coronary",213),
    "QTD000141":("artery_tibial",584),
    "QTD000171":("brain_cortex",205),
}
def source_file(root,study):
    tissue,n=SOURCES[study]
    return root/f"{study}_{tissue}_GRCh38_chr12_111650000_112370000.tsv.gz"
def num(x):
    try:
        v=float(x)
        return v if math.isfinite(v) else None
    except (TypeError,ValueError):return None
def write_tsv(path,data,fields):
    with path.open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",extrasaction="ignore")
        w.writeheader();w.writerows(data)
def process_study(study,source,genes,gwas_data,cs_snps,lifter,out,min_source_bytes=1000):
    tissue,sample_size=SOURCES[study]
    if not source.is_file():raise FileNotFoundError(source)
    if source.stat().st_size<min_source_bytes:raise RuntimeError("Missing/too-small real QTL extract "+str(source))
    match=defaultdict(dict)
    qc=Counter()
    all_cs_present=set()
    gwas_cs_matched=set()
    # Count sources by variant/genes after allele-aware liftOver.
    with gzip.open(source,"rt") as h:
        for line in h:
            if not line.strip():continue
            row=parse_nominal(line)
            qc["source_records"]+=1
            if row["chromosome"]!="12":raise RuntimeError("Unexpected chromosome")
            if row["type"]!="SNP" or len(row["ref"])!=1 or len(row["alt"])!=1:
                qc["non_snp"]+=1;continue
            if row["variant"]!=f'chr12_{row["position"]}_{row["ref"]}_{row["alt"]}':
                raise RuntimeError("QTL variant genomic tuple disagreement")
            locus,state=variant_liftover(row["variant"]+"_b38",lifter)
            qc[state]+=1
            if locus is None:continue
            k,strand=locus
            if k in cs_snps:all_cs_present.add(k)
            gene_id=stable(row["gene_id"])
            if gene_id not in genes:continue
            qc["target_gene_source_records"]+=1
            if k not in gwas_data:
                qc["not_in_EAS_AIS_harmonized_GWAS"]+=1;continue
            if {gwas_data[k]["ref"],gwas_data[k]["alt"]} in ({"A","T"},{"G","C"}):
                qc["palindromic_excluded"]+=1;continue
            beta=num(row["beta"]);se=num(row["se"]);p=num(row["pvalue"])
            maf=num(row["maf"]);ac=num(row["ac"]);an=num(row["an"])
            if None in (beta,se,p,maf) or se<=0 or not (0<=p<=1 and 0<maf<=.5):
                qc["invalid_qtl_effect"]+=1;continue
            if ac is not None and an is not None and (an<=0 or not 0<=ac<=an):
                qc["invalid_allele_count"]+=1;continue
            if ac is not None and an is not None:
                af=ac/an
                if abs(min(af,1-af)-maf)>.02:qc["minor_AF_disagrees_with_AC_AN"]+=1
            else:af=None
            g=gwas_data[k]
            record={
              "source_dataset_id":study,"source_tissue":tissue,
              "gene_id":gene_id,"gene_symbol":genes[gene_id],
              "variant_id_grch37":k,"source_variant_id_grch38":row["variant"],
              "liftover_strand":strand,
              "gwas_beta_alt":g["alt_effect_beta"],"gwas_se":g["se"],
              "gwas_p":g["p"],"gwas_alt_eaf":g["alt_effect_eaf"],
              "qtl_beta_alt":row["beta"],"qtl_se":row["se"],"qtl_p":row["pvalue"],
              "qtl_maf":row["maf"],"qtl_ac":row["ac"],"qtl_an":row["an"],
              "qtl_alt_eaf":round(af,10) if af is not None else "",
              "qtl_study_n_catalogue":sample_size,
              "gwas_95pct_cs_snp":int(k in cs_snps),
              "qtl_effect_allele":"ALT","gwas_effect_allele":"ALT",
              "scientific_gate":"MATCHED_SNP_ONLY_NOT_VALIDATED_COLOCALIZATION"}
            prev=match[gene_id].get(k)
            if prev:
                for fld in ("qtl_beta_alt","qtl_se","qtl_p","qtl_maf"):
                    if prev[fld]!=record[fld]:
                        raise RuntimeError("Same gene+SNP has conflicting QTL effect")
                qc["duplicate_same_gene_snp"]+=1
                continue
            match[gene_id][k]=record
            qc["matched_gene_snp"]+=1
            if k in cs_snps:gwas_cs_matched.add(k)
    out.mkdir(parents=True,exist_ok=True)
    audit=[]
    for gene_id,symbol in sorted(genes.items(),key=lambda a:a[1]):
        r=list(match.get(gene_id,{}).values())
        r.sort(key=lambda x:(int(x["variant_id_grch37"].split(":")[1]),x["variant_id_grch37"]))
        if r:
            write_tsv(out/f"G0022_{symbol}_AIS_{study}_ABF_INPUT.tsv",r,list(r[0]))
        ngwas=min((float(x["gwas_p"]) for x in r),default=None)
        nqtl=min((float(x["qtl_p"]) for x in r),default=None)
        lcs={x["variant_id_grch37"] for x in r if x["gwas_95pct_cs_snp"]}
        has_gws=ngwas is not None and ngwas<=5e-8
        status=("READY_FOR_PILOT_COLOC_DIAGNOSTIC_NOT_VALIDATION" if
                len(r)>=100 and has_gws and nqtl is not None and nqtl<=1e-4 and len(lcs)>0 else
                "MISSING_GWS_OR_CS_AND_OR_QTL_SIGNAL_NOT_NEGATIVE")
        audit.append({
          "dataset_id":study,"tissue":tissue,"gene":symbol,"gene_id":gene_id,
          "matched_gene_snp":len(r),
          "min_EAS_AIS_p":ngwas if ngwas is not None else "",
          "min_qtl_p":nqtl if nqtl is not None else "",
          "cs_matched_distinct_variants":len(lcs),
          "matched_cs_ids":";".join(sorted(lcs)),
          "source_cs_snp_any_gene_total":len(all_cs_present),
          "status":status})
    write_tsv(out/"G0022_MULTITISSUE_QTL_INPUT_AUDIT.tsv",audit,list(audit[0]))
    report={
      "status":"ALLELE_HARMONIZATION_COMPLETE_COLOC_UNVALIDATED",
      "dataset_id":study,"tissue":tissue,"catalogue_sample_size":sample_size,
      "region_grch38":REGION,"raw_path":str(source),"raw_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
      "genes_checked":len(genes),"target_gene_snp_matches":sum(x["matched_gene_snp"] for x in audit),
      "genes_with_matched_snps":sum(x["matched_gene_snp"]>0 for x in audit),
      "genes_ready_for_diagnostic":sum(x["status"].startswith("READY") for x in audit),
      "all_four_GWAS_CS_SNPs":sorted(cs_snps),
      "GWAS_CS_SNPs_in_QTL_extract_any_gene":sorted(all_cs_present),
      "GWAS_CS_SNPs_matched_in_8_target_genes":sorted(gwas_cs_matched),
      "qc":dict(qc),"warning":"No validated causal genes or molecular colocalization; source only regional nominal eQTL"}
    (out/"G0022_MULTITISSUE_QTL_INPUT_SUMMARY.json").write_text(json.dumps(report,indent=2))
    print(json.dumps({k:report[k] for k in (
      "dataset_id","tissue","target_gene_snp_matches","genes_with_matched_snps",
      "genes_ready_for_diagnostic","GWAS_CS_SNPs_in_QTL_extract_any_gene",
      "GWAS_CS_SNPs_matched_in_8_target_genes")},indent=2),flush=True)
    return report

def run(args):
    ref=json.loads((args.full/"G0022_FULL_AIS_SUSIE_SUMMARY.json").read_text())
    if ref["status"]!="FULL_LOCUS_EXPLORATORY_ONLY":raise ValueError("full-locus source QC missing")
    cs={v for arr in ref["credible_set_variant_ids"].values() for v in arr}
    if len(cs)!=4:raise ValueError("Expected four full-locus CS variants")
    genes=candidate_genes(args.annotation);gw=gwas(args.gwas)
    if not cs.issubset(gw):raise ValueError("Full CS SNPs missing from GWAS")
    if LiftOver is None:raise RuntimeError("pyliftover isolated install required")
    converter=LiftOver(str(args.chain))
    datasets=args.datasets.split(",")
    if len(datasets)!=len(set(datasets)) or any(x not in SOURCES for x in datasets):
        raise ValueError("Unknown/repeated QTL Catalogue dataset")
    summaries=[]
    for study in datasets:
        path=source_file(args.raw,study)
        summaries.append(process_study(study,path,genes,gw,cs,converter,args.out/study))
    output={"datasets":len(summaries),"dataset_ids":datasets,
      "total_gene_snp_matches":sum(x["target_gene_snp_matches"] for x in summaries),
      "datasets_with_any_GWAS_CS_variant":sum(bool(x["GWAS_CS_SNPs_matched_in_8_target_genes"]) for x in summaries),
      "scientific_result":"HARMONIZATION_ONLY_NOT_VALIDATED_COLOCALIZATION"}
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"G0022_MULTITISSUE_QTL_HARMONIZATION_MASTER.json").write_text(json.dumps(output,indent=2))
    print(json.dumps(output,indent=2),flush=True)
    return output
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--datasets",default="QTD000131,QTD000136,QTD000141,QTD000171")
    ap.add_argument("--raw",type=Path,default=RAW)
    ap.add_argument("--out",type=Path,default=OUT)
    ap.add_argument("--chain",type=Path,default=CHAIN)
    ap.add_argument("--gwas",type=Path,default=GWAS)
    ap.add_argument("--annotation",type=Path,default=ANNOTATION)
    ap.add_argument("--full",type=Path,default=ROOT)
    run(ap.parse_args())
