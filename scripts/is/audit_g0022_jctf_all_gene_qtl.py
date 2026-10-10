#!/usr/bin/env python3
"""Preserve ALL genes shown for G0022 four-JCTF-variant QTL; no narrow shortlist."""
import argparse,csv,json
from collections import defaultdict
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
SOURCE=ROOT/"jctf_japan_omics_variants_v1/G0022_JCTF_4SNP_JAPANESE_EQTL_PQTL_EVIDENCE.tsv"
def build(source,root):
    with source.open() as f:rows=list(csv.DictReader(f,delimiter="\t"))
    if len(rows)!=72:raise ValueError("Source must be 72 molecular associations across four CS variants")
    groups=defaultdict(list)
    for x in rows:
        k=(x["gene_id"],x["gene_symbol"],x["qtl_category"])
        groups[k].append(x)
    genes=set(x["gene_symbol"] for x in rows)
    data=[]
    for (gid,sym,qtl),group in sorted(groups.items()):
        p=[float(x["jctf_allele_association_p"]) for x in group]
        susie=[float(x["jctf_susie_variant_pip"]) for x in group]
        fm=[float(x["jctf_finemap_variant_pip"]) for x in group]
        top=max(range(len(susie)),key=lambda i:susie[i])
        data.append({"gene_id":gid,"gene_symbol":sym,"molecular_trait":qtl,
         "number_of_GWAS_CS_variants_in_JCTF":len(group),
         "minimum_JCTF_association_p_among_four":min(p),
         "maximum_JCTF_SuSiE_variant_PIP_among_four":susie[top],
         "maximum_JCTF_FINEMAP_variant_PIP_among_four":max(fm),
         "variant_with_max_JCTF_SuSiE_PIP":group[top]["GWAS_variant_GRCh37"],
         "assay_type":"BLOOD_RNA_SEQ" if qtl=="eQTL" else "OLINK_EXPLORE",
         "full_JCTF_QTL_variants_analyzed":False,
         "full_locus_AIS_QTL_coloc":False,
         "protein_assay_epitope_QC_if_pQTL":"NOT_RESOLVED" if qtl=="pQTL" else "NOT_APPLICABLE",
         "report_category":"MOLECULAR_ASSOCIATION_NOT_MEDIATION"})
    root.mkdir(parents=True,exist_ok=True)
    with (root/"G0022_JCTF_ALL_GENE_QTL_SIGNAL_AUDIT.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter="\t");w.writeheader();w.writerows(data)
    result={"source_rows":len(rows),"unique_genes":len(genes),
      "gene_trait_pairs":len(data),
      "genes_with_any_eQTL":len(set(x["gene_symbol"] for x in rows if x["qtl_category"]=="eQTL")),
      "genes_with_any_pQTL":len(set(x["gene_symbol"] for x in rows if x["qtl_category"]=="pQTL")),
      "all_4_snp_same_causal_variant_assumption":False,
      "validated_AIS_gene_mediation":0,
      "preserve_gene_diversity":True}
    (root/"G0022_JCTF_ALL_GENE_QTL_SUMMARY.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    for row in sorted(data,key=lambda x:-x["maximum_JCTF_SuSiE_variant_PIP_among_four"])[:12]:
        print(row["gene_symbol"],row["molecular_trait"],row["minimum_JCTF_association_p_among_four"],
           row["maximum_JCTF_SuSiE_variant_PIP_among_four"])
    return result
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path,default=SOURCE)
    ap.add_argument("--root",type=Path,default=ROOT)
    a=ap.parse_args();build(a.source,a.root)
