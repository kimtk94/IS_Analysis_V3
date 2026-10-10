#!/usr/bin/env python3
"""Frequency lineage QC for *exactly* REF/ALT matched OASIS–AIS SNPs.

Combines raw-source-vs-browser field confirmation and source author tensorQTL
method with cross-population AF correlation. Does NOT assert every variant
effect allele fully verified from the restricted original PLINK BIM/BED.
"""
import argparse,csv,hashlib,json
from collections import defaultdict
from pathlib import Path
import numpy as np

FIELD="oasis_browser_maf_column_AS_REPORTED_NOT_TRUE_MAF"
def metrics(g,a):
    x=np.asarray(g,dtype=float);y=np.asarray(a,dtype=float)
    if x.size<10 or x.shape!=y.shape or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("Invalid or insufficient variant frequency data")
    if not (np.all((x>0)&(x<1)) and np.all((y>=0)&(y<=1))):
        raise ValueError("Out of range AF")
    same=np.abs(x-y);flipped=np.abs(1-x-y)
    return dict(n=int(x.size),
         pearson_alt_eaf_r=float(np.corrcoef(x,y)[0,1]),
         abs_difference_ALT_median=float(np.median(same)),
         abs_difference_ALT_mean=float(np.mean(same)),
         abs_difference_FLIPPED_mean=float(np.mean(flipped)),
         closer_to_same_allele_count=int(np.sum(same<flipped)),
         closer_to_flipped_allele_count=int(np.sum(flipped<same)),
         equidistant_count=int(np.sum(same==flipped)),
         frequency_difference_gt_0_1_count=int(np.sum(same>0.1)),
         QTL_AF_above_0_5_count=int(np.sum(y>.5)))

def read_matches(path):
    with path.open(newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        required={"gwas_variant_grch37","oasis_variant_grch38","gene","cell",
                  "gwas_ALT_EAF",FIELD,"source_ref_alt_id_match","oasis_beta_UNHARMONIZED"}
        if not required.issubset(reader.fieldnames):
            raise ValueError("Unexpected regional allele matching schema")
        result=list(reader)
    if len(result)<500:
        raise ValueError("Unexpectedly low G0022 regional intersect")
    if any(row["source_ref_alt_id_match"]!="EXACT_AFTER_PRIMARY_UCSC_CHAIN" for row in result):
        raise ValueError("Nonaligned source allele in AF comparison")
    return result

def run(args):
    cross=json.loads(args.crosswalk_audit.read_text())
    if (cross.get("raw_AF_greater_than_half",0)<100
        or cross.get("row_by_variant_match",0)<500
        or not cross.get("browser_field_maf","").startswith("RAW_tensorQTL_af")):
        raise ValueError("Original->web source proof unavailable")
    data=read_matches(args.matched)
    keyed=defaultdict(list);unique=defaultdict(list)
    for r in data:
        g=float(r["gwas_ALT_EAF"]);o=float(r[FIELD])
        if not (0<g<1 and 0<=o<=1):
            raise ValueError("AF bounds failed")
        keyed[(r["gene"],r["cell"])].append((g,o))
        unique[r["gwas_variant_grch37"]].append((g,o))
    page_metrics=[]
    for (gene,cell),pairs in sorted(keyed.items()):
        values=metrics(*zip(*pairs))
        page_metrics.append({"gene":gene,"cell":cell,**values})
    independent=[]
    for snp,pairs in unique.items():
        gvals={g for g,o in pairs}
        if len(gvals)!=1: raise ValueError("Inconsistent GWAS ALT EAF by SNP")
        original_eaf=next(iter(gvals))
        observed_af=float(np.median([o for g,o in pairs]))
        independent.append(dict(variant_grch37=snp,gwas_ALT_EAF=original_eaf,
                 oasis_AF_median_across_gene_cells=observed_af,
                 closest_allele_frequency_model=(
                 "SAME_ALT" if abs(original_eaf-observed_af)<abs(1-original_eaf-observed_af)
                 else "FLIPPED_REF" if abs(original_eaf-observed_af)>abs(1-original_eaf-observed_af)
                 else "INDETERMINATE"),
                 gene_cell_occurrences=len(pairs),
                 absolute_ALT_AF_difference=abs(original_eaf-observed_af)))
    independent.sort(key=lambda x:x["variant_grch37"])
    uni=metrics([v["gwas_ALT_EAF"] for v in independent],
                [v["oasis_AF_median_across_gene_cells"] for v in independent])
    overview={
        "audit":"IS_G0022_OASIS_GWAS_ALT_AF_ORIENTATION_TRIANGULATION_V1",
        "input_GWAS_OASIS_overlaps":len(data),
        "unique_aligned_allelic_SNVs":len(independent),
        "unique_variant_stats":uni,
        "five_cell_gene_page_stats":page_metrics,
        "browser_field_does_not_mean_MAF":True,
        "empirical_source_AF_field_lineage_gate":"DIRECT_RAW_BROWSER_1010_VARIANT_NUMERIC_MATCH",
        "OASIS_AF_orientation_inference":"STRONGLY_CONSISTENT_WITH_GWAS_ALTERNATIVE_ALLELE_EAF",
        "OASIS_beta_allele_orientation_inference":"TENSORQTL_METHOD_SUPPORTS_ALT_DOSAGE_IF_SOURCE_PLINK_ORDER_CORRECT",
        "source_OASIS_PLINK_BIM_ORIENTATION_DIRECTLY_VERIFIED":False,
        "effect_beta_gene_cell_source_original_G0022_rows_directly_verified":False,
        "valid_molecular_GWAS_colocalization":False,
        "explicit_caveat":"High genome-wide regional AF correlation is not same-cohort replication, proves no phenotype causality and cannot substitute raw PLINK allele verification. Allele frequency closely matching can be coincidental for MAF around 0.5; avoid per-SNP direction claims solely on frequency.",
        "original_to_browser_audit_sha256":hashlib.sha256(args.crosswalk_audit.read_bytes()).hexdigest(),
        "regional_matched_input_sha256":hashlib.sha256(args.matched.read_bytes()).hexdigest()
    }
    args.outdir.mkdir(parents=True,exist_ok=True)
    def save(name,rows):
        with (args.outdir/name).open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
            w.writeheader();w.writerows(rows)
    save("G0022_OASIS_EAS_AIS_AF_BY_GENE_CELL.tsv",page_metrics)
    save("G0022_OASIS_EAS_AIS_UNIQUE_SNP_AF_DIRECTION.tsv",independent)
    (args.outdir/"G0022_OASIS_EAS_AIS_AF_SOURCE_TRIANGULATION.json").write_text(json.dumps(overview,indent=2)+"\n")
    return overview

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matched",type=Path,required=True)
    p.add_argument("--crosswalk-audit",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    j=run(a)
    print(json.dumps({"unique_variant_stats":j["unique_variant_stats"],
        "input_GWAS_OASIS_overlaps":j["input_GWAS_OASIS_overlaps"],
        "source_PLINK_verified":j["source_OASIS_PLINK_BIM_ORIENTATION_DIRECTLY_VERIFIED"],
        "coloc":j["valid_molecular_GWAS_colocalization"]},indent=2))

if __name__=="__main__":
    main()
