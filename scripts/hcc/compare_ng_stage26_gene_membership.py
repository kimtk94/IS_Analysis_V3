#!/usr/bin/env python3
"""Compare actual publisher Figure2b 5481 ENSG identities to 5343 current
strict reviewed UniProt 1:1 RNA/protein IDs; report only aggregate overlap.
"""
import argparse,hashlib,json,pathlib,re
ENSG=re.compile(r"^ENSG\d{11}$")
def read_ids(filepath, expected):
    p=pathlib.Path(filepath)
    if not p.is_file() or p.stat().st_size>300000:
        raise ValueError("runner-private source file missing or unusually large")
    labels=[line.strip() for line in p.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(labels)!=expected or len(set(labels))!=expected:
        raise ValueError("source original gene set expected exactly 5481 unique identifiers")
    if any(not ENSG.fullmatch(g) for g in labels):
        raise ValueError("unrecognized original Ensembl gene")
    return set(labels)
def read_current_pairs(filepath):
    p=pathlib.Path(filepath)
    if not p.is_file() or p.stat().st_size>300000:
        raise ValueError("runner-private current crosswalk not present")
    lines=[line for line in p.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines)!=5343:raise ValueError("current source one-to-one policy cardinality drift")
    pairs=[row.split("\t") for row in lines]
    if any(len(row)!=2 for row in pairs):raise ValueError("not tab-separated strict pairs")
    genes=[g for g,p in pairs]
    features=[p for g,p in pairs]
    if len(set(genes))!=5343 or len(set(features))!=5343:
        raise ValueError("not strictly one-to-one")
    if any(not ENSG.fullmatch(g) for g in genes):
        raise ValueError("current Ensembl gene not validated")
    return set(genes)
def compute(original,current):
    overlap=original & current
    result={
        "schema":"hcc-stage27-original-fig2b-5481-vs-current-5343-gene-membership-v3",
        "published_Ng_Fig2b_2022_ENSG_n":len(original),
        "current_Suter_strict_1to1_ENSG_n":len(current),
        "exact_Ensembl_gene_ID_overlap_n":len(overlap),
        "original_2022_genes_not_in_current_strict_n":len(original-current),
        "current_strict_genes_not_in_original_2022_n":len(current-original),
        "overlap_fraction_of_original_5481":round(len(overlap)/len(original),6),
        "overlap_fraction_of_current_5343":round(len(overlap)/len(current),6),
        "jaccard":round(len(overlap)/len(original|current),6),
        "original_gene_id_set_sha256":hashlib.sha256("\n".join(sorted(original)).encode()).hexdigest(),
        "current_gene_id_set_sha256":hashlib.sha256("\n".join(sorted(current)).encode()).hexdigest(),
        "shared_gene_id_set_sha256":hashlib.sha256("\n".join(sorted(overlap)).encode()).hexdigest(),
        "original_gene_set_semantics_reviewer_verified":True,
        "original_Figure2b_5481_gene_membership_verified":True,
        "source_sheet_label":"Fig2b",
        "official_Figure2b_source_relationship_verified":True,
        "exact_matching_of_original_2022_molecular_quantitation":False,
        "original_ng_2022_expression_RSEM_and_normal_reference_reproduced":False,
        "2022_Ensembl_vs_current_UniProt_mapping_versions_identical":False,
        "biological_significance_or_replication_tested":False,
        "patient_feature_or_original_gene_ID_lists_exported":False,
        "interpretation":"Exact ENSG string membership from official published Figure2b correlation Source Data vs 2026 reviewed-UniProt strict 1:1 subset. Original paper RNA/protein numerical calculations, annotation versions, and original sample cohort remain NOT reproduced.",
    }
    if result["exact_Ensembl_gene_ID_overlap_n"]==0:
        raise ValueError("unexpected empty overlap")
    return result
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--published-private",required=True)
    ap.add_argument("--current-private",required=True)
    ap.add_argument("--out",required=True)
    a=ap.parse_args()
    original=read_ids(a.published_private,5481)
    current=read_current_pairs(a.current_private)
    result=compute(original,current)
    p=pathlib.Path(a.out);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
if __name__=="__main__":main()
