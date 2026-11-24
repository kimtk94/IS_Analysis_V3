#!/usr/bin/env python3
"""Fail-closed verifier for *selected reduced* Ng/Suter HCCinfo gene mapping."""
import argparse,json,pathlib,re

SRC="438f52271d0a86a6701baa7a96e80c516ba575c2"
def check(rec):
    if rec.get("status")!="SOURCE_REDUCED_GENE_CROSSWALK_PRESENT":
        raise ValueError("source-selected gene crosswalk not verified")
    if rec.get("source_commit")!=SRC or rec.get("source_mapping_object_class") not in ("bnInfo","bnInfo_like_allnames"):
        raise ValueError("unpinned or unrecognized mapping source")
    for k in ("HCCdata_sha256","HCCinfo_sha256"):
        if not re.fullmatch(r"[0-9a-f]{64}",str(rec.get(k,""))):
            raise ValueError("missing source SHA256")
    if rec.get("samples_per_modality")!=50 or rec.get("transcript_feature_total")!=188 or rec.get("protein_feature_total")!=116:
        raise ValueError("wrong reduced matrix geometry")
    for n,total in (("T",188),("P",116)):
        m=rec.get(n+"_source_feature_keys_matched_n")
        u=rec.get(n+"_unique_gene_symbols_mapped_n")
        if not isinstance(m,int) or not isinstance(u,int) or not (0<=u<=m<=total):
            raise ValueError("invalid source gene label mapping count")
    overlap=rec.get("T_P_shared_gene_symbols_n")
    unique=rec.get("T_P_unambiguous_one_to_one_genes_n")
    if not all(isinstance(v,int) for v in (overlap,unique)) or not(0<unique<=overlap<=116):
        raise ValueError("no unambiguous gene-level crossomics anchors")
    for k in ("full_feature_Ensembl_protein_mapping_verified",
              "source_gene_mapping_independently_curated",
              "clinical_effects_verified",
              "patient_ids_feature_names_expression_values_exported"):
        if rec.get(k) is not False:raise ValueError("unsupported scientific or privacy claim")
    return {"status":"PASS_SOURCE_REDUCED_GENE_CROSSWALK_GEOMETRY",
            "shared_genes_in_selected_188T_116P":overlap,
            "unambiguous_one_to_one_selected_gene_anchors":unique,
            "full_quantitative_gene_wide_mapping":"HOLD",
            "clinical_or_causal_significance":"NOT_EVALUATED",
            "identifiers_or_per_sample_values_released":False}
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--source",required=True);args=ap.parse_args()
    p=pathlib.Path(args.source)
    if not p.is_file() or p.stat().st_size>100000:raise SystemExit("source missing/oversized")
    print(json.dumps(check(json.loads(p.read_text(encoding="utf-8"))),indent=2))
if __name__=="__main__":main()
