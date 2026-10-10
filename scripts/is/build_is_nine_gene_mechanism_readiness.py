#!/usr/bin/env python3
"""Descriptive nine-gene coloc-prior + six-donor healthy brain reference evidence.

Do not interpret absence from feature matrix as no expression; do not interpret
healthy localization as stroke disease-cell-specific expression or mediator.
The chosen nine genes are historic exploratory comparisons, NOT nine validated
causal genes and NOT a filter on the other 2,216 positional candidates.
"""
import argparse,csv,hashlib,json,math
from pathlib import Path
from collections import Counter,defaultdict
NINE=("FGF5","ALDH2","SH3PXD2A","COL4A2","COL4A1",
      "CALHM2","NEURL1","C4orf22","INA")
def read(path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f,delimiter="\t"))
def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(4194304),b""):h.update(b)
    return h.hexdigest()
def numeric(s):
    try:
        a=float(s)
        return a if math.isfinite(a) else None
    except (TypeError,ValueError):
        return None
def audited(genes,cell,features,donor,grid,human_manifest):
    if len(genes)!=2225 or len({r["gene_id_stable"] for r in genes})!=2225:
        raise ValueError("2,225 positional Ensembl genes required")
    if len(cell)!=9 or len(features)!=9 or len(grid)!=3230:
        raise ValueError("Nine-gene donor / 646 x five-prior model inputs required")
    if human_manifest["donor_analysis_unit"]!="Patient" or human_manifest["reference_type"]!="ADULT_CONTROL_TEMPORAL_LOBE_VASCULAR_PERIVASCULAR":
        raise ValueError("Human reference differs from frozen adult healthy control")
    bycell={r["gene"]:r for r in cell}
    byfeat={r["gene"]:r for r in features}
    if set(NINE)!=set(bycell) or set(NINE)!=set(byfeat):
        raise ValueError("Expected exact historic gene reference panel")
    geneidx={r["gene_id_stable"]:r for r in genes}
    symbol_map={}
    for name in NINE:
        alias=("NEURL" if name=="NEURL1" else name)
        hits=[g for g in genes if g["canonical_GENCODE_v19_gene_symbol"]==alias]
        if len(hits)!=1:raise ValueError("Ambiguous canonical stable ID "+name)
        if name=="NEURL1" and hits[0]["gene_id_stable"]!="ENSG00000107954":
            raise ValueError("Legacy NEURL1 / NEURL mapping drift")
        symbol_map[name]=hits[0]
    # Archive GTEx ABF for the nine candidates only, exact stable ID, never
    # intermix healthy reference gene symbols with unrelated QTL gene IDs.
    bygid=defaultdict(list)
    for r in grid:
        bygid[r["gene_base"]].append(r)
    pair=defaultdict(list)
    for r in donor:
        pair[r["gene"]].append(r)
    out=[]
    for name in NINE:
        g=symbol_map[name]
        feature=byfeat[name]["feature_status"]
        present=bycell[name]["present"]=="TRUE"
        if present!=(feature=="FEATURE_PRESENT"):
            raise ValueError("Healthy donor presence disagrees with feature QC for "+name)
        selected=bygid.get(g["gene_id_stable"],[])
        baseline=[r for r in selected if abs(float(r["p12"])-1e-5)<1e-13]
        high=[r for r in selected if abs(float(r["p12"])-1e-4)<1e-13]
        strict=[r for r in selected if abs(float(r["p12"])-1e-6)<1e-13]
        if len(baseline)!=len(high) or len(baseline)!=len(strict):
            raise ValueError("Coloc prior grid incomplete for "+name)
        if len(baseline)!=int(g["new_DIRECT_prior_grid_assays"]):
            raise ValueError("Archived QTL assay count mismatch for "+name)
        best=max(baseline,key=lambda r:float(r["PP_H4"]),default=None)
        best_h3=float(best["PP_H3"]) if best else None
        best_h4=float(best["PP_H4"]) if best else None
        bestkey=(best["locus"],best["dataset_key"],best["gene_base"]) if best else None
        matchstrict=next((r for r in strict if (r["locus"],r["dataset_key"],r["gene_base"])==bestkey),None)
        matchhigh=next((r for r in high if (r["locus"],r["dataset_key"],r["gene_base"])==bestkey),None)
        pair_list=pair.get(name,[])
        passing=[r for r in pair_list if r["qc_status"]=="DESCRIPTIVE_PAIRED_QC_PASS"]
        if len(pair_list)!=len(passing):raise ValueError("Healthy donor QC audit incomplete for "+name)
        out.append({
          "historic_symbol":name,
          "canonical_ensembl_stable_id":g["gene_id_stable"],
          "current_GENCODE_v19_symbol":g["canonical_GENCODE_v19_gene_symbol"],
          "n_positional_regions":g["n_positional_regions"],
          "n_BBJ_GTEx_v8_archived_QTL_assays":len(baseline),
          "best_baseline_tissue":best["dataset_key"] if best else "NO_ARCHIVED_ASSAY",
          "best_baseline_BBJ_locus":best["locus"] if best else "",
          "best_H3_original_p12_1e5":best_h3 if best else "",
          "best_H4_original_p12_1e5":best_h4 if best else "",
          "same_assay_H4_strict_p12_1e6":float(matchstrict["PP_H4"]) if matchstrict else "",
          "same_assay_H4_high_p12_1e4":float(matchhigh["PP_H4"]) if matchhigh else "",
          "best_baseline_H4_gt_H3":best_h4>best_h3 if best else "NOT_TESTED",
          "H4_ge_0p8_original":best_h4>=.8 if best else "NOT_TESTED",
          "healthy_adult_reference_feature_status":feature,
          "healthy_reference_top_celltype":bycell[name]["top_celltype"] if present else "NOT_ASSESSABLE",
          "healthy_reference_top_broad_class":bycell[name]["top_broad_class"] if present else "NOT_ASSESSABLE",
          "healthy_reference_top_detected_fraction":bycell[name]["top_detection_fraction"] if present else "",
          "healthy_reference_top_CPM":bycell[name]["top_pseudobulk_cpm"] if present else "",
          "healthy_reference_total_cells_detected":bycell[name]["total_detected_cells"] if present else "",
          "n_preselected_donor_pair_comparisons":len(pair_list),
          "n_descriptive_donor_QC_PASS":len(passing),
          "donor_comparison_scope":"ADULT_CONTROL_REFERENCE_ONLY",
          "candidate_causality":"NOT_ESTABLISHED",
          "stroke_patient_control_DGE":"NOT_PERFORMED",
          "BBJ_GWAS_effect_allele_reference":"NATIVE_8956_AUDIT_SEPARATELY",
          "GTEx_effect_allele_per_variant_validity":"OFFICIAL_ALT_EFFECT_BY_CATALOGUE; RAW_GENOTYPE_NOT_ATTESTED",
          "matched_LD_multisignal":"NOT_VERIFIED",
        })
    if len(out)!=9 or set(x["historic_symbol"] for x in out)!=set(NINE):
        raise ValueError("Nine-gene curated matrix incomplete")
    return out

def main():
    a=argparse.ArgumentParser()
    for arg in ("genes","human-summary","feature-status","donor-pairs","prior-grid","human-manifest","out-dir"):
        a.add_argument("--"+arg,required=True,type=Path)
    v=a.parse_args()
    files={"genes":v.genes,"cell":v.human_summary,
           "features":v.feature_status,"donor":v.donor_pairs,
           "grid":v.prior_grid,"human_manifest":v.human_manifest}
    inp={k:json.loads(p.read_text()) if k=="human_manifest" else read(p) for k,p in files.items()}
    result=audited(**inp)
    if v.out_dir.exists() and any(v.out_dir.iterdir()):
        raise FileExistsError("Refuse to overwrite results")
    v.out_dir.mkdir(parents=True,exist_ok=True)
    with (v.out_dir/"IS_NINE_GENE_QTL_PRIOR_HUMAN_DONOR_CONVERGENCE.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter="\t")
        w.writeheader();w.writerows(result)
    summary={
       "status":"DESCRIPTIVE_EXISTING_DATA_ONLY_NINE_GENE_MECHANISM_READINESS",
       "nine_historic_genes":[x["historic_symbol"] for x in result],
       "positional_gene_universe_kept":2225,
       "n_with_archived_GTEx_ABF":sum(x["n_BBJ_GTEx_v8_archived_QTL_assays"]>0 for x in result),
       "n_healthy_brain_reference_feature_present":sum(x["healthy_adult_reference_feature_status"]=="FEATURE_PRESENT" for x in result),
       "n_feature_not_in_matrix":sum(x["healthy_adult_reference_feature_status"]=="FEATURE_NOT_IN_MATRIX" for x in result),
       "n_high_baseline_H4_0p8":sum(x["H4_ge_0p8_original"] is True for x in result),
       "source_SHA256":{k:digest(p) for k,p in files.items()},
       "no_stroke_case_control_expression_test":True,
       "no_causal_gene_promotion":True,
    }
    (v.out_dir/"IS_NINE_GENE_MECHANISM_READINESS_MANIFEST.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("IS_NINE_GENE_DESCRIPTIVE_CONVERGENCE_PASS",json.dumps({k:summary[k] for k in ("n_with_archived_GTEx_ABF","n_healthy_brain_reference_feature_present","n_feature_not_in_matrix","n_high_baseline_H4_0p8")}),flush=True)
    for r in result:
        print("GENE",r["historic_symbol"],r["best_H4_original_p12_1e5"],
              r["healthy_reference_top_celltype"],r["healthy_adult_reference_feature_status"])
if __name__=="__main__":main()
