#!/usr/bin/env python3
"""Gate B: GSE313736 post-exposure cell-state eligibility and author DEG-list provenance.

Read ONLY official public cell metadata or Series SOFT and frozen aggregate
edgeR results. Never publish donor IDs, barcodes, single-cell expression,
individual genotypes, demographics, or per-donor count vectors.

This is a METADATA sensitivity gate; it does not refit edgeR or establish
that gene effects survive an alternative post-treatment selection threshold.
"""
import argparse,csv,gzip,io,json,math,re,tempfile
from collections import Counter,defaultdict
from pathlib import Path
from statistics import median

from ocular_p1_integrity import BASE,download
from ocular_p0_audit import fetch,root

METADATA_URL=BASE+"series/GSE313nnn/GSE313736/suppl/GSE313736_HTMC_scRNA_cell_metadata.tsv.gz"
CELL_MODES={
 "BeamA_score0.50":("BeamA",0.50),
 "BeamA_score0.75":("BeamA",0.75),
 "JCT_score0.50":("JCT",0.50),
 "JCT_score0.75":("JCT",0.75),
 "scpred_TM_Fibro_unscored":("scpred",None),
}
MIN_CELLS=(30,50)
def pool_name(lib):
    m=re.search(r"pool[-_](?:dex)?([A-I])(?:_|$)",lib)
    return m.group(1) if m else None
def output_json(out,name,obj):
    (out/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf8")
def output_tsv(out,name,rows):
    if not rows:raise ValueError("EMPTY_AGGREGATE_OUTPUT")
    with (out/name).open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader();w.writerows(rows)
def metaaudit(out):
    with tempfile.TemporaryDirectory(prefix="ocular-gateb-meta-") as tmp:
        fp=Path(tmp)/"cellmeta.gz"
        digest=download(METADATA_URL,fp,12*1024*1024)
        source=defaultdict(Counter) # anonymous in-memory biological donor x treatment
        all_counts=defaultdict(Counter)
        pools={}
        labels=Counter()
        nrows=0
        with gzip.open(fp,"rt",encoding="utf8",newline="") as f:
            rdr=csv.DictReader(f,delimiter="\t")
            required={"donor_id_demuxlet","DEX Treatment","library_id","droplet_type",
                      "celltypist_chen2025_fine","celltypist_chen2025_fine_score",
                      "scpred_vanzyl2022"}
            if not required.issubset(rdr.fieldnames or []):
                raise ValueError("CELL_METADATA_REQUIRED_ANNOTATION_COLUMNS_MISSING")
            for r in rdr:
                nrows+=1
                if r["droplet_type"]!="SNG":continue
                donor=r["donor_id_demuxlet"].strip()
                arm=r["DEX Treatment"].strip().lower()
                pool=pool_name(r["library_id"].strip())
                if not donor or arm not in ("dex","no-dex") or pool is None:
                    raise ValueError("UNEXPECTED_DONOR_ARM_OR_POOL")
                if donor in pools and pools[donor]!=pool:
                    raise ValueError("DONOR_OCCURS_IN_MULTIPLE_POOLS")
                pools[donor]=pool
                all_counts[donor][arm]+=1
                label=r["celltypist_chen2025_fine"].strip()
                alternate=r["scpred_vanzyl2022"].strip()
                try:score=float(r["celltypist_chen2025_fine_score"])
                except (TypeError,ValueError):score=float("nan")
                for k,(want,cutoff) in CELL_MODES.items():
                    selected=(alternate=="TM_Fibro" if want=="scpred" else label==want and math.isfinite(score) and score>=cutoff)
                    if selected:
                        source[(donor,arm)][k]+=1
                        labels[(arm,k)]+=1
        if nrows!=169398 or len(pools)!=59 or len({x[0] for x in source})>59:
            raise ValueError("EXPECTED_PUBLIC_META_SNG_169398_AND_59_DONORS")
        if set(all_counts)!={d for d in pools}:
            raise ValueError("UNEXPECTED_SINGLETON_ONLY_DONORS")
        result=[]
        perpool=[]
        donor_sets={}
        arm_total={"dex":sum(c["dex"] for c in all_counts.values()),
                   "no-dex":sum(c["no-dex"] for c in all_counts.values())}
        for mode in CELL_MODES:
            for cells in MIN_CELLS:
                eligible={d for d in pools
                   if source[(d,"dex")][mode]>=cells and source[(d,"no-dex")][mode]>=cells}
                donor_sets[(mode,cells)]=eligible
                retained={a:sum(source[(d,a)][mode] for d in eligible) for a in arm_total}
                fractions=[]
                for d in eligible:
                    na=all_counts[d]["dex"];nb=all_counts[d]["no-dex"]
                    if na>0 and nb>0:
                        fractions.append(source[(d,"dex")][mode]/na - source[(d,"no-dex")][mode]/nb)
                record={"annotation_mode":mode,"minimum_cells_per_arm":cells,
                   "eligible_matched_donors":len(eligible),
                   "selected_metadata_dex_cells":retained["dex"],
                   "selected_metadata_vehicle_cells":retained["no-dex"],
                   "selected_cells_dex_all_donors":labels[("dex",mode)],
                   "selected_cells_vehicle_all_donors":labels[("no-dex",mode)],
                   "median_donor_fraction_difference_dex_minus_vehicle":round(median(fractions),7) if fractions else "",
                   "anonymous_donor_IDs_released":False,
                   "raw_count_barcode_joins_performed":False,
                   "effect_direction_QC_performed":False}
                result.append(record)
                for p in "ABCDEFGHI":
                    s={d for d in eligible if pools[d]==p}
                    perpool.append({"annotation_mode":mode,"minimum_cells_per_arm":cells,"pool":p,
                       "eligible_pairs":len(s),"total_metadata_donors_in_pool":sum(x==p for x in pools.values())})
        if len(donor_sets[("BeamA_score0.75",30)])!=57:
            raise ValueError("P4_EXPECTED_57_STRICT_PAIRS_DRIFT")
        if not donor_sets[("BeamA_score0.75",30)].issubset(donor_sets[("BeamA_score0.50",30)]):
            raise ValueError("NESTED_BEAMA_THRESHOLD_BROKEN")
        if not donor_sets[("BeamA_score0.75",50)].issubset(donor_sets[("BeamA_score0.75",30)]):
            raise ValueError("NESTED_MINIMUM_CELL_THRESHOLD_BROKEN")
        overlap=donor_sets[("BeamA_score0.75",30)]&donor_sets[("scpred_TM_Fibro_unscored",30)]
        summary={"status":"PASS_METADATA_ONLY_CELL_STATE_SELECTION_GRID",
            "study":"GSE313736","cell_metadata_sha256":digest["sha256"],
            "metadata_cells":nrows,"biological_donors_in_both_arms":len(pools),
            "per_arm_total_sng_cells":arm_total,
            "n_modes":len(CELL_MODES),"n_cutoff_grid":len(result),
            "strict_57_BeamA_score_075_min30":57,
            "score050_additional_donors_vs_score075":len(donor_sets[("BeamA_score0.50",30)]-donor_sets[("BeamA_score0.75",30)]),
            "alternative_scpred_TM_Fibro_overlap_strict_BeamA_donors":len(overlap),
            "score050_075_subsets_nesting":"PASS",
            "cell_annotation_was_applied_after_DEX":True,
            "score_change_does_not_establish_CPM_or_edgeR_gene_stability":True,
            "per_donor_vectors_and_donor_IDs_not_published":True,
            "clinical_delta_IOP_measured":False}
        output_tsv(out,"GateB_cell_state_10_threshold_metadata_grid.tsv",result)
        output_tsv(out,"GateB_cell_state_per_pool_grid.tsv",perpool)
        output_json(out,"GateB_cell_state_selection_metadata_gate.json",summary)
        print("GATEB_METADATA",json.dumps(summary,ensure_ascii=False),flush=True)

def authorsaudit(out,primary_file):
    raw,tx=fetch(root("GSE313736")+"soft/GSE313736_family.soft.gz",limit=9*1024*1024)
    series_files=[]
    per_sample_files=[]
    title=""
    abstract=""
    relation=[]
    with gzip.open(io.BytesIO(raw),"rt",encoding="utf8",errors="replace") as f:
        for s in f:
            line=s.rstrip("\r\n")
            if line.startswith("!Series_supplementary_file"):
                series_files.append(line.split("=",1)[-1].strip())
            elif line.startswith("!Sample_supplementary_file"):
                per_sample_files.append(line.split("=",1)[-1].strip())
            elif line.startswith("!Series_title"):
                title=line.split("=",1)[-1].strip()
            elif line.startswith("!Series_summary"):
                abstract+=(line.split("=",1)[-1].strip()+" ")
            elif line.startswith("!Series_relation"):
                relation.append(line.split("=",1)[-1].strip())
    with gzip.open(primary_file,"rt",encoding="utf8",newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        records=list(reader)
    if len(records)!=16050:raise ValueError("P5B_EXPECTED_16050_GENES")
    high=[]
    for r in records:
        logfc=float(r["logFC"]);fdr=float(r["FDR"])
        if fdr<.05 and abs(logfc)>=3.:
            high.append({"gene_id":r["gene_id"],"gene_symbol":r["gene_symbol"],
                         "our_48h_log2FC":round(logfc,6),"our_edgeR_BH_FDR":fdr,
                         "direction":"UP" if logfc>0 else "DOWN",
                         "original_paper_individual_gene_overlap":"NOT_DETERMINED_NO_AUTHORS_LIST"})
    up=sum(x["direction"]=="UP" for x in high)
    down=sum(x["direction"]=="DOWN" for x in high)
    if (len(high),up,down)!=(83,67,16):
        raise ValueError("P5_83_LARGE_EFFECT_SOURCE_CHANGED "+str((len(high),up,down)))
    # Do not infer that a GEO file list (or article abstract) contains a full
    # authors' gene-level DEG table, and never fabricate a gene-name overlap.
    output_tsv(out,"GateB_our_83_high_effect_genes_NOT_author_overlap.tsv",sorted(high,key=lambda x:x["our_48h_log2FC"],reverse=True))
    summary={"status":"PASS_OFFICIAL_GEO_SOURCE_AUDIT_HOLD_AUTHOR_LIST_COMPARISON",
       "geo_study":"GSE313736","series_title":title,
       "soft_transfer_sha256":tx.get("sha256"),"series_supplementary_file_count":len(series_files),
       "sample_supplementary_file_entry_count":len(per_sample_files),
       "series_supplement_names_only":[x.rsplit("/",1)[-1] for x in series_files],
       "author_abstract_claim_up":63,"author_abstract_claim_down":17,
       "author_abstract_downregulated_threshold_text_AMBIGUOUS":True,
       "our_high_effect_up_BH_FDR_lt_05_log2FC_ge_3":up,
       "our_high_effect_down_BH_FDR_lt_05_log2FC_le_minus3":down,
       "author_vs_our_numeric_high_effect_net_difference":len(high)-80,
       "author_exact_signed_gene_list_accessed":False,
       "gene_level_author_vs_our_overlap_n":"NOT_VERIFIABLE",
       "do_not_treat_83_vs_80_as_three_discordant_genes":True,
       "comparison_blockers":["authors 80 signed gene IDs not available in verified SOFT listing",
                              "original author analysis donor/cell threshold/count protocol not matched",
                              "GEO abstract downregulated threshold '<3' may be a typo; cannot assume '<-3'"],
       "requires_original_author_supplement_or_communication":True}
    output_json(out,"GateB_authors_63_17_exact_list_provenance_gate.json",summary)
    print("GATEB_AUTHOR",json.dumps(summary,ensure_ascii=False),flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",required=True,choices=("eligibility","author"))
    ap.add_argument("--primary",type=Path)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    if a.mode=="eligibility":metaaudit(a.out)
    else:
        if a.primary is None:ap.error("author audit requires --primary")
        authorsaudit(a.out,a.primary)
if __name__=="__main__":main()
