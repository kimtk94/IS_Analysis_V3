#!/usr/bin/env python3
"""Read-only IS backlog inventory using ONLY already acquired local research outputs.

Never downloads data, changes original files or promotes an association
to an independent replication, a causal gene, or validated fine-mapping.
All counts are derived from the actual TSV/JSON inputs supplied.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

FILES = {
    "expanded": "stage5_functional/broad_discovery_v2_ancestry/IS_ANCESTRY_GENE_EVIDENCE_V3.tsv",
    "region_summary": "stage5_functional/broad_discovery_v2_ancestry/IS_ALL_GENE_UNIVERSE_SUMMARY.json",
    "coloc": "stage5_functional/phase9c_convergence/COLOC_ABF_MASTER_ANNOTATED_V2.tsv",
    "coloc_replay": "audits/legacy_coloc_646_batch_cache30_20261010_v1/IS_LEGACY_646_SNP_REPLAY_STATUS.tsv",
    "priority8": "audits/is_priority8_source_verified_20261010_v2/IS_PRIORITY8_SOURCE_VERIFIED.tsv",
    "human_manifest": "stage5_functional/phase9f_e_human_r3_7_donor_import_20261010/HUMAN_RUN_MANIFEST.json",
    "donor_audit": "audits/phase9f_r3_7_donor_reaudit_20261010_v1/R3_7_DONOR_REAUDIT.json",
    "p12_grid": "audits/legacy_coloc_p12_sensitivity_20261010_v1/IS_LEGACY_TOP_GENE_P12_GRID.tsv",
    "crossancestry": "stage4_cross_eas/GIGASTROKE_BBJ_HARMONIZED_VARIANTS.tsv",
}
ACTIONS = [
    ("P0","EVIDENCE_UNIVERSE","80 regions / 2,225 gene IDs: preserve full positional universe; link EAS/EUR, QTL, literature and cell-type evidence with missing/not-tested flags"),
    ("P0","COLOC_PROVENANCE_RECOVERY","Recover locally archived per-SNP inputs for 616 MISSING_INPUT pairs before expanding exact SNP-level ABF replay"),
    ("P1","COLOC_SENSITIVITY","Independently audit 30 PASS SNP-level replays and selected 8 priority gene-tissue pairs against source hashes, sample-N and prior p12 grid"),
    ("P1","DONOR_LEVEL_REFERENCE","Reassess healthy adult brain reference cell-type pseudobulk, donor coverage and feature presence; not stroke differential expression"),
    ("P1","CROSS_ANCESTRY_SENSITIVITY","Check existing BBJ/GIGASTROKE allele harmonization, direction, ancestry heterogeneity and overlapping-cohort caveats"),
    ("P1","COMPETING_MECHANISMS","Build gene-wise alternative mechanism/evidence gaps for FGF5/CALHM2/NEURL1/C4orf22/INA/SH3PXD2A/COL4A1/COL4A2/ALDH2"),
    ("P2","PAPER_AND_DASHBOARD","Produce Figure-ready reproducible outputs and past-study comparison; distinguish computational completion from scientific verification"),
]
def read_tsv(path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b""):h.update(chunk)
    return h.hexdigest()

def audit(results_root):
    paths={k:results_root/rel for k,rel in FILES.items()}
    absent=[k for k,p in paths.items() if not p.is_file()]
    if absent:raise FileNotFoundError(f"Missing required pre-existing IS outputs: {absent}")
    expanded=read_tsv(paths["expanded"])
    regions=json.loads(paths["region_summary"].read_text())
    coloc=read_tsv(paths["coloc"])
    replay=read_tsv(paths["coloc_replay"])
    priority=read_tsv(paths["priority8"])
    human=json.loads(paths["human_manifest"].read_text())
    donor=json.loads(paths["donor_audit"].read_text())
    p12=read_tsv(paths["p12_grid"])

    ids={r["gene_id_stable"] for r in expanded if r.get("gene_id_stable")}
    locus_ids={(r["locus"],r["dataset_key"],r["gene_base"]) for r in coloc}
    replay_ids={(r["locus"],r["dataset_key"],r["gene_base"]) for r in replay}
    if len(expanded)!=2426 or len(ids)!=2225 or regions["region_gene_associations"]!=2425:
        raise ValueError("Expanded universe drift: inspect gene/region provenance")
    if len(coloc)!=646 or len(locus_ids)!=646 or len(replay)!=646 or locus_ids!=replay_ids:
        raise ValueError("Legacy coloc versus SNP replay key mismatch")
    replay_counts=Counter(r["status"] for r in replay)
    if any(x not in ("PASS","MISSING_INPUT") for x in replay_counts):
        raise ValueError("Unexpected SNP replay status needs scientific review")
    if len(priority)!=8:
        raise ValueError("Priority gene-tissue source replay count drift")
    for r in priority:
        if r["status"]!="ORIGINAL_FIVE_HYPOTHESES_SNP_REPLAY_PASS_NO_CAUSAL_CLAIM":
            raise ValueError("Unexpected priority 8 source replay classification")
    if human["donor_analysis_unit"]!="Patient" or donor["status"]!="DESCRIPTIVE_RECALC_PASS_NO_CAUSAL_INFERENCE":
        raise ValueError("Unexpected donor-scientific-status drift")

    return {
        "status":"READ_ONLY_LOCAL_DATA_AUDIT_NO_NEW_DB",
        "gene_universe":{"region_count":regions["regions"],
            "unique_positional_gene_ids":len(ids),
            "gene_region_pairs":regions["region_gene_associations"],
            "expanded_evidence_rows":len(expanded),
            "ancestry_rows":dict(Counter(x["ancestry"] or "LEGACY_ANCHOR_ONLY" for x in expanded)),
            "eqtl_status":dict(Counter(x["eqtl_status"] for x in expanded)),
            "sqtl_status":dict(Counter(x["sqtl_status"] for x in expanded)),
            "pqtl_status":dict(Counter(x["pqtl_status"] for x in expanded)),
            "causal_genes_established":0,
            "note":"Expanded table has a legacy-anchor-only row; regional associations=2425, evidence rows=2426"},
        "coloc":{"original_abf_assays":len(coloc),
            "SNP_replay_status":dict(replay_counts),
            "priority8_source_verified":len(priority),
            "priority_genes":[r["gene"] for r in priority],
            "prior_sensitivity_records_for_selected_genes":len(p12),
            "interpretation":"No cohort-matched molecular QTL LD verified; ABF posterior is not a causal gene guarantee"},
        "human_reference":{"donors":donor["donor_n"],"cells":donor["n_cells"],
            "pseudobulk_rows":donor["pseudobulk_n"],
            "paired_comparisons":donor["paired_n"],"paired_QC_PASS":donor["paired_pass_n"],
            "source_reference_type":human["reference_type"],
            "clinical_stroke_differential_expression":False},
        "tasks":[{"priority":prio,"workstream":name,"ready_action":desc}
                 for prio,name,desc in ACTIONS],
        "excluded_external_dependencies":["CKB encrypted summary-stat ZIP key",
             "TPMI restricted 433.21 GWAS download",
             "Japanese cohort-matched GWAS LD",
             "new paid/controlled individual-level data"],
        "scientific_guards":["Preserve 2225-gene universe: rank evidence, never hard-exclude positional candidates",
            "No reuse of BBJ and GIGASTROKE EAS as independent replication",
            "Do not label a human healthy reference as stroke lesion-vs-control differential expression",
            "Explicitly separate MISSING_INPUT from negative coloc",
            "ALDH2 alcohol fine-mapping BLOCKED; ADH1B EXPLORATORY"],
        "source_sha256":{k:sha256(p) for k,p in paths.items()},
        "source_paths":{k:str(p) for k,p in paths.items()}
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--results-root",type=Path,default=Path("/srv/is-analysis/results/is"))
    p.add_argument("--out-dir",required=True,type=Path)
    a=p.parse_args()
    result=audit(a.results_root)
    if a.out_dir.exists() and any(a.out_dir.iterdir()):
        raise FileExistsError("Refuse to overwrite previous audit")
    a.out_dir.mkdir(parents=True,exist_ok=True)
    (a.out_dir/"IS_NO_NEW_DB_BACKLOG_AUDIT.json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    with (a.out_dir/"IS_NO_NEW_DB_WORKSTREAMS.tsv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=["priority","workstream","ready_action"],delimiter="\t")
        writer.writeheader()
        writer.writerows(result["tasks"])
    print("IS_NO_NEW_DB_BACKLOG_PASS")
    print("REGIONS",result["gene_universe"]["region_count"],
          "GENES",result["gene_universe"]["unique_positional_gene_ids"],
          "PAIRS",result["gene_universe"]["gene_region_pairs"])
    print("COLOC_REPLAY",json.dumps(result["coloc"]["SNP_replay_status"],sort_keys=True))
    print("DONORS",result["human_reference"]["donors"],
          "CELLS",result["human_reference"]["cells"],
          "PAIRED_QC_PASS",result["human_reference"]["paired_QC_PASS"])
    print("BACKLOG",str(a.out_dir))

if __name__=="__main__":
    main()
