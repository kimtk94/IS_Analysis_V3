#!/usr/bin/env python3
"""Complete-universe, noncausal IS 2225-gene evidence readiness audit.

Priority means existing-data work feasibility only, never gene causality.
Every stable Ensembl ID is retained, including rows without QTL evidence.
"""
import argparse,csv,hashlib,json
from collections import Counter
from pathlib import Path

def load(path):
    with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))

def run(gene_path, region_path, manifest_path,out):
    genes=load(gene_path);regions=load(region_path);m=json.loads(manifest_path.read_text())
    if not (len(genes)==m["stable_gene_ids"]==2225 and len(regions)==m["gene_region_pairs"]==2425):
        raise ValueError("source matrix count does not agree with manifest")
    ids=[x["gene_id_stable"] for x in genes]
    if len(set(ids))!=len(ids) or any(not s.startswith("ENSG") for s in ids):raise ValueError("invalid stable ID universe")
    support=Counter(x["gene_id_stable"] for x in regions)
    if set(support)!=set(ids):raise ValueError("Orphan or missing regional gene")
    if any(int(g["n_positional_regions"])!=support[g["gene_id_stable"]] for g in genes):
        raise ValueError("gene region multiplicity not reconciled")
    cases=[]
    for g in genes:
        def n(key):return int(g[key] or 0)
        p=float(g["minimum_region_lead_p"]) if g["minimum_region_lead_p"] else None
        if p is not None and not 0<p<=1:raise ValueError("Invalid GWAS p")
        eas=n("n_EAS_Japanese_regions");eur=n("n_EUR_regions")
        replay=n("archived_EAS_GTEx_SNP_replay_PASS")
        missing=n("archived_EAS_GTEx_SNP_replay_MISSING_INPUT")
        h=n("human_reference_donor_paired_QC_PASS")
        lite=n("linked_literature_gene_locus_records")
        pr=n("SNP_replayed_priority8_evidence")
        tier=("A_EXISTING_SNP_REPLAY_READY" if replay or pr else
              "B_RECOVER_ARCHIVED_ABF_INPUT" if missing else
              "C_EAS_GWAS_POSITIONAL_NO_QTL" if eas else
              "D_EUR_GWAS_POSITIONAL_NO_QTL")
        # This ranking is a workflow sequence, NOT evidence of mediation.
        score=(40 if replay or pr else 25 if missing else 10 if eas else 0)
        score+=min(eas,2)*5+min(eur,2)*2
        score+=8 if h else 0
        score+=min(lite,3)*3
        score+=min(n("archived_EAS_GTEx_ABF_assays"),3)
        row={"gene_id_stable":g["gene_id_stable"],
             "symbol":g["canonical_GENCODE_v19_gene_symbol"],
             "region_ids":g["region_ids"],
             "workflow_priority_tier":tier,
             "workflow_priority_score_NONCAUSAL":score,
             "ancestry_context":"EAS_AND_EUR" if eas and eur else "EAS_ONLY" if eas else "EUR_ONLY",
             "n_EAS_regions":eas,"n_EUR_regions":eur,
             "archived_ABF_assays":n("archived_EAS_GTEx_ABF_assays"),
             "archived_replay_PASS":replay,"archived_replay_MISSING_INPUT":missing,
             "human_reference_paired_QC_PASS":h,
             "linked_literature_rows":lite,
             "recommended_next_action":(
               "CHECK_ALREADY_REPLAYED_QTL_SNP_SIGNAL_PRIOR_AND_LD" if tier.startswith("A_") else
               "RESTORE_EXISTING_NOMINAL_QTL_SNP_INPUT_WITH_PROVENANCE" if tier.startswith("B_") else
               "REVIEW_EAS_LOCUS_QTL_SOURCE_COVERAGE_FIRST" if tier.startswith("C_") else
               "WAIT_FOR_VALIDATED_EUR_ANCESTRY_SPECIFIC_LD_AND_QTL"),
             "direct_causal_gene_status":"NOT_ESTABLISHED",
             "not_assessed_is_not_negative":True}
        cases.append(row)
    cases.sort(key=lambda r:(-r["workflow_priority_score_NONCAUSAL"],r["gene_id_stable"]))
    tier_counts=Counter(x["workflow_priority_tier"] for x in cases)
    for a,b in {"A_EXISTING_SNP_REPLAY_READY":14,"B_RECOVER_ARCHIVED_ABF_INPUT":29}.items():
        if tier_counts[a]!=b:raise ValueError(f"readiness category drift {a} {tier_counts[a]}")
    # Avoid G0022 dominating next steps: explicit secondary broad work queue.
    other=[x for x in cases if "G0022" not in x["region_ids"]]
    top=other[:40]
    out.mkdir(exist_ok=True,parents=True)
    def tsv(name,rows):
        with (out/name).open("w",newline="") as f:
            wr=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
            wr.writeheader();wr.writerows(rows)
    tsv("IS_2225_NONCAUSAL_WORKFLOW_PRIORITY.tsv",cases)
    tsv("IS_TOP40_NON_G0022_EXISTING_DATA_QUEUE.tsv",top)
    report={"audit":"IS_2225_COMPLETE_NONCAUSAL_EXISTING_DATA_PRIORITY_V1",
        "total_source_unique_genes":len(cases),
        "total_source_gene_region_pairs":len(regions),
        "source_ancestry_pairs":m["ancestry_pairs"],
        "tier_counts":dict(tier_counts),
        "non_G0022_genes":len(other),
        "TOP40_non_G0022_workflow_genes":[{"symbol":x["symbol"],"id":x["gene_id_stable"],"tier":x["workflow_priority_tier"],"score":x["workflow_priority_score_NONCAUSAL"]} for x in top],
        "source_manifest_sha256":hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "source_gene_matrix_sha256":hashlib.sha256(gene_path.read_bytes()).hexdigest(),
        "source_region_matrix_sha256":hashlib.sha256(region_path.read_bytes()).hexdigest(),
        "causal_genes_established":0,
        "research_wide_hard_candidate_filter_applied":False,
        "scientific_caveat":"Workflow scores measure audit/action readiness and are not Bayesian causal priors, gene biology importance or independent GWAS discoveries. Source ABF H4 requires overlap and LD verification."}
    (out/"IS_2225_WORKFLOW_PRIORITY_MANIFEST.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    return report

if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--root",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    p=a.parse_args()
    print(json.dumps(run(p.root/"IS_2225_GENE_EVIDENCE_MATRIX.tsv",
      p.root/"IS_2425_GENE_REGION_EVIDENCE.tsv",
      p.root/"IS_GENE_EVIDENCE_MATRIX_MANIFEST.json",p.out),indent=2))
