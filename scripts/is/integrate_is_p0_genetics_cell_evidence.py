"""Reproducible genetic × human-reference status overlay.

This is a provenance-checked *evidence ledger*, NOT a causal gene rank.
Reads historical H4 ABF summary, coloc prior-reweighted diagnostics, R3_7
verified donor comparisons, and feature status. Writes into a fresh folder.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

FOCUS=("FGF5","CALHM2","SH3PXD2A","COL4A2","COL4A1","ALDH2","NEURL1","C4orf22","INA")
KEY_P12=(0.000001,0.00001,0.0001)


def load(p):
    with p.open(newline="",encoding="utf8") as f:
        return list(csv.DictReader(f,delimiter="\t"))


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(master,priors,donors,features,out):
    if out.exists() and any(out.iterdir()):
        raise FileExistsError("Refusing nonempty output")
    a=load(master)
    p=load(priors)
    d=load(donors)
    f=load(features)
    if len(a)!=646 or len(p)!=3230 or len(d)!=10 or len(f)!=9:
        raise ValueError("Unexpected original/conditioned/donor counts")
    if set(FOCUS)-set(x["gene"] for x in f):
        raise ValueError("Gene feature status missing")
    fs={x["gene"]:x["feature_status"] for x in f}
    if fs["FGF5"]!="FEATURE_NOT_IN_MATRIX":
        raise ValueError("FGF5 feature changed; must reassess")
    if any(x["model_note"]!="SUMMARY_ONLY_CONDITIONAL_PRIOR_REWEIGHT" for x in p):
        raise ValueError("Not a summary-only prior sensitivity provenance")
    if any(abs(float(x["old_p12"])-1e-5)>1e-12 for x in p):
        raise ValueError("Assumed original p12 unexpectedly changed")
    audit=[]
    for gene in FOCUS:
        rows=[x for x in a if x["gene_symbol"]==gene and x["status"]=="PASS"]
        if not rows:
            raise ValueError("Missing legacy gene in source ABF: "+gene)
        top=max(rows,key=lambda x:float(x["PP.H4"]))
        qs=[x for x in p if x["gene_symbol"]==gene and x["dataset_key"]==top["dataset_key"]
            and x["locus"]==top["locus"] and x["gene_base"]==top["gene_base"]]
        if len(qs)!=5:raise ValueError("Expected five conditional p12 values")
        prior={float(x["conditional_p12"]):x for x in qs}
        if set(prior)!=set((1e-6,3e-6,1e-5,3e-5,1e-4)):
            raise ValueError("Unrecognized prior grid")
        if abs(float(prior[1e-5]["PP_H4"])-float(top["PP.H4"]))>1e-10:
            raise ValueError("Legacy baseline posterior not conserved")
        pairs=[x for x in d if x["gene"]==gene]
        n_present=sum(x["qc_status"]=="DESCRIPTIVE_PAIRED_QC_PASS" for x in pairs)
        n_unanimous=sum(int(x["n_positive"])==int(x["n_paired_donors"]) for x in pairs)
        missing=fs[gene]=="FEATURE_NOT_IN_MATRIX"
        status="REFERENCE_NOT_ASSESSABLE" if missing else (
               "DONOR_PAIRED_DESCRIPTIVE_ONLY" if n_present else
               "FEATURE_PRESENT_BUT_NO_PREDECLARED_PAIR")
        audit.append(dict(
            gene=gene,locus=top["locus"],selected_tissue=top["dataset_key"],
            n_legacy_assays=len(rows),nsnps=int(top["nsnps"]),
            ABF_H4_baseline=float(top["PP.H4"]),ABF_H3_baseline=float(top["PP.H3"]),
            H4_over_H3_baseline=float(top["PP.H4"])/float(top["PP.H3"]) if float(top["PP.H3"])>0 else None,
            H4_p12_1e_minus_6=float(prior[1e-6]["PP_H4"]),
            H4_p12_1e_minus_5=float(prior[1e-5]["PP_H4"]),
            H4_p12_1e_minus_4=float(prior[1e-4]["PP_H4"]),
            reference_RNA_feature_status=fs[gene],
            donor_comparisons_qc_pass=n_present,
            donor_comparisons_unanimous_higher_CPM=n_unanimous,
            reference_cell_evidence_status=status,
            joint_mechanism_verdict="NOT_ESTABLISHED",
            model_source_status="ASSUMED_BASELINE_PRIORS_UNCONFIRMED",
            ancestry_LD_status="NOT_RECONFIRMED",
            independent_molecular_QTL_status="PENDING",
            disease_cell_QTL_status="PENDING"
        ))
    out.mkdir(parents=True,exist_ok=True)
    cols=list(audit[0])
    with (out/"IS_INTEGRATED_P0_GENE_EVIDENCE.tsv").open("w",newline="",encoding="utf8") as h:
        w=csv.DictWriter(h,fieldnames=cols,delimiter="\t")
        w.writeheader();w.writerows(audit)
    (out/"IS_INTEGRATED_P0_GENE_EVIDENCE.json").write_text(json.dumps({
        "schema":"IS_EVIDENCE_OVERLAY_V1",
        "claims":"PRIOR_SENSITIVITY_AND_HEALTHY_REFERENCE_ONLY",
        "legacy_coloc_tests":646,"conditional_prior_rows":3230,
        "donor_paired_comparisons":10,"feature_targets":9,
        "sources_SHA256":{str(q):digest(q) for q in (master,priors,donors,features)},
        "genes":audit
    },ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    md=["# IS genetic and reference evidence integration (exploratory)",
       "","**No new causal genes established.**",
       "",
       "Original four-locus GTEx ABF signals and R3_7 healthy donor expression",
       "are biologically compatible but not evidence for an allele-specific or",
       "disease-specific regulatory pathway. The p12 grid is summary-only,",
       "assuming that old p1=p2=1e-4 and p12=1e-5, not confirmed from the",
       "source SNP-level run log. The best H4 tissue is selected post hoc.",
       "",
       "| Gene | Best historical tissue | Baseline H4 | H4 (p12=1e-6) | H4 (p12=1e-4) | Donor comparisons with higher CPM | Cell evidence |",
       "|---|---|---:|---:|---:|---:|---|"]
    for x in audit:
        md.append(f'| {x["gene"]} | {x["selected_tissue"]} | {x["ABF_H4_baseline"]:.3f} | {x["H4_p12_1e_minus_6"]:.3f} | {x["H4_p12_1e_minus_4"]:.3f} | {x["donor_comparisons_unanimous_higher_CPM"]}/{x["donor_comparisons_qc_pass"]} | {x["reference_cell_evidence_status"]} |')
    md+=["","## Prepublication blockers",
          "",
          "1. Verify the original ABF prior settings and region-level H0–H4 from archived scripts, run logs, and full variant-level inputs.",
          "2. Harmonize alleles, chromosome build, GWAS/EAS QTL SNP overlap, sample sizes, LD ancestry and cis-window definitions.",
          "3. Run genuine SNP-level prior sensitivity and multi-signal SuSiE on matched QTL/GWAS signals; summary reweight is only a conditional diagnostic.",
          "4. Obtain independent molecular QTL and disease cell-type QTL or validated functional perturbation if claiming a causal regulatory mechanism.",
          "5. Correct for tissue/gene selection and multiplicity before manuscript statements.",
          "6. Preserve 80 provisional regions / 2,225 positional genes; original 646 tests are a selected four-locus analysis, not an expanded genome-wide QTL assay.",
          "",
          "## FGF5 caution",
          "",
          "FGF5 has high historical cerebellar H4 under baseline p12 but falls substantially",
          "under conservative p12. The gene is absent as an RNA feature in this",
          "particular healthy temporal-lobe processed Seurat object (not biological zero),",
          "so this R3_7 dataset cannot supply independent FGF5 cell-localization evidence.",
          "",
          "## SH3PXD2A and COL4A2 caution",
          "",
          "Positive paired donor expression contrast in normal tissue is a reference-cell",
          "localization observation; it does not close the eQTL–GWAS causal chain.",
          ""]
    (out/"IS_INTEGRATED_P0_INTERPRETATION.md").write_text("\n".join(md)+"\n",encoding="utf8")
    print("IS_P0_GENETIC_CELL_OVERLAY=PASS")
    print("FOCUS_GENES=",len(audit),"GENES_WITHOUT_REGULATORY_CAUSAL_PROOF=",sum(x["joint_mechanism_verdict"]=="NOT_ESTABLISHED" for x in audit))
    for x in audit[:6]:
        print(x["gene"],f'{x["ABF_H4_baseline"]:.3f}',
              f'{x["H4_p12_1e_minus_6"]:.3f}',
              x["reference_cell_evidence_status"])
    return audit


def main():
    ap=argparse.ArgumentParser()
    for flag in ("master","priors","donors","features","out"):
        ap.add_argument("--"+flag,type=Path,required=True)
    x=ap.parse_args()
    run(x.master,x.priors,x.donors,x.features,x.out)


if __name__=="__main__":
    main()
