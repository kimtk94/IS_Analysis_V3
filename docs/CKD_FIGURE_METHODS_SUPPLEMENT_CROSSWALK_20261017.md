# CKD manuscript Figure and Methods crosswalk (2026-10-17)

This document maps the four panels to original server aggregates, evidence claims, and required Supplementary items. Figure rendering success is **software reproducibility**, not final causal inference validation.

| Figure | Input and method | Source/frozen contract | Permitted claim | Must not claim | Supplement |
|---|---|---|---|---|---|
| **Figure 1** Study workflow | UKB-PPP cis-pQTL, EUR kidney-trait GWAS, MR/coloc/SuSiE and EAS/tissue evidence; schematic | Integrated CKD nine-gene TSV, SHA256 `32fa0475461b206c48295833acfc3075bfb24d89c322b5f9be554b1883a2f3d3` | Overview of existing analysis stages and **reported nine-gene panel** | Quantified candidate attrition from full pQTL universe, unless separately source-audited | Full sample definitions, cis window, complete inclusion/exclusion counts and 9-gene identity |
| **Figure 2** EAS single-anchor Wald MR | Nine EAS eGFR anchor beta/SE; 95% CI = beta ±1.96SE | Same SHA-pinned nine-gene table | Reported estimate and uncertainty for each anchor | Korean longitudinal causal replication, direct drug effect, independent validation | Harmonized allele direction, F statistics, source N, multiple-testing family, weak IV/pleiotropy and overlap |
| **Figure 3** EUR ABF vs SuSiE | EUR pQTL vs eGFR regional ABF PP.H4 and SuSiE max PP.H4 | Same SHA-pinned table | Agreement/discordance of **source statistical model estimates** | Proved causal protein/gene; H4=0 as substitute for nonconverged SuSiE | ABF priors (p1,p2,p12), H3/H4, overlap SNP count, LD provenance, sensitivity, convergence, credible sets |
| **Figure 4** Kidney QTL and atlas | Kidney significant-only eQTL hit counts (meta/tubule/glomerulus), annotation labels | Same SHA-pinned table | Tissue/cell **availability and localization** | Direct cross-study effect-size comparisons, causal cell-type perturbation | Tissue source/sample sizes, detection/coverage, atlas version, QTL vs mere expression, published target overlap |

**Figure source scripts**: `scripts/plot_ckd_manuscript_fig1.py`, `scripts/plot_ckd_manuscript_fig2_fig4.py`, `scripts/plot_ckd_manuscript_fig3.py`. All require source SHA256, refuse nonempty output directories, and use source aggregate records rather than fabricated biological estimates. Figure 1 is a conceptual schematic; arrows are conceptual workflow, **not real candidate attrition counts**.

**Original source**: `/srv/is-analysis/results/ckd/integrated_evidence/CKD_9GENE_EVIDENCE_MATRIX.tsv`; source SHA256 `32fa0475461b206c48295833acfc3075bfb24d89c322b5f9be554b1883a2f3d3`. GitHub Actions runner fetches existing read-only server data and deletes temporary copies. Artifacts expire after 30 days.

**Primary literature alignment:** Won-lab dyslipidemia integrative target prioritization (*Cell Reports Medicine*, 2023), problematic alcohol use MR+coloc (*eBioMedicine*, 2025); Hirohama et al. kidney proteogenomics (*Nature Medicine*, 2025) and Sheng et al. kidney cell type (*Nature Genetics*, 2021).

### Outstanding science gates
- Numeric EAS MR output orientation/scale, F and phenotype coding; effect-axis is in **original model units** only.
- EUR QTL/GWAS allele matching, SNP counts, priors, LD reference and SuSiE convergence, especially UMOD/INHBC.
- Kidney tissue eQTL source significance/absence vs detection limits; dominant cell expression != regulatory colocalization.
- KoGES individual-level longitudinal validation remains **NOT COMPLETED**; do not imply otherwise.
- Before main merge, PR #34 must be reviewed for code/data reproducibility and claims; no blanket science PASS.
