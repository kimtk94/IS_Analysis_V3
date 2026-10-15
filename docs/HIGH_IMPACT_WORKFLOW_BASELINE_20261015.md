# High-impact journal workflow baseline — IS / CKD (2026-10-15)

**Policy:** Match the **question-relevant analytic workflow** of representative Nature / Nature Medicine / Nature Genetics research papers. Do not add new analytic layers after coverage is supported by audited outputs. A journal's JIF is a journal-level metric; it is not a substitute for methodologic fitness, independent validation, or data integrity.

## Journal-level benchmark

Publisher-reported 2025 JIF (not the article's historical IF): Nature 56.1; Nature Medicine 52.5; Nature Genetics 25.5.
Source: https://www.nature.com/nature-portfolio/about-journals/journal-metrics

## Papers and exact roles

| Research | Paper | DOI | Match only |
|---|---|---|---|
| IS | Mishra et al., Nature (2022), *Stroke genetics informs drug discovery and risk prediction across ancestries* | 10.1038/s41586-022-05165-3 | Stroke GWAS signals, cross-ancestry comparisons, SuSiE fine-mapping, QTL/TWAS/PWAS-supported gene prioritization, relevant cell annotation |
| IS | Wälchli et al., Nature (2024), *Single-cell atlas of the human brain vasculature across development, adulthood and disease* | 10.1038/s41586-024-07493-y | Reuse published brain endothelial/perivascular cell annotations; not rebuild an atlas |
| CKD | Hirohama et al., Nature Medicine (2025), *The proteogenomic landscape of the human kidney and implications for cardio-kidney-metabolic health* | 10.1038/s41591-025-03872-8 | GWAS + tissue/plasma pQTL/eQTL; MR; coloc; integrated kidney protein/tissue evidence |
| CKD | Sheng et al., Nature Genetics (2021), *Mapping the genetic architecture of human traits to cell types in the kidney identifies mechanisms of disease and potential treatments* | 10.1038/s41588-021-00909-9 | Assign GWAS/QTL genes to biologically appropriate kidney cell types using existing published reference data |

The reference papers represent different primary research questions and data collection designs. Their entire experiments do **not** need to be copied. For example, Mishra et al. also investigated drug target repositioning and polygenic prediction, whereas our IS question is locus-to-molecular mechanism. Hirohama et al. generated primary kidney multi-omic data, whereas our CKD work can reuse publicly available summary-level results and supplements.

## IS — minimum completed workflow

1. **GWAS / candidate universe:** harmonize ancestry-specific BBJ / published stroke GWAS by assembly, effect allele, phenotype, and P; run a reproducible **full eligible-locus screen** before displaying the four priority genes. Retain non-promoted loci in a secondary evidence table. The highlighted four are not an initial discovery restriction.
2. **Cross-ancestry check:** for measurable loci compare EAS vs EUR direction, effect, and locus patterns. SuSiE 95% credible sets when trustworthy ancestry-matched LD is available; report 1KG / JPT reference-based results as sensitivity analyses, not a BBJ cohort-matched fine-map. Report nonconvergence.
3. **Molecular prioritization:** tissue-relevant eQTL, pQTL and colocalization (including multi-signal analysis where warranted); for a TWAS-supported claim, inspect the Nature stroke paper's transcriptome-wide correction and PP4 >=0.75 criterion. Avoid upgrading a candidate from PP4 alone. Coding ALDH2/rs671 is a distinct coding/protein/metabolic branch; do not require positive bulk eQTL for it.
4. **Cell context:** annotate prioritized candidates in **published** brain-vascular sc/snRNA and, only where needed, published regulatory/ATAC evidence. Distinguish expression-only evidence from allele-to-regulatory-element functional proof.
5. **Claim audit / figures:** 1 study flow; 2 EAS/EUR GWAS + representative credible-set locus; 3 coloc/MR evidence matrix with failure labels; 4 annotated cell evidence + prioritized mechanism summary. State the limitation of unavailable BBJ LD.
6. **STOP** once these results, QC/provenance, failure states, and figures are auditable. No mandatory new MENTR model, de novo atlas, drug repurposing, PGS, or wet-lab experiments.

**Current evidence (documentation-level, not an execution PASS):**
- `docs/IS_MASTER.md`: BBJ GWAS, FGF5/ALDH2/SH3PXD2A/COL4A2 branches and ABF/SuSiE interpretation are documented.
- **PARTIAL**: verify that discovery screening covers all eligible signals, not only four finalized genes.
- **PARTIAL**: source-matched BBJ LD unavailable; report surrogate LD sensitivity honestly.
- **PARTIAL**: cell-context and ATAC-to-variant link evidence are Phase 10 targets.
- Validate all of the above with **dated server outputs** before asserting COMPLETE.

## CKD — minimum completed workflow

1. **Outcome + exposure QC:** harmonize UKB-PPP cis-pQTL instruments and eGFRcrea primary GWAS; supporting CKD/BUN/UACR and EAS outcome direction; perform allele, build, EAF, LD and weak-instrument checks.
2. **Causal / shared-signal inference:** screen with appropriate MR estimator for available SNP count; report MR estimate + confidence interval and multiplicity control; regional coloc.abf; use coloc.susie when multiple signals are present and reliable LD / convergence are available. Distinguish PP4 support from proven causality.
3. **Kidney plausibility:** published kidney pQTL/eQTL and corresponding supplement-based candidate overlaps, plus kidney sc/snRNA cell-localization annotation. Do not treat absence from significant-only eQTL hit files as a measured null.
4. **Independent support:** compare EUR and EAS effect directions; use KoGES individual-level longitudinal validation only when access exists and the resulting claim requires that validation. If controlled-access data remain unavailable, state **NOT VALIDATED IN KOGES** rather than blocking the public-data thesis.
5. **Claim audit / figures:** 1 population and instrument flow; 2 proteome-wide MR / candidate table; 3 representative colocalization and locus evidence; 4 kidney tissue/cell contextualization with cross-ancestry evidence and stated limitations.
6. **STOP** once auditable results cover the above. No mandatory new tissue proteomics, new scRNA experiments, model stacking, broad extra drug screens, or forcing SuSiE convergence with arbitrarily many iterations.

**Current evidence (documentation-level, not an execution PASS):**
- `docs/ckd_stage1_mr.md`, `docs/ckd_stage2c_susie.md`, `docs/ckd_stage3a_kidney.md` describe the relevant analysis stages.
- **PARTIAL**: candidate-wise tissue evidence and cross-ancestry validation require dated output audit.
- **PARTIAL/UNRESOLVED**: known difficult loci (e.g. UMOD SuSiE failure) should be labeled, not made to PASS by brute-force retries.
- **BLOCKED BY DATA, NOT SOFTWARE**: KoGES controlled-access participant-level longitudinal validation, if claimed.

## Completion definition

A stage is **PASS** only when **(a)** original output and QC/provenance exist, **(b)** statistical assumptions/data matches are documented, **(c)** interpretation distinguishes association, shared-signal support, and causal mechanism, and **(d)** results feed a final manuscript table or figure. Code existence and successful CI are necessary engineering signals, not evidence of biological validity.

A stage is **PARTIAL** if only code, a snapshot, or incomplete evidence exists; **BLOCKED** when a required reference input is unavailable. Do not spend research time on data acquisition when the data cannot realistically be obtained; qualify the claim instead.

**Change policy:** new analysis may be added only to close a concrete gap in one of the four reference workflows or repair an identified validity problem. Any added stage must state the exact paper/method and the unsupported manuscript claim it repairs. Otherwise it is out of scope.

## Won lab alignment and final thesis direction

See [RESEARCH_DIRECTION_WON_HIGH_IMPACT_20261015.md](./RESEARCH_DIRECTION_WON_HIGH_IMPACT_20261015.md) for Hong-Hee Won publication workflows, IS/CKD-specific questions, statistical criteria, evidence tiers, figure plans, data-access blockers, and analysis stop gates. This document remains the minimal high-impact benchmark; the Won lab study structures are selected by research-question relevance rather than treated as extra mandatory analyses.
