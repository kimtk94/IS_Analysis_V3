# IS EAS novelty and feasibility audit (2026-11-09)

**Purpose:** research design/evidence log only. Not submission text. The proposal template supplied for the class prohibits generative-AI authored submissions.

## Decision (conditional GO)

**Do not claim that an EAS-only or AF→cardioembolic stroke (CE) GWAS/eQTL analysis is novel.** Published work already combines Japanese AF genetics, EUR/EAS left-atrial eQTLs, atrial sn-multiome, AF genetics with CE outcomes, and EUR/EAS proteome-wide MR for ischemic stroke. A narrow, testable gap to audit is **whether ancestry-matched regulatory information changes locus-to-gene prioritization and the interpretation of transferable ischemic stroke (IS) susceptibility loci**. Novelty is *hypothetical* until a systematic locus/analysis comparison and feasible QTL access are established.

### GIGASTROKE ancestry-specific GWAS (known literature values; NOT CI-verified source files)

| Outcome | EAS accession | EAS cases | EAS controls | Intended use |
|---|---|---:|---:|---|
| Any stroke (AS) | GCST90104544 | 27,413 | 237,242 | sensitivity |
| Ischemic stroke (IS) | GCST90104545 | 19,032 | 237,242 | preferred primary EAS outcome |
| Cardioembolic stroke (CE) | GCST90104546 | 926 | 237,242 | exploratory only |
| Large artery stroke (LAS) | GCST90104547 | 1,735 | 237,242 | exploratory |
| Small-vessel stroke (SVS) | GCST90104548 | 5,532 | 237,242 | secondary prespecified |

**Do not confuse:** GCST90104544 is *any stroke*; GCST90104545 is *ischemic stroke*. EUR cardioembolic outcome GCST90104541 has 10,804 cases and 1,234,808 controls. These are counts, not an independent validation cohort. Some GIGASTROKE ancestry/subtype scans share people.

**Phenotype ascertainment update (2026-11-11):** [Stroke diagnosis and cohort audit](stroke_diagnosis_pheno_ancestry_protocol_2026-11-11.md) finds GIGASTROKE preprint appendix descriptions of physician/TOAST-adjudicated BBJ and Hisayama cohorts versus ICD/claims/registry-driven CKB and KCPS2. This supports ascertainment heterogeneity but does NOT yet verify exact study-to-accession representation in GCST90104544–48. External differences are not interpretable as ancestry-specific biology before this gate. Novelty competition also includes Surakka et al. Cell Genomics 2023 doi:10.1016/j.xgen.2023.100345 (ABO ancestry heterogeneity, ALDH2 ancestry–sex, and IS locus eQTL colocalization).

## Candidate comparisons, ranked

1. **Conditional GO: IS regulatory portability (EUR discovery, EAS outcome + ancestry-matched eQTL comparison).** Research question: At preselected IS-associated loci, does ancestry-matched blood cis-eQTL versus EUR-biased eQTL materially change colocalization support, lead variant interpretation, or gene prioritization? Primary estimands: proportion of analyzable loci with robust shared signals per ancestry; per-locus concordance after locus/LD coverage matching. Require EAS and EUR ancestry-specific LD for fine-mapping and colocalization, not cross-ancestry substitution. Do not interpret low EAS detection as true biological difference.
2. **Exploratory: CE/AF atrial mechanism.** EUR CE genetic signals interpreted with AF GWAS and atrial eQTL/snRNA, plus EAS CE *directional only*. EAS CE n=926 implies poor power; neither shared variant nor AF-independent causal pathway follows automatically. Do not call this EAS mechanistic replication.
3. **Alternative: IS versus SVS specificity in EAS.** SVS n=5,532 is less underpowered than CE n=926 but partial phenotype overlap makes contrast statistics and diagnostic definitions important. Direct subtype-case comparisons need compatible independent individual-level or robust summary methods; do not infer differential effects from significance in one and nonsignificance in another.

## Hypotheses and disconfirmation

H1: Ancestry-matched whole-blood eQTL (if full summary access is permitted) changes prioritized IS gene assignments relative to EUR-based eQTL at a measurable fraction of loci, after matched SNP/allele/QC/LD coverage. **Null:** no robust locus-level changes.
H2 (secondary): Putative shared causal variants at IS loci demonstrate directionally consistent effects in EUR and EAS where both are adequately powered. **Null / non-identifiability:** wide EAS CIs or insufficient colocalization information prevent assessment.
H3 (exploratory): Some IS loci have arterial/atrial/immune cell-type transcript expression consistent with a plausible mechanism. This is **biological annotation, not patient-level validation or causal proof**.

## Protocol guardrails

- Outcome is **genetic liability to IS** (case-control); not treatment response, recurrence or 6-month progression.
- Define GWAS discovery and EAS checks, study ancestry, sample overlap, covariates, TOAST phenotype, SNP build, allele handling, MAF and INFO thresholds *before* running tests.
- Select loci using discovery data only; do not train/reselect on same EAS evaluation set and call it independent replication.
- Harmonize to GRCh37/38, identify strand ambiguous SNPs and incompatible alleles, handle multi-allelic variants and population-specific allele frequencies. Confirm colocalization input region, shared SNP coverage and ancestry-appropriate LD.
- Prefer multi-signal fine-mapping + coloc-SuSiE when multiple signals are present; report PP(H3) versus PP(H4) and posterior sensitivity (prior settings); conditional colocalization is not proof of a gene's causal mediation.
- A gene's expression in atrial cells is *tissue context*, not an eQTL signal in those cells. Bulk blood and heart atrium are not biologically interchangeable.
- Correct multiple testing and report confidence intervals/effect estimates, not significant-vs-not-significant as evidence of heterogeneity. Assess heterogeneity formally, with attention to subtype GWAS power and sample overlap.
- Do not promise individual-level clinical AUC, prospective prediction, patient-matched omics, or wet-lab evidence with summary statistics only.
- For the 3–12 month proposal: scope one primary IS analysis and optional single tissue-context secondary analysis, *not* pQTL + MR + TWAS + single-cell + PRS all as equal aims.
- EAS MAEEA whole blood eQTL resource (2026, n=2,024) is published; actual **full summary-stat access, allele/LD fields and licensing remain unverified** and are a GO gate.

## Evidence and competition (publication, not proof of novel gap)

- Mishra et al., Nature 2022 GIGASTROKE, data availability GCST90104534–GCST90104563: https://doi.org/10.1038/s41586-022-05165-3
- EAS subtype counts and GCST90104544–48: https://www.mdpi.com/2308-3425/11/8/247
- 2023 Japanese AF multi-ancestry genetics and cardioembolic outcomes: https://doi.org/10.1038/s41588-022-01284-9
- 2024 cross-ancestry atrial eQTL/sn-multiome, 14 genes at nine loci: https://doi.org/10.1016/j.isci.2024.110660
- 2025 AF pathway association with cardioembolic stroke: https://pubmed.ncbi.nlm.nih.gov/40526779/
- 2025 EUR/EAS ischemic stroke proteome MR: https://pubmed.ncbi.nlm.nih.gov/40304040/
- 2025 AF cross-population GWAS and proteomics: https://doi.org/10.1038/s41467-025-61720-2
- 2026 MAEEA EAS whole blood eQTL (n=2,024): https://doi.org/10.1093/nsr/nwag566
- 2024 stroke subtype cell-specific eQTL coloc: https://doi.org/10.1016/j.compbiolchem.2024.108086
- GWAS Catalog REST v2 documentation (v1 retired Aug 2026): https://www.ebi.ac.uk/gwas/docs/programmatic-access/rest-api/

## Decision gates

- **G0 Source/phenotype:** five EAS IDs and GWAS Catalog trait metadata have passed CI; now retrieve full GWAS bytes/headers and **map per-cohort clinical/TOAST vs ICD/claims definitions from final supplementary materials**. Preprint examples are partial evidence, not phenotype-level validation.
- **G1 Regulatory source:** full eQTL summary-stat access confirmed for EAS and comparator, including SNP/build/effect/SE/MAF/N, plus matched LD.
- **G2 Statistical feasibility:** estimate local power / credible set overlap and define a small prespecified discovery-locus set; cannot label low powered EAS CE null as biological difference.
- **G3 Novelty:** publishable locus-wise novelty matrix proving candidate interpretation not already reported; if not, STOP or change focus.
- **G4 Execution:** CI analytic tests with real loci (no fixture data passed off as biological results), independent evidence and complete audit trail.

Current grade: **CONDITIONAL_GO / DATA_NOT_AUDITED**. CI metadata ping does not count as biological analysis or verify summary-stat bytes.
