# MUSCLE MASTER — Precision Resistance Training

## Research question
Which genetic loci associated with inter-individual resistance-training response show convergent functional support in human skeletal muscle multi-omics?

## Scope
This project treats **resistance-training trainability** as the target phenotype. Sarcopenia, appendicular lean mass and grip strength may be used later as convergent population evidence, but they are not equivalent to training response.

## Stage plan
1. **Stage 0 — feasibility/audit:** literature, data availability, QC risks.
2. **Stage 1 — GWAS locus audit:** harmonize reported RT-response loci and recalculate significance tiers.
3. **Stage 2 — variant-to-gene:** skeletal-muscle eQTL/sQTL, LD proxies, regulatory annotation.
4. **Stage 3 — acute perturbation:** MoTrPAC resistance-exercise RNA/protein/phosphoprotein/chromatin evidence.
5. **Stage 4 — chronic exercise replication:** MetaMEx and OMAx.
6. **Stage 5 — hypertrophy responder validation:** GSE277819 low/medium/high responder transcriptomics.
7. **Stage 6 — convergent population evidence:** ALM/grip-strength GWAS and KoGES where appropriate.

## Stage 1 v1.1 design

### Primary RT discovery set
Two published RT cohorts are retained as **suggestive discovery evidence**:

- Yang et al. 2024, *Physiological Genomics*
  - RT n=179
  - phenotype: percent change in rectus femoris muscle thickness after 12-week RT
  - 11 reported lead SNPs
  - Methods threshold: P < 1e-5
- Gu et al. 2026, *Journal of Cachexia, Sarcopenia and Muscle*
  - RT n=187
  - phenotype: percent change in lean body mass after 12-week RT
  - 9 reported independent lead SNPs
  - Methods threshold: P < 1e-5
  - Abstract states P < 5e-8, which is inconsistent with the reported lead SNP P values and the Methods threshold

The primary audit therefore contains **20 RT loci**.

### HIIT comparator
The Yang et al. 2024 HIIT cohort contributes 8 lead SNPs as an **exercise-mode comparator only**. These loci are not counted as primary RT discovery loci.

## Statistical interpretation
The project uses two separate concepts:

- **Conventional GWAS significance:** P < 5e-8
- **Paper-specific exploratory threshold:** the threshold actually used in the source Methods section

Recalculated tiers:
- Tier A: P < 5e-8
- Tier B: 5e-8 <= P < 1e-5
- Tier C: P >= 1e-5

A locus can pass a paper's Methods threshold while remaining Tier B by conventional GWAS standards.

## Stage 1 v1.1 provenance fields
The audit records, separately:

- `methods_threshold`
- `abstract_claimed_threshold`
- `paper_calls_genomewide`
- `recalculated_genomewide`
- `methods_threshold_pass`
- `abstract_threshold_pass`
- `paper_internal_threshold_discordance`
- `nonstandard_genomewide_threshold`

This prevents terminology in the source papers from being conflated with conventional GWAS significance.

## Functional triage
Stage 1 functional evidence points are only a queue for Stage 2. They are not causal evidence.

Current hand-curated support includes:
- rs74038095 / SV2B: skeletal-muscle GTEx support reported in Gu et al. 2026
- rs4665972 / SNX17: skeletal-muscle GTEx support reported in Yang et al. 2024
- rs7924637 / CCDC15: skeletal-muscle GTEx support reported in Yang et al. 2024
- rs62149957 / LIMS1: skeletal-muscle GTEx support reported in Yang et al. 2024
- rs10212396 / ROBO2: CADD/RegulomeDB regulatory evidence reported in Gu et al. 2026

These annotations will be independently re-queried in Stage 2.

## Stage 1 outputs
- `MUSCLE_STAGE1_GWAS_AUDIT.tsv`
- `MUSCLE_STAGE1_SUMMARY.json`
- `MUSCLE_STAGE1_REPORT.md`
- `MUSCLE_STAGE1_HIIT_COMPARATOR_AUDIT.tsv`
- `MUSCLE_STAGE1_HIIT_COMPARATOR_SUMMARY.json`

## Sources
- Yang X et al. Genome-wide association study of exercise-induced skeletal muscle hypertrophy and the construction of predictive model. *Physiol Genomics*. 2024. DOI: 10.1152/physiolgenomics.00019.2024
- Gu Z et al. Genome-Wide Association Study of Lean Body Mass Response to Resistance Training in Young Asians. *J Cachexia Sarcopenia Muscle*. 2026. DOI: 10.1002/jcsm.70347
