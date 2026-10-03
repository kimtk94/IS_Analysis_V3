# MUSCLE MASTER — Precision Resistance Training

## Research question
Which genetic loci associated with inter-individual resistance-training response show convergent functional support in human skeletal muscle multi-omics?

## Scope
This project treats resistance-training trainability as the target phenotype. Sarcopenia, appendicular lean mass and grip strength may be used later as convergent population evidence, but they are not equivalent to training response.

## Stage plan
1. **Stage 0 — feasibility/audit:** literature, data availability, QC risks.
2. **Stage 1 — GWAS locus audit:** harmonize reported RT-response loci and recalculate significance tiers.
3. **Stage 2 — variant-to-gene:** skeletal-muscle eQTL/sQTL, LD proxies, regulatory annotation.
4. **Stage 3 — acute perturbation:** MoTrPAC resistance-exercise RNA/protein/phosphoprotein/chromatin evidence.
5. **Stage 4 — chronic exercise replication:** MetaMEx and OMAx.
6. **Stage 5 — hypertrophy responder validation:** GSE277819 low/medium/high responder transcriptomics.
7. **Stage 6 — convergent population evidence:** ALM/grip-strength GWAS and KoGES where appropriate.

## Stage 1 QC rule
Published labels are not accepted at face value. Loci are reclassified from the numeric P value:
- Tier A: P < 5e-8
- Tier B: 5e-8 <= P < 1e-5
- Tier C: P >= 1e-5

The 2026 JCSM paper reports nine loci as genome-wide significant in the abstract, while Table 1 gives P values from 9.28e-7 to 8.71e-6. The seed audit therefore treats all nine as **suggestive Tier B** until full summary statistics or an author correction supports a different interpretation.

## 2026 seed study
- Gu et al. 2026, *Journal of Cachexia, Sarcopenia and Muscle*
- DOI: 10.1002/jcsm.70347
- PMC: PMC13371584
- n=187 Chinese Han young adults
- phenotype: percentage change in lean body mass after 12-week resistance training
- use: discovery/suggestive loci only; not independent replication

## Stage 1 outputs
`MUSCLE_STAGE1_GWAS_AUDIT.tsv`: row-level recalculated QC and triage.

`MUSCLE_STAGE1_SUMMARY.json`: machine-readable counts and discordance list.

`MUSCLE_STAGE1_REPORT.md`: human-readable audit report.

## Interpretation guardrails
Functional-priority points are used only to decide what to validate next. They do not establish causality. Population muscle-mass or grip-strength associations are convergent evidence rather than replication of exercise trainability.
