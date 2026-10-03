# MUSCLE MASTER — Precision Resistance Training / Skeletal Muscle Trainability

## Core question
Which genetic loci are associated with inter-individual hypertrophy response to resistance training, and which loci have convergent functional support in human skeletal-muscle regulatory and exercise-response omics?

## Evidence architecture
1. **Stage 0 — Feasibility/QC:** literature, data access, overlap risk, public-resource audit.
2. **Stage 1 — RT-GWAS audit:** independently reclassify reported loci and lock variant/build/allele metadata.
3. **Stage 2 — Variant-to-gene mapping:** skeletal-muscle eQTL/sQTL, LD, chromatin and coding/regulatory annotation.
4. **Stage 3 — Acute perturbation:** MoTrPAC resistance exercise RNA/protein/phosphoprotein/ATAC evidence; compare resistance vs endurance where available.
5. **Stage 4 — Chronic exercise replication:** MetaMEx and OMAx.
6. **Stage 5 — Hypertrophy responder validation:** GSE277819 low/medium/high responder transcriptomic evidence.
7. **Stage 6 — Convergent population evidence:** ALM/grip-strength GWAS and Korean cohort phenotypes where appropriate. These are not direct trainability replication.

## Statistical guardrails
- Do not call a locus genome-wide significant unless its numeric P is independently verified below 5e-8.
- Treat P < 1e-5 but >= 5e-8 as suggestive.
- Do not equate nearest gene with causal gene.
- Do not run MR/coloc without sufficiently complete summary statistics, alleles, genomic build, effect sizes and SE/variance information.
- Track ancestry and potential cohort overlap explicitly.

## Current gate
Stage 1 must produce a verified lead-locus table before Stage 2 begins.
