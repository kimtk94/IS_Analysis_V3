# P5 pre-registered statistical plan: steroid-response BeamA pseudobulk (2026-10-30)

## Primary analysis question
Does 48-hour dexamethasone treatment (nominal 500 nM as GEO record) alter human cultured TM **BeamA** transcript abundance reproducibly across biological donors? Is there substantial, reproducible inter-donor heterogeneity after QC and library normalization?

This is **not** an IOP, glaucoma incidence, drug-toxicity or clinical pharmacogenomics predictor. GSE313736 donors have no verified clinical dex-induced ΔIOP measurement.

## Eligibility frozen before inspecting gene-level DE p-values
- GSE313736 cell_metadata, assigned `SNG`; `celltypist_chen2025_fine == BeamA` and `celltypist_chen2025_fine_score >= 0.75`; >=30 high-confidence matched raw-matrix barcodes in **both** no-dex and DEX arms per biological donor.
- Source data span 9 donor pools, each represented by distinct vehicle and DEX 10x libraries. Do not treat individual cells or 10x droplets as independent replicates. Every donor has one pseudobulk profile per condition after summing **raw integers**.
- Public metadata identified 59 paired donors across all TM cells; P4 preflight predicted **57** BeamA-qualified donor pairs. A donor failing a threshold is excluded as a *pair*, not condition-specifically.
- Barcode recall >=90% per library; <95% is flagged for sensitivity. GSM9373291 preflight recall=92.7427% warrants leave-one-pool-G-out, reassess no-dex/dex cell capture and per-gene dropout.

## Differential-expression model (planned; not yet executed)
- Count input: 36,601 genes x 2N donor-condition pseudobulk integer counts, identical gene-ID order verified across 18 matrices; sum sparse counts **only** over eligible cells.
- In R edgeR: `DGEList(counts=raw_counts)`; `filterByExpr` with paired design, `calcNormFactors(method="TMM")`, `design <- model.matrix(~ donor + treatment)`, `stopifnot(qr(design)$rank == ncol(design))`, `estimateDisp(..., robust=TRUE)`, `glmQLFit(..., robust=TRUE)`, `glmQLFTest(coef="treatmentDEX")`.
- BH-FDR<0.05 is the primary genome-wide significance screen. Effect sizes (log2FC) and uncertainty should be shown for all tested genes. Pre-specify practical prioritization `abs(log2FC)>=0.25` alongside FDR, **not as an implicit statistical alpha replacement**.
- *Donor* captures donor-persistent differences; pool labels are nested within donor, so donor and pool fixed effects are linearly dependent. Batch confounding of treatment vs library handling cannot be eliminated by the paired donor term; assess consistency **across independent donor pools** and external experiment.
- Diagnostics: genes kept after filtering, library sizes, TMM factors, per-group zero counts, PCA/MDS, residual leverage, dispersion curves, sample outliers and signed log2FC direction.
- Sensitivity: leave-one-donor-out for prioritized genes; leave-one-pool-out for all 9 pools (not an independent validation cohort); alternative high-confidence annotation and different minimum cell criteria must be clearly exploratory.
- Heterogeneity: group-specific paired delta log-CPM; per-donor effect variance, direction concordance, treatment×covariate exploratory tests only when powered. Stable donor covariates (age, PRS) cannot be simultaneously fit as main effects with per-donor fixed effects due collinearity.

## Validation and novelty
- A 2026 TM/SCE dex time-course study (PMID 41636427; 100 nM DEX, 1h,6h,48h; paired bulk DESeq2) already discovered 857 TM and 2,086 SCE DEGs at 48h, including FKBP5 and FAM107A. This is a positive-control concordance source, *not* a new result.
- The original GSE313736 deposit already claims glucocorticoid DEGs and sc-eQTLs; basic DEG reproduction is low novelty.
- Clinical GWAS PMID 42388255 already reports UBE2E3, SSFA2, GUCY2D regulation near steroid-induced ΔIOP loci. Do not claim these targets are novel. Clinical association requires controlled dbGaP or independent approved data.
- PRJNA1243143 ENA provides 46 paired-end runs/BioSamples, but run count != donor count and automated keyword matching does not verify DEX/no-DEX/timepoint. External replication remains HOLD until samples map to actual treatment, donor and dose/time.

## GO / HOLD / STOP
- GO technical: >=50 paired donors passing raw count and barcode QC; 18/18 matrices, NNZ+CRC, gene order, total UMI conservation.
- GO inference: paired count-based NB QL model passes rank, dispersion, library, FDR diagnostics, effect-size stability, enough independent pool support.
- HOLD novelty: no contribution beyond the original sc-eQTL/DEX paper and published TM time-course; independent cohort not mapped.
- STOP clinical prediction: cultured-cell data alone cannot establish actual individual-level steroid IOP adverse effect.
- No release of donor ID, raw donor-level expression matrices, controlled genotypes or patient data in repository or public CI artifacts. Gene-level aggregate DE statistics and figures may be considered only after data-use and disclosure assessment.
