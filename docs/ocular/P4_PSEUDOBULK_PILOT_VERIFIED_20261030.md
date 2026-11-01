# Ocular P4: Actual raw-count paired BeamA pseudobulk — verification and decision

Date: 2026-10-30 (KST). Disease context: glaucoma / in-vitro dexamethasone response of human trabecular-meshwork cells.

## Executed Pool C evidence
- **[Actions #41067820710](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41067820710)** PASS. [Aggregate evidence ZIP](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41067820710/artifacts/12983121752).
- Original files: GSM9373278 vehicle matrix 343,745,686 bytes, GSM9373287 DEX 375,519,107 bytes. All gzip CRC, MatrixMarket dimensions, NNZ and SHA-256 checks passed.
- Celltypist Chen2025 **fine BeamA**, confidence >= **0.75**, minimum **30 cells per donor per arm**; pool C final **6 paired donors**, not 7 (one donor does not pass stringent cell-state criterion).
- Vehicle: selected 3,727 matrix barcodes from 3,736 metadata candidate cells (99.759% recall), 84,295,022 integer UMI.
- DEX: 3,495 matched barcodes from 3,496 candidate cells (99.971% recall), 86,604,299 UMI.
- Raw UMI sums conserved exactly (170,899,321 in both arms combined). 29,597 genes detectable in either condition; 21,898 detected in >=3 donors in **each arm**.
- Per-donor median total UMI: vehicle 11,422,950.5, DEX 11,997,151. Per-donor median detected genes: vehicle 21,859, DEX 22,734.
- **No differential-expression p-values, false discovery rates, therapeutic response predictor, or causal inference were computed.**

## Full cohort design
- [P4 full cohort workflow](../../.github/workflows/ocular-p4-fullcohort.yml) streams all 18 public GSM raw 10x matrices **sequentially** on runner. The initial automated trial may be reconfigured to manual-only after evidence PASS. It does not upload or commit raw count matrices, donor IDs, or any donor-level gene matrices.
- Paired donor inclusion is defined **before examining gene expression**: Celltypist Chen2025 fine BeamA >=0.75; >=30 eligible SNG cells in each DEX/no-DEX arm; no post-hoc exclusion on DE results.
- 59 matched individuals are available in all cell metadata. P3 cell-state readiness suggests 57 BeamA-eligible donor pairs at cutoff .75. Final physical barcode joins may alter the count.
- Prespecified primary model after count QC: raw counts with donor blocking `~ donor + treatment`, with appropriate offset/library size and NB dispersion model (edgeR or DESeq2), gene filter, BH-FDR <0.05, global QQ/dispersion plots, donor leave-one-out, pool sensitivity and annotation sensitivity.
- JCT remains exploratory since only 8 pairs at high confidence .75. Schlemm endothelial DE is not inferentially powered.
- No distributional assumptions or EDA evidence support claims before the model is actually run.

## Novelty / paper precedence
1. GSE313736 deposit already reports 63 DEX-up and 17 DEX-down genes and sc-eQTL mapping. Do not call reproducing those DEGs novel ([dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736)).
2. Mehrotra et al., *IOVS* 2026, PMID 41636427, used paired bulk DESeq2 across 1 h, 6 h and 2 days at 100 nM DEX; found 857 TM and 2,086 SCE 2-day DEGs and shared ECM/adhesion pathways ([paper](https://pubmed.ncbi.nlm.nih.gov/41636427/)). Our GSE313736 is 500 nM 48h; only direction/pathway replication is reasonable without dose harmonization.
3. Lama et al., *Ophthalmology Science* 2026, PMID 42388255, clinical 6-month delta-IOP GWAS already reported UBE2E3, SSFA2 and GUCY2D colocalization signals. **No clinical steroid-response outcome is linked to the GSE313736 cultured-donor experiment** ([paper](https://pubmed.ncbi.nlm.nih.gov/42388255/)). Avoid pharmacogenomic prediction claims without obtaining clinical data.
4. A potential original result would be **donor-to-donor response heterogeneity with cell-state specificity, independent molecular replication and a mechanistic explanation**, not generic gene list overlap.

## Decision
- **GO** technical feasibility of pooled raw-UMI count construction in a paired six-donor pilot.
- **GO** whole-dataset QC scaling if all 18 matrix transfer and matched metadata joins succeed.
- **HOLD** gene-level DE statistics, novelty, independent replication and steroid-induced IOP risk claims pending direct analysis and validation.
