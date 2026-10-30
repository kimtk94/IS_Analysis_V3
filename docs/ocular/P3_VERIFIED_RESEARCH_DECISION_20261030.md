# Ocular P3 verified integrative decision · 2026-10-30

## Executed evidence
- P0/P1/P2 public dataset audit and 59 paired donors: see `P1_P2_VERIFIED_RESULTS_20261030.md`.
- P3 full matrix CI [run 41064078334](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41064078334) passed 8/8 independent pool jobs and the aggregate check (artifact 12981837085). Together with previously verified pool C (run 41062752244), 18/18 10x matrices passed complete gzip CRC, MatrixMarket dimensions, and declared vs actual sparse NNZ; 16 newly checked matrices and 2 previously checked.
- Recovered filtered 10x matrices represent 350,171 barcodes (204,070 DEX + 146,101 vehicle), 1,751,655,185 nonzero records and 6.506 GB compressed bytes across all 18. **Do not equate filtered barcodes with the final 169,398 singlet metadata rows.** Mismatch/attrition and barcode selection to be documented.
- All 9 pools from P2 still have 59 unique paired donors by final metadata, not 59 clinical treatment response outcomes.
- P3 cell-state metadata [run 41064298979](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41064298979) passed the 59 donor and 169,398 cell metadata checks. Celltypist-fine score >= 0.5 and >=30 cells per arm: BeamA 57 matched donors, JCT 30, Ciliary_Muscle 8, Schlemm_Endothelium 0. With score >=0.75 and >=30 cells: BeamA 57, JCT 8. Alternative scpred TM_Fibro 41 donors at >=0.5, 21 at >=0.75; these are **competing annotations, not automatically interchangeable anatomical cell populations**.
- **Gene-order caveat**: 18 compressed features.tsv.gz files have differing compressed SHA-256, so *gene feature identity/order is not yet confirmed*. A separate decompressed gene-order CI is in `ocular-p3b-features.yml`; until its artifact passes, prohibit concatenation by row index.

## Proposed manuscript-worthy question and novelty protection
### Primary hypothesis (conditional)
**Paired-donor variability of high-confidence BeamA glucocorticoid response is associated with cell-state regulation and reproducible extracellular-matrix / aqueous-outflow biology.**

Avoid claiming this is a steroid-induced ocular hypertension risk model. Human eye cell cultures were DEX-exposed in vitro; no post-DEX ocular pressure phenotype for these 59 donors has been verified.

### Prior analysis requiring explicit coverage comparison
- GSE313736 authors already describe 74 donors, 316,799 singlets, 63 upregulated / 17 downregulated glucocorticoid response genes, and sc-eQTL mapping ([GEO GSE313736](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736)). These numbers differ from final downloadable post-QC metadata; not yet reconciled.
- Mehrotra et al., IOVS 2026 ([PMID 41636427](https://pubmed.ncbi.nlm.nih.gov/41636427/)): 100 nM DEX, 1h / 6h / 2d, 9 TM and 4 SCE QC-passing strains, 857 TM and 2,086 SCE DEGs at 2d. Separate PRJNA1243143 ENA registry has 46 runs/BioSamples, **not confirmed 46 donors or matched drug/control pairs**. GSE313736 uses 500 nM DEX for 48h per GEO descriptions: assess dose and time differences before replication claims.
- Singh et al., 2026 ([GSE301517](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301517)): dex-induced chromatin compartment / looping integration with IOP/POAG loci, 26 IOP +52 POAG candidate genes. Do not claim generic GWAS-overlap or DEX DEG enrichment as novel.
- New native human TM glaucoma nucleus atlas (2026) adds tissue context but published sequencing file accessibility remains to be independently verified.

### Statistical analysis GO/HOLD gates
1. **GO** data acquisition: all 18 matrices physically available through public GEO and verified in transient CI. Raw files are not persisted to Git or Drive.
2. **HOLD** integrated count matrix: wait for 18/18 *uncompressed gene IDs in identical order* or implement explicit feature harmonization.
3. **GO-design** paired pseudobulk: `~ donor + dex_treatment` by each prespecified cell state with raw **gene counts**; donor is biological replicate (59 possible), not 169,398 cells. Donor fixed effect subsumes pool in this nested design; do not add a collinear `pool` fixed effect. A model with gene expression outcome, treatment indicator and donor block must be tested for identifiability and low-count filtering.
4. **Primary cell-state** BeamA only under stringent label confidence (score >=0.75, >=30 cells per treatment and donor), re-run with alternative labels and score thresholds. JCT is exploratory: 30 donors at threshold 0.5 drops to 8 at 0.75. Exclude sparse Schlemm endothelial cell-type DE inference.
5. **HOLD statistical results**: count-level pseudobulk QC, library sizes, gene coverage, differential expression FDR, sample-wise outliers and independent concordance not yet executed.
6. **HOLD clinical drug safety hypothesis**: controlled-access clinical steroid-related delta IOP genetic phenotype from dbGaP phs004037 is separate and unverified. No patient-level risk predictions.
7. **HOLD external validation**: confirm PRJNA1243143 individual strain, treatment and timepoint layout. Cross-study 100 nM vs 500 nM dose and design mismatch mean not an exact replication.

## Resource governance
- Bulk P1/P3 transfer CI workflows are `workflow_dispatch` only to avoid re-downloading ~6.5GB on every PR edit. Small audit workflows are also made opt-in after their initial evidence passes. Trigger intentionally when actual evidence must be refreshed.
- Do not commit personal metadata, individual donor identifiers, genotype records, raw matrices, clinical information, or patient-level results. GitHub only includes scripts, MD and aggregate QC checks.
- The pull request is pending merge. No final paper-level GO claim yet.
