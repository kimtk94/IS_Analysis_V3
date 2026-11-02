# Ocular P1 — count matrix integration and independent replication

Objective: establish whether **GSE313736** paired dexamethasone / vehicle human trabecular meshwork scRNA-seq can support biological-donor-level pseudobulk; do **not** infer the individual steroid-induced IOP phenotype from cultured cells.

## Files verified by P0
- GSE313736 cell metadata: 169,398 rows, 59 donors observed in both dex and no-dex conditions; original paper/record 74 donors and 316,799 cells, GEO descriptions 70 individuals. QC attrition must be explained before manuscript.
- Metadata includes `library_id`, `barcode`, `donor_id_demuxlet`, `DEX Treatment`. Representative pool C vehicle GSM9373278, DEX GSM9373287 has barcodes and demux.
- P0 cannot conclude raw count matrix is readable because count matrices not transferred.

## P1 automated test
Run `.github/workflows/ocular-p1-matrix-integrity.yml`, source `scripts/ocular/ocular_p1_integrity.py`.
- Download complete **public** GSM9373278 and GSM9373287 pool C `matrix.mtx.gz` plus barcodes/features/demux files, streamed to transient GitHub runner disk (620 MiB per matrix cap). No FASTQ, genotypes or clinical microdata.
- Validate gzip CRC by streaming every decompressed coordinate row; compare claimed nonzero count and 10x matrix dimensions to barcodes/features.
- Identify metadata library via the strongest within-treatment barcode match. Require >90% metadata and demux joins. Count biological donors and donor overlap. Do not output donor identifiers or expression values.
- Probe ENA independent registry for PRJNA1243143; presence of SRA runs does not establish paired treatment metadata or any successful replication.
- The workflow's `success` means the script ran, **not** that full matrix verification succeeded. Only JSON `status=GO_POOL_C_FULL_COUNT_MATRIX_INTEGRITY` is a positive integrity gate.

## GO/HOLD/STOP thresholds
- **GO matrix pilot**: both matrices byte-downloaded; CRC/NNZ and dimensions validated; >=90% metadata/demux barcode join; >=1 shared biological donor in pool C.
- **HOLD full analysis**: the above is insufficient for generalization; run all 9 donor pools with UMI/celltype QC, >=10 evaluable donors with adequate per-condition cells for primary cell types, donor-level pseudobulk, formal QC and independent cohort.
- **HOLD clinical risk**: dbGaP phs004037 controlled genotype+steroid IOP phenotype, not provided by GSE313736.
- **STOP naive novelty**: 2026 time-course and multiome papers already cover simple dex DEGs & gene enrichment. Need novelty via testable cell-specific regulatory mechanisms or external replication.
- **HOLD DME**: anti-VEGF patient×protein matrix unavailable to this audit.
