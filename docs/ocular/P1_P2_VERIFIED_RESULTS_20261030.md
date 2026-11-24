# Ocular multiomics P1–P2 verified audit · 2026-10-30

## Final evidence and scope
- P1 [full paired count matrix CI](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41062752244) **success**, [artifact 12980741486](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41062752244/artifacts/12980741486).
- P2 [all-pool metadata CI](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41063135483) **success**, [artifact 12982425241](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41063135483/artifacts/12982425241).
- Clinical interpretation: **in-vitro DEX effect** in primary human trabecular-meshwork cells, **not** patient-observed steroid-induced ocular hypertension or risk prediction.
- Public raw matrix files were transiently downloaded into GitHub runners and **not committed**. Artifacts contain only SHA256 and aggregate summary, not donor names/genotypes/UMIs.

## Pool C actual 10x matrices
| Check | Vehicle: GSM9373278 | DEX: GSM9373287 |
|---|---:|---:|
| Features | 36,601 | 36,601 |
| Matrix columns/barcodes | 16,980 | 22,011 |
| Matrix Market declared and actual NNZ | 92,732,638 | 100,908,670 |
| Full decompression / gzip CRC | PASS | PASS |
| Matched metadata barcode fraction | 97.52% | 99.78% |
| Demux barcode coverage | 100% | **77.47%** |
| Biological donors in both arms | 7 | 7 |

DEX demux coverage is less than 90% but this does not invalidate the count matrix; inspect SNP read support and demux missingness before treatment comparisons. Exact donor assignment uses the separately published cell metadata.

## Nine-pool (18 libraries) metadata checks
All **18/18** libraries passed >=90% metadata-to-barcode matching; unique **59/59 donor IDs** pair between vehicle and dex, split: A=8, B=7, C=7, D=7, E=7, F=6, G=7, H=6, I=4. Final public metadata has **169,398 SNG-classified cells**. DEX pools B (demux 88.67%) and C (77.47%) warrant coverage QC. Original study abstract 74 donors / 316,799 singlets and GEO sample description 70 donor IDs vs final 59-cell-metadata donors remain an **unexplained attrition/mapping discrepancy** until paper QC is reconciled.

## External validation
Independent cohort PRJNA1243143: ENA **46 read runs**, **46 BioSamples**, all paired-end, with 46 FASTQ FTP fields. These are *runs/samples, not 46 independent donors*. No treatment pairing or timepoint mapping has yet been validated.

## GO/HOLD/STOP research gate
- **GO**: Pool C full raw count matrix/metadata integrity; metadata and barcode pairing for all 9 pools; implement analysis planning and download manifests.
- **HOLD**: Actual differential expression across all donors until the other 16 compressed matrices are verified and donor×condition×cell-state pseudobulk can be constructed. Independent replication design must be verified.
- **HOLD**: Clinical steroid ΔIOP genetic association (dbGaP phs004037) requires approved clinical phenotype/genotype. DME anti-VEGF patient×protein response matrix is not available in prior P0 attempts.
- **STOP**: Claiming a clinically validated pharmacogenomic predictor or novel simple DEX DEGs based only on these in vitro data.

## Novelty filter
Existing GSE313736 authors already reported 63 up-/17 downregulated dex response genes and performed sc-eQTL work. A 2026 multi-omics IOP/POAG chromatin paper reported 78 candidate genes; an independent 2026 TM/SCE time-course paper identified 857 TM and 2,086 SCE DEGs after 2 days in its QC-passing cell strains. Do not recycle basic DEG/pathway figures as novel. Advance only cell-state regulatory, donor-response heterogeneity, or clinically validated adverse-effect hypotheses with explicit independent confirmation and adequate power.

## Next
1. Batch transfer/verify 16 remaining `matrix.mtx.gz` files; record accession, file sizes, SHA-256 and integrity but do not git commit matrices.
2. Count donor×treatment×cell subtype, freeze exclusion criteria, and pseudobulk by donor (do not treat cells as independent biological replicates).
3. Verify PRJNA1243143 sample/timepoint/treatment mapping and independently reproduce directions, accounting for dose/time differences.
4. Build side-by-side paper figure matrix, specify what is already claimed vs what is still an original hypothesis.
