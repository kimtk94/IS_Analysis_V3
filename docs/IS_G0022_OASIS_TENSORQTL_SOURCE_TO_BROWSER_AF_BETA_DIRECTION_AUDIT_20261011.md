# IS G0022 — OASIS original-to-browser tensorQTL AF and effect-size lineage

**Executed:** 2026-10-11 KST
**Principal result:** Source-level Bokeh "maf" = original tensorQTL "af" to display precision; very strong EAS ALT AF cross-study orientation agreement; valid biological mediator/colocalization not established.

## 1. Original 39GB DDBJ dataset: real first 256KB only

Accessed the official public OASIS E-GEAD-1054 processed.zip through a strict HTTP 206 first-262144-byte request (not a full 39GB download). Decompressed just ZIP raw-DEFLATE → first TAR file B_Activated_PC15_MAF0.05.cis_nominal.txt.gz → first gzip text prefix, yielding exactly **1,010 complete LINC01409 GRCh38 cis-SNP association records** from the first 110000 output bytes. The matching official Japan Omics gene page OASIS-B_Activated-L2 for ENSG00000237491 contains 2,492 Bokeh SNPs.

**Direct original-to-browser audit:**
- **1,010/1,010 source SNP IDs occur in the browser data.**
- Raw pval_nominal ↔ Bokeh pval_nominal; raw slope ↔ Bokeh effect_size; raw slope_se ↔ Bokeh effect_size_SE; raw af ↔ Bokeh field named maf.
- All four numeric columns agree with a maximum absolute discrepancy 0.0004998978, consistent with public browser three-significant-digit rounding.
- Raw AF >0.5: **311 original records**. Browser field maf >0.5 for exactly those 311.
- This is direct proof the Bokeh **maf name is a display label for the original af, not mathematical minor allele frequency** in the one directly traced gene-cell raw source.

The original source is preserved as a bounded 262KB byte snapshot with SHA256, and the browser HTML is also independently snapshotted with SHA256. 1,010 association-by-SNP evidence rows are preserved in TSV. Do not extrapolate that the *entire* 40-cell-type original dataset has identical sampled-SNP coverage from one partial TAR entry.

## 2. OASIS authors' method scripts (higher-grade provenance)

The original OASIS authors' public repository exposes the exact source mapping script:

- cis-eQTL mapping: https://github.com/REdahiro/OASIS_project/blob/main/scripts/cis-eQTL_mapping/tensorQTL_cis_nominal.py — GitHub blob SHA 5158f6c4429b1e5c8f43400664ce02b4afb94a7b.
- It calls tensorQTL cis.map_nominal with PLINK input genotype, maf_threshold=0.05, and exports the original columns pval_nominal, slope, slope_se, af, ma_samples, ma_count in cell-specific compressed text.
- OASIS coloc example: https://github.com/REdahiro/OASIS_project/blob/main/scripts/colocalization/coloc.R — blob SHA 26159f9681df73fbf72547893b7183da3a1297f5. This author script explicitly maps beta=slope and MAF=af. Here MAF is a misnamed field: converting AF into actual MAF would require min(af,1-af). However many biallelic variance formulas depend on p(1-p), which is symmetric; the naming alone does NOT prove a colocalization implementation error. Full input quality/prior diagnostics would still be needed.

The public tensorQTL source defines af as genotype dosage allele frequency, and the PlinkReader implementation contains a deliberate 2 - genotype transformation with an explanation about PLINK REF coding. This strongly supports the authors' *intended* dosage-ALT effect convention:
- https://github.com/broadinstitute/tensorqtl/blob/master/docs/outputs.md
- https://github.com/broadinstitute/tensorqtl/blob/master/tensorqtl/core.py
- https://github.com/broadinstitute/tensorqtl/blob/master/tensorqtl/genotypeio.py

**Remaining build/allele qualification:** the authors' actual private PLINK BIM/BED allelic order cannot yet be independently inspected. Therefore a browser G>A variant label plus matching AF/source code establishes strong orientation consistency, not a universal variant-by-variant certificate that the original genotype count is ALT for every SNP.

## 3. Independent internal orientation sanity check at the actual G0022 locus

In previous independent read-only work, 5,789 source-QC GIGASTROKE East Asian acute ischemic stroke (AIS) GCST90104545 GRCh37 variants were lifted via UCSC official chr12 chain to GRCh38, and 8,921 exact coordinate+REF/ALT OASIS Bokeh gene-cell-SNV association records were joined. There are **2,866 unique AIS SNPs**, not 8,921 independent SNPs.

For each unique matched variant, the GIGASTROKE ALT EAF was compared with the median observed OASIS source-mapped AF across available gene-cell records. This is a *descriptive cross-population frequency audit*, **not independent-donor replication**.

| Metric, 2,866 unique matched variant identities | Result |
|---|---:|
| Pearson r, GWAS ALT EAF vs OASIS source-like AF | **0.9958820914** |
| Median absolute ALT frequency difference | **0.0173** |
| Mean absolute ALT frequency difference | **0.0237549** |
| Mean absolute difference to complemented (1 − GWAS EAF) | **0.5480730** |
| Same-allele frequency closer than flipped frequency | **2,846 / 2,866** |
| Flipped allele frequency closer | **19 / 2,866** |
| Ties (e.g. ~0.5 EAF) | **1** |
| ALT AF difference >0.1 | **0** |
| 2,866-median OASIS AF observations >0.5 | **1,180** |

In all five populated gene-cell strata separately the GWAS-vs-OASIS allele-frequency r is 0.9953–0.9964. These strong data-level results reinforce that, in this regional source, the "maf" Bokeh column is ALT-consistent AF and not necessarily <0.5. Minor exceptions are very plausibly due to stochastic between-cohort AF differences near 0.5 and do NOT justify flipping those specific alleles without source genotype reference confirmation.

**Science gate updated from UNKNOWN_FIELD_NAME to SOURCE_VERIFIED_AF_LINEAGE; genome-wide/in-study effect-allele reference genotype verification is still NOT COMPLETE.**

## 4. Scientific impact and remaining barriers

This audit advances frequency and effect column semantic confidence, but it DOES NOT establish disease–QTL causal colocalization:

1. The direct 1,010-variant raw-source check is on **LINC01409/activated B L2 (chr1)**, not original ALDH2/BRAP/RPH3A source TAR rows. Generalizing source field labels across same browser rendering pipeline is a strong inference, not a directly repeated original file audit per cell type.
2. The other 3 G0022 AIS 95%-CS SNPs (rs671, rs11066132, rs77768175) are **not displayed** in the five relevant Bokeh source gene pages despite being inside gene cis spans; neither donor genotyping absence nor a biological null is established.
3. The full-tested cis SNP universe, origin cell-type-specific participant N/covariates/INFO, and in-cohort signed QTL LD remain unavailable. Three-gene competitive causal comparisons and SuSiE-coloc stay CLOSED.
4. OASIS Japanese COVID/control PBMC eQTL signals and stroke GWAS sample ancestry/covariate composition are not independent disease-mechanism confirmation. Keep the 2,225-gene IS candidate scope open.

Use author methods and these **source-backed column correspondences** in a future properly complete GWAS–QTL SNP intersect, but do not prematurely report colocalization PP.H4, MR or clinical therapeutic effects.

## 5. Code, tests, artifacts

Git additions:
- scripts/is/audit_g0022_oasis_tensorqtl_source_semantics.py
- scripts/is/audit_g0022_oasis_allele_frequency_concordance.py
- tests/test_is_g0022_oasis_tensorqtl_source_semantics.py
- tests/test_is_g0022_oasis_allele_frequency_concordance.py

Read-only standalone artifacts:
- /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/oasis_raw_browser_semantics_v1/
- /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/oasis_af_beta_semantics_v1/

The source-to-browser proof report and G0022 frequency-concordance report retain SHA256 provenance, 1,010 direct row-level raw-vs-Bokeh comparisons and 2,866 unique variant AF direction rows.

**No large DDBJ source bulk download, source original modifications, reference-genotype access to restricted PLINK, or production-service writes performed.**
