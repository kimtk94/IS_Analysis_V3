# IS P3A — Original GTEx/ROSMAP sample-level QC summary counts, Release7 source crosswalk and GIGASTROKE GWAS Catalog v2 metadata

**2026-11-23 (KST). MASTER DEGREE / ischemic stroke only. Scientific status: HOLD.** All **36 authentic gene×tissue QTL groups (20 unique genes, ONE chr10 G0015 stroke locus)** were joined to original public metadata. No private or pseudonymous individual identifiers were uploaded; original raw sample IDs stayed within transient GitHub Actions memory during aggregation.

## Verified source files and published provenance

Original [eQTL Catalogue 2021 processed RNA count metadata Zenodo deposit](https://zenodo.org/records/4678936), publisher API [record 4678936](https://zenodo.org/api/records/4678936), original published record JSON **SHA256 818a5ebd7f448dd685f3f23d90d42f719a520a8fc27db233cb624bee21619682**, 32 declared public files. The two small original source metadata files were actually downloaded and verified against their **original publisher file MD5 and exact byte sizes**; no big count matrices were downloaded.

| Source | Exact original file, real published MD5 | Original bytes | Actually observed SHA256 | Total published metadata rows | Distinct pseudonymous source IDs in file (aggregate only) |
|---|---|---:|---|---:|---:|
| GTEx v8 | GTExV8.tsv; MD5 **0e29ac0664ad10403edf7a833d0a4503** | **1,769,117** | **9645c9754feeb51e600bcaafea47e99bc76c16d42a4bfb2710f3012ba91805c2** | **17,350** | **17,350** |
| ROSMAP | ROSMAP.tsv; MD5 **038f2fc144129d401144f3924109e1df** | **53,064** | **19a3718f2b6f3469e6218f41fa330c05e96fc363f86ada8800c041d7a825f69b** | **581** | **581** |

Source metadata have **13 named columns**, including pseudonymous ID, tissue qtl_group, sex, cell type, condition, **rna_qc_passed** and **genotype_qc_passed**, plus library and study annotations. Full source column names (but never any original row IDs) are included in the original CI JSON. The test workflow rejects malformed rows, duplicated column names, unexpected original source checksums and unknown Boolean flags in the three target groups.

### Independently counted RNA + genotype QC intersections in original public 2021 metadata

For each tissue the original pseudonymous samples were counted **once each**, QC flags classified separately from nominal QTL association data, and jointly passing samples compared with the independent Release7 eQTL Catalogue source tissue sample-N. The 2021 source has RNA QC PASS for all these target-group rows, while genotype QC FAIL explains GTEx source selection. **The 2021 source release and Release7 analysis-sample lineage are not thereby proven identical.**

| Study and public metadata qtl_group | Original 2021 tissue metadata rows / unique IDs | 2021 RNA PASS + Genotype FAIL | **2021 both QC PASS** | Independent Release7 nominal tissue N | Comparison |
|---|---:|---:|---:|---:|---|
| QTD000131 GTEx v8 artery_aorta | **432 / 432** | **45** | **387** | **387** | **Exact count agreement** |
| QTD000171 GTEx v8 brain_cortex | **255 / 255** | **50** | **205** | **205** | **Exact count agreement** |
| QTD000434 ROSMAP brain_naive (DLPFC) | **581 / 581** | **12** | **569** | **560** | **9 QC-pass samples not reconciled** |

**21/36 original gene×tissue groups** (10 GTEx Aorta and 11 GTEx Cortex) have original 2021 RNA+genotype QC pass counts numerically equal to the independently verified eQTL Catalogue Release7 tissue-wide sample N. **15/36 ROSMAP gene×tissue groups** retain the documented **569 vs 560 = 9 sample discrepancy**. The difference is not presumed to represent 9 excluded individuals for any specific reason: the source stage, genotype matching, molecular phenotype completeness, covariates and Release7 inclusion filters must be checked.

This audit changes the interpretation of earlier source AN/2 observations:

- Previous original eQTL variant-wise source AN/2 for 148,344 valid nonpal SNP×gene-assay associations: GTEx Aorta **366–387**, GTEx Cortex **186–205**, ROSMAP **560**. The 2021 source QC pass counts **387**, **205** and **569** are informative about tissue cohorts but not proof of each gene's actual regression sample inclusion; even though ROSMAP original variant AN/2=560, those 9 differing pseudonymous source samples cannot be assigned an unverified exclusion reason.
- The previous [primary complete source cis nominal vs permutation comparison](IS_P3A_36_PRIMARY_PERMUTED_CIS_NVARIANT_RECONCILIATION.md) verified exact source tested variant counts **36/36**. This is a variant denominator, not a sample N certificate.
- The [original eQTL Catalogue qtlmap specification](https://github.com/eQTL-Catalogue/qtlmap/blob/master/docs/usage.md) requires a **normalized phenotype input matrix**, phenotype metadata, genotype VCF and sample_meta, with **sample intersection** by QTL group. The 2021 deposited metadata and processed RNA counts are not automatically the **exact Release7 normalized qtlmap matrix** or a measured original gene-specific phenotype **sdY**.
- **Source raw or processed counts ≠ inverse-normal-transformed per-gene regression trait variance**. The source residual proxy from previous SE×sqrt(2f(1−f)AN/2) diagnostics is not measured sdY.

### Actual original 2022 GIGASTROKE study metadata, migrated to active GWAS Catalog API v2

The previous 2026 code path /gwas/rest/api/studies/ was retired and returned **HTTP 410**. Official [GWAS Catalog August 2026 migration documentation](https://www.ebi.ac.uk/gwas/docs/news/rest-api-v2-migration-guide/) requires **/gwas/rest/api/v2/studies/{accession_id}**. The new version was **actually queried successfully for both accession IDs**, their original study name and PMID, and initial/replication sample descriptions.

| Original GIGASTROKE AIS study | Original discovery cases | Original discovery controls | Official Catalog metadata replication stage cases | Replication controls |
|---|---:|---:|---:|---:|
| [EAS GCST90104545](https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90104545) | **19,032** | **237,242** | **8,304** | **174,855** |
| [EUR GCST90104540](https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90104540) | **62,100** | **1,234,808** | **59,911** | **665,176** |

These discovery counts independently agree with Mishra et al. 2022 primary original Supplementary Table 1 cohort-row arithmetic. The v2 public study record original SHA256 fingerprints were **329c7d49e6566236f31aef90525a7866450b1d0eebf5947c55c7f7ec22ceeb69 (EAS)** and **2f2a84b0b04eef5acf3302769abd79fdaf7a0fa92883bf7fe59837dd2c784630 (EUR)**; original PMID **36180795**, disease trait "Ischemic stroke." The Catalog reports a **replication stage sample description**, **not an automatically available replication GWAS summary-stat dataset or proof of participant-level independence**. GIGASTROKE original EAS discovery already contains **BIOBANK_JAPAN*** with 17,651 cases and 159,610 controls, so BBJ cannot be called independent replication of this outcome study without a BBJ-excluded discovery analysis.

### Exact GWAS beta effect-scale approval is still unresolved

The official v2 study records describe the disease and sample sizes, but do **not directly establish that beta in the exact original released AIS summary-stat files GCST90104545 and GCST90104540 is uniformly natural-log odds ratio**. The original 2022 publication describes **logistic regression OR Cox proportional hazards in some original longitudinal cohorts**, followed by fixed-effect inverse-variance weighted meta-analysis; not a direct proof that all final beta values have exactly the same log-OR interpretation.

GWAS-SSF [summary-statistics definitions](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format/) distinguish beta, odds_ratio and hazard_ratio. A submitted column labelled **beta** with standard_error is source evidence for effect-size/uncertainty availability, not proof of exact link function or odds-ratio conversion. The valid [coloc beta/varbeta type=cc pathway](https://chr1swallace.github.io/coloc/articles/a02_data.html) need not invent GWAS N_eff/s, but scientific approval of an appropriate outcome effect size, actual QTL phenotype sdY or N/MAF route and complete shared cis locus must precede calculation.

## Exact 36-group scientific gate

All real gene×tissue groups retain:
- 2021 source tissue sample and direct pseudonymous unique-ID count **only as aggregated numerical data**;
- group RNA+genotype dual-QC count and comparison with official Release7 nominal source N;
- GTEx 21/36 count-aligned or ROSMAP 15/36 count-mismatch flags;
- **Release7 gene-specific regression phenotype N not proven**, actual source processed sample intersection not independently reconstructed;
- **direct normalized QTL sdY not measured**;
- **original EAS/EUR GWAS file beta log-OR unit not confirmed**, even though official 2026 GWAS v2 metadata is now checked;
- **0/36 formal coloc.abf eligible/executed, 0/36 coloc.susie eligible/executed**; no posterior PP.H3/PP.H4, invented N_eff, fake sdY or causal gene attribution.

## Next source-driven scientific actions

**Priority 1:** Resolve the **nine ROSMAP QC-PASS vs Release7 N=560** samples by verifying qtl_group-specific normalized qcnorm counts, sample_meta/genotype intersection and original study-specific expression/covariates, without committing sample identifiers or guessing an exclusion criterion.

**Priority 2:** Obtain or reproduce actual **36-gene × three QTL study final normalized expression vectors and sdY** (or justify the exact genotype/phenotype model-based N+MAF approach), then review original GIGASTROKE outcome meta-analysis beta effect-scale from the original meta-analysis documentation. Only after both science approvals may the entire 36-assay batch be admitted to formal coloc.abf and subsequent sensitivity/LD-based follow-up. Investigate published GWAS Catalog **replication stage** source files/participants separately, without claiming independence from the stage summary alone.

Reproducibility: scripts/is/audit_p3a_36_original_zenodo_metadata_gwas_beta_model.py; tests/test_is_p3a_36_original_zenodo_metadata_gwas_beta_model.py; .github/workflows/is-p3a-36-2021-small-metadata-gwas-beta-audit.yml. Original aggregate evidence stored only as Actions artifact; [Issue #307](https://github.com/kimtk94/IS_Analysis_V3/issues/307) remains OPEN until scientific prerequisites satisfied. No CKM, Sales or other research components modified.
