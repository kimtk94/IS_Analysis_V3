# IS P3A — Official ROSMAP R6/R7/R8 release lineage and qtlmap input inventory

**2026-11-23 KST — ischemic stroke only.** 36 original gene-by-tissue assays, 20 genes, one chr10 G0015 stroke locus. Scientific status HOLD; formal coloc.abf and coloc.susie 0/36.

## Research question

Does 2021 ROSMAP metadata RNA/genotype QC PASS N569 versus eQTL Catalogue QTD000434 nominal N560 arise *newly* in Release7? Does the official public QTL FTP path table directly publish the Release7 normalized qtlmap expression, sample metadata or genotype inputs?

## Independent source evidence

| Exact authoritative dataset | ROSMAP QTD000434 brain_naive N | Meaning |
| --- | ---: | --- |
| Release6 official dataset_metadata_r6.tsv | 560 | Published QTL sample N already 560 in R6 |
| Release7 official dataset_metadata_r7.tsv | 560 | No *numerical* change from R6 |
| Release8 official dataset_metadata_r8.tsv | 560 | Same nominal N persisted |
| Official FTP tabix_ftp_paths.tsv | 560 | Association and SuSiE outputs listed |
| 2021 Zenodo ROSMAP.tsv joint QC PASS | 569 | Earlier independently verified QC metadata, not R7 final association sample identities |

The 569 - 560 = 9 discrepancy **cannot be attributed specifically to the Release 7 transition**. In the earlier PR #320, the independent historical official QC report also declared 581 before and 576 after QC, but these reported sample counts are from multiple distinct source stages. The differences are not certified sequential exclusions or exact sample-level matches.

## What the original public FTP inventory covers

[Official data access page](https://www.ebi.ac.uk/eqtl/Data_access/) links the [Release7 metadata](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/data_tables/dataset_metadata_r7.tsv) and the [official QTL output manifest](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/tabix_ftp_paths.tsv). For QTD000434 the exact output fields are ftp_path (full-cis association summary statistics), ftp_cs_path (SuSiE credible sets) and ftp_lbf_path (SuSiE log Bayes factors). The manifest does not offer input fields named count_matrix, sample_meta, pheno_meta, vcf or covariates.

The [official qtlmap input format](https://github.com/eQTL-Catalogue/qtlmap/blob/master/docs/usage.md) independently requires the *normalized* expression count_matrix, phenotype metadata, sample metadata (sample ID/genotype ID/qtl group), matching genotype VCF, and the intersection of the two sample sets. Published association results alone do not measure gene-specific phenotype sdY.

**Critical negative-evidence boundary:** This public table is an inventory of association/fine-mapping *outputs*, not an exhaustive inventory of source studyFiles or all available preprocessed sample-level inputs. Its lack of input columns **does not prove** the exact qtlmap files are unavailable elsewhere. The [official ROSMAP study page](https://www.ebi.ac.uk/eqtl/Studies/) points to original raw study data in Synapse syn3219045; this is not the exact reconstructed Release7 normalized input or a license to use donor-level data.

The [official release notes](https://www.ebi.ac.uk/eqtl/Release_notes/) describe R6 (April 2023), R7 (June 2024) and R8 pre-release (January 2026). Equal N across R6/R7/R8 does not prove sample IDs, genotypes, covariates or phenotypes were identical.

## Reproducible pipeline

- Source script: scripts/is/audit_p3a_rosmap_release_input_inventory.py
- Synthetic negative tests: tests/test_is_p3a_rosmap_release_input_inventory.py
- GitHub Actions: .github/workflows/is-p3a-rosmap-r6-r8-input-inventory.yml
- Live original HTTPS reads: R6, R7, R8 metadata, official FTP path manifest, official qtlmap input specification. Source SHA-256 and byte counts preserved in aggregate JSON.
- 36 distinct gene-by-tissue groups retained, 15 ROSMAP, 20 unique genes.
- No original sample IDs, raw genotype, donor-level data or expression values stored.
- All 36 groups scientifically HOLD. No exact Release7 sdY measured, no gene-specific regression N approved, no GIGASTROKE final log(OR) units confirmed, no official coloc posterior computed.

## Next source evidence gate

Find the authentic Release7 ROSMAP studyFile with exact normalized expression matrix, phenotype metadata, sample_meta/genotype intersection and covariate model. Or, where scientifically justified, establish and document an equivalent quantitative-QTL N + MAF pathway. Separately review original GIGASTROKE AIS beta units. Do not promote technical CI success to scientific Coloc eligibility.
