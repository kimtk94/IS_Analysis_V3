# IS P3A — official publisher Ischemic Stroke YAML and QTL allele-number sample audit

Date: 2026-11-22. Project MASTER DEGREE / ischemic stroke. **36 original gene×tissue QTL groups, 20 genes, ONE chr10 locus; coloc 0/36 HOLD.**

## Actual official-source execution

[GitHub Actions source audit #45033915115](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45033915115): **SUCCESS**. 15/15 focused negative tests, repository regression 335 tests (334 passed, one unrelated Rscript skipped).

Inputs were **original GIGASTROKE GWAS Catalog sidecar YAML files** for EAS AIS GCST90104545 and EUR AIS GCST90104540, and actual original eQTL Catalogue Release7 indexed GRCh38 nominal cis studies QTD000131/171/434 read over chr10:101837472–105353833. The exact 36 authentic assay identities and the 148,344 canonical nonpalindromic SNP×group record denominator match the previous original source audit; no 11GB download.

### Original GIGASTROKE GWAS published YAML evidence

| GWAS accession | Publisher trait_description | Build | Publisher original file MD5 | sample_size | case_count | control_count |
|---|---|---|---|---|---|---|
| EAS AIS GCST90104545 | Ischemic stroke | GRCh37 | 904e12e33ec3c2eeefe4dc93bec4995c | **MISSING** | **MISSING** | **MISSING** |
| EUR AIS GCST90104540 | Ischemic stroke | GRCh37 | 1bfe4ae8a5042fb24cb0a562b05b0d2f | **MISSING** | **MISSING** | **MISSING** |

Original publisher YAML entries also omit case_control_study and ancestry_category. Each source is identity verified against the whole-file MD5 previously computed. The original YAML files have date_metadata_last_modified=2025-01-05. Their own SHA-256 digests in this run are:
- EAS source YAML: 7383b859ce9bb4811c15bc634592ee1544788536797a5c663644e1aa648f1d3e
- EUR source YAML: 1c82cf26ec0554958ab4caebcce3083743077e4eb11e24830d5dbdc0305df247

The separate literature numbers EAS 19,032/237,242 cases/controls and EUR 62,100/1,234,808 may be used **descriptively**, but are **not substantiated by the accompanying original YAML**. Refer to the original GIGASTROKE 2022 Nature Supplementary Table 1 and source cohorts to independently approve the exact AIS discovery cases, controls and BBJ overlap; do not silently assign model case fraction or N_eff.

### Original QTL variant-specific allele count (an/2) vs nominal donor N

The original eQTL Catalogue pinned [Tabix column documentation](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/07918a267afdc1c84b5444bea9bd902e7201728f/tabix/Columns.md) states that allele number (an) is the number of observed alleles, usually two times the sample size for autosomal variants. ac is ALT count, ma_samples is minor-allele carrier donors, and maf is a separately reported minor frequency. This is **not** automatic exposure regression effective N approval.

| Original QTL source | Tissue | Source valid nonpal SNP×gene-assay rows | Valid even positive an | Source variant an/2 range | Catalogue sample N | At catalogue N | Odd/missing an | Invalid ac/an or ma_samples | MAF vs minor ac/an > 0.01 |
|---|---|---:|---:|---|---:|---:|---:|---:|---:|
| QTD000131 | GTEx aorta | 44,693 | 44,693 | 366–387 | 387 | 44,185 (98.86%) | 0 | 0 | 0 |
| QTD000171 | GTEx cortex | 45,570 | 45,570 | 186–205 | 205 | 45,132 (99.04%) | 0 | 0 | 0 |
| QTD000434 | ROSMAP DLPFC | 58,081 | 58,081 | 560 | 560 | 58,081 (100%) | 0 | 0 | 45 |
| **Total** | **36 assays** | **148,344** | **148,344** | — | — | **147,398** | **0** | **0** | **45** |

Repeated original rsID rows are deduplicated within each assay by original chromosome/position/ref/alt/molecular trait. The 45 ROSMAP MAF discrepancies are association rows rather than necessarily 45 distinct variants. Their original frequency representation may differ, and no alleles or frequencies are corrected without further source evidence. Full source count 148,344 matches the previous exact real QTL–GWAS matching audit.

The **GTEx observed per-SNP donor allele-N varies**. Applying a single study-level 387 or 205 to every SNP would conflate genotyped coverage with regression effective N. The technical per-SNP an/2 count is now observed; **per-QTL exposure model effective N remains unapproved**. For ROSMAP all 58,081 original canonical nonpal records have an/2=560, but that still does not independently certify model assumptions.

### Scientific approval unchanged

1. The original GWAS phenotype field is *Ischemic stroke*, build GRCh37, and original datafile is identity matched, but original GWAS YAML **does not provide cases, controls, sample size or effective N**. Primary original Supplementary Table 1 is the next evidence route, with BBJ overlap separated from independent replication.
2. The previous 3.52Mb QTL positional reference union and 36 real assays are present. Exact original per-assay phenotype coordinate/version and truly tested cis variants are not independently approved.
3. Source QTL an/2 count is **not** silently adopted as a regression effect model N. No unapproved sample harmonization, palindromic/indel variant imputation, signed LD or causal gene claims.
4. Formal coloc.abf 0/36 eligible, 0 executed. coloc.susie 0/36 eligible, 0 executed. No PP.H4 computed.

Scripts: scripts/is/audit_p3a_primary_gwas_yaml_qtl_an_36.py; tests/test_is_p3a_primary_gwas_yaml_qtl_an_36.py; workflow .github/workflows/is-p3a-gwas-primary-yaml-qtl-an.yml. Artifact includes public aggregate source metadata only, never original SNP beta/SE rows or private genotypes.
