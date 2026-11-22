# IS P3A — original GIGASTROKE 2022 AIS discovery cohorts and BBJ independence gate

**Date:** 2026-11-22. **Scope:** MASTER DEGREE ischemic stroke P3A, not Sales/CRM, not CKM code modifications. **36 genuine QTL gene×tissue assays from 20 unique genes, ONE chr10 G0015 locus**.

## Original publisher evidence and real GitHub Actions

The original GIGASTROKE Nature 2022 paper is [Mishra et al., Stroke genetics informs drug discovery and risk prediction across ancestries](https://doi.org/10.1038/s41586-022-05165-3). The original publicly released Supplementary Table 1, file **41586_2022_5165_MOESM4_ESM.xlsx**, was downloaded directly from the Springer static publisher site on the GitHub Actions runner. **756,910 bytes**, original source workbook SHA-256 **7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4**. The job parsed the original OOXML sheet **Suppl.Table1** without a third-party spreadsheet transformation or changes to other research code.

[Actual source verification CI #45035576990](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45035576990) **SUCCESS**: 15/15 focused controls, complete research code regression **350 tests executed, 349 PASS and one skipped (Rscript unavailable)**.

**Critical source-boundary rule:** The workbook includes discovery cohort summaries, followed by secondary sections where EUR appears again. The source rows **5–49 EUR** and **50–65 EAS** form the first discovery cohort table, ending at **Total ALL discovery (row 86)**. Repeated ancestry headers in later evaluation/PGS sections **are not included**. Any mapping that combines those sections falsely alters discovery sample size. The case column is **K, AIS**, not J (any stroke); discovery controls are **P**, not post-discovery evaluation samples.

### Independent original-row sum reconciliation

| Publisher-original source and coordinates | AIS cases in K | Discovery controls in P | Nominal total | Actual original published K/P sum equals recomputed cohort rows |
|---|---:|---:|---:|---|
| EAS discovery, source rows 50–65 | **19,032** | **237,242** | **256,274** | **YES** |
| EUR discovery, source rows 5–49 | **62,100** | **1,234,808** | **1,296,908** | **YES** |

The EAS source has 15 cohort rows in the first EAS block, with **one incident-only entry at row 60 excluded from the combined AIS discovery summary**; 14 contribute the published discovery total. EUR has **44 published study rows**. For cases and controls the original published total row equals the sum of included study rows. Some studies have missing individual controls or use shared discovery control definitions; **cohort-row arithmetic is not verification of individual participant-level nonoverlap** or of per-SNP meta-analysis effective N.

### Exact BBJ role in EAS discovery: independent replication is invalid

Original source EAS AIS case-contributing studies are:

| EAS discovery contributing cohort | Original AIS cases K | Discovery controls P |
|---|---:|---:|
| **BIOBANK_JAPAN*** | **17,651** | **159,610** |
| KCPS2_ASN** | 46 | 12,366 |
| HISAYAMA_FSR | 1,113 | 901 |
| INTERSTROKE_ASN | 222 | 698 |
| **Total four AIS-case contributors** | **19,032** | **173,575** |

The remaining EAS controls are supplied by **10 China Kadoorie Biobank (CKBB_ASN) regional rows** without AIS case counts, adding **63,667 controls**, so the original EAS discovery control total is **237,242**.

**BioBank Japan is a primary discovery contributor of GIGASTROKE EAS AIS**, not an external independent replication cohort. BBJ supplies **17,651 / 19,032 (~92.7%) of EAS discovery AIS cases** and **159,610 / 237,242 (~67.3%) of controls**. Consequently, testing BioBank Japan AIS independently against this GIGASTROKE EAS GWAS cannot establish independent external replication; use a BBJ-excluded outcome GWAS or an unrelated cohort after checking shared subjects and genotype/phenotype construction.

**Do not confuse** the source cohort named BIOBANK_JAPAN* with an independently retrieved BBJ GWAS of another phenotype or release. Even a cohort name/size difference does not certify absence of participant overlap.

### EUR original cohort composition

Original Table 1 comprises **44 EUR study rows** (including UK Biobank, METASTROKE, FinnGen hospital/population contributions and other cohorts) and published case/control totals above. The machine report retains each source study label and original K/P values with their source row numbers. EUR is **not** BBJ-based in this GIGASTROKE discovery table, but this does not by itself certify independence of other specific validation GWAS.

### Connection to the 36 real QTL source groups

The pipeline requires actual immutable original 36 QTL assay identifiers and the [independent QTL variant allele-number and original GWAS YAML evidence](IS_P3A_PRIMARY_GWAS_QTL_VARIANT_AN_AUDIT_20261122.md). Each one of 36 original groups is annotated with EAS discovered cohort and BBJ dependency, and explicit science blockers. It does not infer SNP-level N from the aggregate cohort sum.

- **Publisher original YAML files GCST90104545 EAS AIS and GCST90104540 EUR AIS** identity/trait/build/MD5 verified, but the YAML source files have NO case count, controls or sample size fields. The original supplemental source now independently confirms published **discovery** cohort totals at the cohort level; **exact publisher summary-stat file N applicability and SNP-specific effective meta-analysis N are still not verified**.
- Original eQTL Catalogue QTL genotyped allele count an/2 was independently checked for all **148,344 original nonpalindromic SNP×gene/tissue source associations**, varying GTEx aorta **366–387**, GTEx cortex **186–205**, ROSMAP DLPFC **560**. This is **genotyped donor evidence, not automatically the exposure regression N used for ABF**.
- The previous [SNP-level GIGASTROKE QTL harmonization](IS_P3A_36_SNP_ORIENTATION_GWAS_AUDIT_20261122.md) confirms EAS **89,597** and EUR **128,282** aligned original allele-oriented SNP×assay records, but all 36 share one locus. The remaining allele-bias exceptions, palindromes and INDELs are held separately.
- The full gene-specific nominal tested cis universe, exact molecular phenotype-version matching and independently justified per-group QTL regression N remain needed. Nonoverlap/LD checks for advanced conditional models remain distinct.

## Approval boundary

**Primary discovery aggregate EAS/EUR AIS cohort counts: technically VERIFIED from original publisher spreadsheet**. **BBJ included in EAS discovery: technically VERIFIED and therefore independent BBJ validation PROHIBITED**. This can support source description of the original discovery design, but **does not automatically authorize coloc ABF input N or case fraction**, nor any independent validation claim.

**Formal coloc.abf eligible 0/36, executed 0. Coloc.susie eligible 0/36, executed 0. No PP.H0–H4 computed, no causal stroke gene conclusion or independent external replication asserted.**

This is an original-source manuscript/QC-ready cohort table, not an authorization to start hypothesis testing. Open scientific Issue #307 remains active.

Source audit: scripts/is/audit_p3a_gigastroke_original_supp_cohort_36.py. Synthetic tests: tests/test_is_p3a_gigastroke_original_supp_cohort_36.py. Workflow: .github/workflows/is-p3a-original-2022-cohort-table1.yml. Original XLSX source SHA-256 and only cohort-level public source values and 36 identity/HOLD rows are in the CI artifact, not original individual SNP statistics, private genotypes or patient-level records.
