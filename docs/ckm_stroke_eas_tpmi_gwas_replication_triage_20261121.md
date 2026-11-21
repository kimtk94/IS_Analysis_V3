# East Asian AIS GWAS replication source audit

Date: 2026-11-21. Technical source QC PASS. Scientific status HOLD. Draft PR #284 remains open and unmerged.

## Verified public data sources

### Taiwan Precision Medicine Initiative (TPMI) PheWeb, GRCh38

The official site https://pheweb.ibms.sinica.edu.tw/about states it provides GWAS summary statistics and PGS model resources from the TPMI Han Chinese program with GRCh38 variant coordinates. Official phenotype pages and case/control totals were independently retrieved in [GitHub Actions #44871010169](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44871010169), source artifact [#14694240592](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44871010169/artifacts/14694240592); 10/10 offline fail-closed tests passed.

| PheWeb phenotype | Observed source cases | Controls | Clinical interpretation |
| --- | ---: | ---: | --- |
| 433.21 Cerebral artery occlusion, with cerebral infarction | 9,249 | 304,660 | Highest-priority AIS-like infarction phenotype; not necessarily TOAST-imaging-confirmed |
| 433.2 Occlusion of cerebral arteries | 12,463 | 304,660 | Broader occlusion, can include no definite infarction |
| 433.3 Cerebral ischemia | 12,468 | 304,660 | Broad cerebral ischemia, not equivalent to acute ischemic stroke |
| China Kadoorie Biobank I63 | 14,302 | 67,954 | Ischemic-stroke code but participant overlap with earlier stroke consortia remains a problem |

Public pages: https://pheweb.ibms.sinica.edu.tw/pheno/433.21 , https://pheweb.ibms.sinica.edu.tw/pheno/433.2 , https://pheweb.ibms.sinica.edu.tw/pheno-filter/433.3 and https://pheweb.ckbiobank.org/pheno/i63 .

**Important**: no complete original outcome GWAS was downloaded; variant-level beta/SE, effect allele, imputation INFO, per-SNP effective N or full-file checksum were not checked. Merely finding a download link does not establish data access, complete variant schema, or license compliance. 2026 TPMI published original GWAS is not assumed identical to these public 2025 cohort PheWeb phenotypes.

### Distinguish original association statistics from polygenic score weights

PGS Catalog https://www.pgscatalog.org/score/PGS019250/ reports a score named TPMI_433.21_Lassosum2 and https://www.pgscatalog.org/score/PGS019234/ reports TPMI_433.2_SBayesR, both GRCh38. These are **polygenic score weights, not original SNP-level outcome GWAS beta/SE/P summary statistics**. The Catalog lists 30,664 association source cases for PGS019250, which differs from the live TPMI 433.21 phenotype page's 9,249 case count. Treat this as a source version/phenotype-cohort mapping issue to reconcile; do not conflate PGS evaluation, training, association and public PheWeb case-control cohorts.

### 2026 Taiwanese early-onset ischemic stroke study

Neurology 2026, DOI https://doi.org/10.1212/WNL.0000000000218126 , reports **5,546 cases and 143,017 controls** in a TPMI-based early-onset IS analysis and loci including rs541118668 near CYP4F3. This is **not established as the same summary-statistics dataset as PheWeb 433.21**. Its original full downloadable GWAS output, detailed cohort participant disjointness and SNP-wise N_eff remain unverified. The separate clinical TOAST-characterized patient group used for follow-up does not automatically redefine the public ICD-code-based discovery GWAS as TOAST-adjudicated.

### Other candidates and study origin risk

BBJ: OpenGWAS bbj-a-129 https://gwas.mrcieu.ac.uk/datasets/bbj-a-129/ lists ischemic stroke with 17,671 cases, 192,383 controls and 8,885,705 variants, but OpenGWAS VCF/API access can require a token and BioBank Japan participates in earlier major stroke consortia. CKB: official PheWeb I63 has 14,302 cases and 67,954 controls, but CKB is also represented in published EAS stroke meta-analysis; it cannot be called an independent GIGASTROKE replication without disjointness evidence.

GBMI 2023 https://doi.org/10.1016/j.xgen.2023.100345 has a 16-biobank meta-analysis of ischemic stroke (60,176 cases/1,310,725 controls) with BBJ and other overlapping data sources. The paper explicitly used a BBJ-excluded GBMI dataset when combining with MEGASTROKE due to overlap. A meta-analysis from GBMI or MEGASTROKE is not automatically independent from the original GIGASTROKE GWAS.

Published EUR GIGASTROKE original corrected supplement specifically includes FinnGen-named contributing cohorts (8,851 original AIS cases from two named source cohorts, or 14.25282% of reported EUR AIS cases); this is study-origin evidence and not a count of overlapping participants in FinnGen R13. See [FinnGen cohort audit](ckm_stroke_finngen_r9_r10_r13_overlap_audit_20261121.md).

## Source and scientific GO / HOLD matrix

| Gate | Status |
| --- | --- |
| Public TPMI 433.21/433.2/433.3 and CKB I63 phenotype metadata | PASS |
| PGS weights distinguished from original GWAS associations | PASS |
| Exact 2026 Taiwan early-onset GWAS equals 2025 TPMI PheWeb phenotype data | HOLD |
| Original whole-outcome summary statistics downloaded, source allele/beta/SE/P verified | NOT RUN |
| Original phenotype case definition independently validated against TOAST/MRI/CT | HOLD |
| Individual-level disjointness between TPMI and GIGASTROKE | HOLD |
| Original GIGASTROKE EUR effective N per variant | HOLD |
| New LDSC h2/intercept/rg, MR, colocalization | NOT RUN |

## Next practical scientific gates

1. Verify a usable official TPMI 433.21 association summary statistics download path, permissions/terms, original filename and header using a bounded source audit. The published site warns of download limits. Never repeatedly bulk-download or bypass restrictions.
2. Map public TPMI case ICD codes, controls and missingness for 433.21 to clinically defined AIS and ensure the 2026 Neurology GWAS is a different analysis if counts/build/source differ.
3. Confirm study disjointness from GIGASTROKE by comparing original contributing cohorts, calendar dates and overlap disclosures; do not claim that different biobank brand names mean disjoint participant IDs.
4. Prepare GRCh38-to-GRCh37 variant normalization, effect allele alignment, EAS-specific LD reference and candidate cis-pQTL MR/coloc gate separately. Do not use an arbitrary pooled case-control N for primary EUR LDSC.

Automation code: scripts/ckm_stroke/audit_eas_tpmi_ais_phenotype_sources.py; tests/test_ckm_eas_tpmi_ais_phenotype_sources.py; .github/workflows/ckm-eas-tpmi-ais-phenotype-source.yml. Draft PR #284 remains unmerged.


## Live official TPMI 433.21 original GWAS HEAD-only accessibility gate

A second source-backed [GitHub Actions #44871483498](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44871483498) completed successfully with **14/14 negative and positive contract tests PASS**; [artifact #14694043597](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44871483498/artifacts/14694043597). Only official original PheWeb phenotype HTML was downloaded (small pages), and **HTTP HEAD** rather than GET for the actual large GWAS result files.

| PheWeb 433.21 original association download | Advertised original source URL | HEAD | Reported compressed bytes | Result body downloaded? |
| --- | --- | --- | ---: | --- |
| Summary Statistics (PLINK2) | https://pheweb.ibms.sinica.edu.tw/download/gwas_plink2/433.21 | 200 | **58,302,393** | NO |
| Summary Statistics (SAIGE GWAS) | https://pheweb.ibms.sinica.edu.tw/download/saige_gwas/433.21 | 200 | **816,992,833** | NO |

Both HEAD requests remained on the official source hostname. HTTP 200 proves an advertised endpoint responded, **not** the file's internal format, full integrity, per-variant effect and sample-size columns, legal access authorization, or consistency with the 2026 Neurology early-onset IS discovery sample.

Important distinction: 2025 TPMI PheWeb **433.21 9,249 cases and 304,660 controls** versus 2026 published TPMI early-onset IS **5,546 and 143,017**. The PGS Catalog uses yet another association count for a score of the same named phecode; none is promoted as original paper-to-portal equivalence.

**Recommended next actual analysis gate:** obtain at most the smaller 58 MB PLINK2 source once, strictly respecting TPMI published download limits and terms, record response SHA256, verify the header and GRCh38 REF/ALT/BETA/SE/P/AF/missingness conventions, and compare the 433.21 ICD case definition. If the user opts to proceed with allowed downloads, carry out one reproducible CI or Drive-backed ingest. Do not run LDSC h2/rg or claim independent replication until overlap provenance and SNP-specific effective N are justified. The 817 MB SAIGE file is not required just to inspect access; avoid redundant downloads and do not bypass site throttling.
