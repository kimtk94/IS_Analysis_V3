# IS P3A — original 36-assay quantitative QTL variance diagnostics and case-control Coloc input pathways

**Date:** 2026-11-23 (KST). Scope **MASTER DEGREE / ischemic stroke**, **36 original gene×tissue molecular QTL assays, 20 unique genes, one chr10 G0015 GWAS region**. These are not 36 independent genomic signals. Scientific Coloc is still **0/36 eligible and 0/36 executed**.

## Original source and actual independent CI

[GitHub Actions original source CI #45038801844](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45038801844): **SUCCESS**, 15/15 focused tests, repository 379 total tests (378 PASS, 1 Rscript SKIP).

The CI independently re-read the original indexed eQTL Catalogue **QTD000131 GTEx Artery Aorta, QTD000171 GTEx Brain Cortex, QTD000434 ROSMAP DLPFC** in the original 20-gene Ensembl105 reference union chr10:101837472–105353833 (GRCh38). It independently verified all three complete source-region stream SHA-256 hashes against the earlier main verified source, exactly reproducing the **148,344** bona fide unique-within-gene-assay nonpalindromic SNP observations, their full source **AN/2 per-SNP histograms, MAF discrepancy counts**, and all 36 original phenotype/gene IDs. No 11GB full-QTL bulk download or original per-SNP effect file has been retained.

Earlier evidence incorporated and enforced by the workflow:

- [Original complete cis-permutation n_variants CI #45037127536](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45037127536): **36/36 gene-assay tested-cis n_variants counts matched**, all original SNP/INDEL tests tracked.
- [Original GWAS–QTL GRCh37/38 SNP signed-allele CI #45032828465](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45032828465): EAS **89,597** and EUR **128,282** SNP×gene-assay signed allele matches (descriptive records, not distinct SNPs).
- [Independent QTL AN/2 CI #45034385821](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45034385821): GTEx aorta **366–387**, GTEx cortex **186–205**, ROSMAP **560**, 45 original ROSMAP MAF-vs-AC/AN difference records remain scientific REVIEW.
- [Original GIGASTROKE 2022 AIS discovery Supplementary Table 1 CI #45035814621](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45035814621): EAS AIS 19,032 cases / 237,242 controls, EUR AIS 62,100 cases / 1,234,808 controls, and BBJ within EAS discovery (17,651 cases and 159,610 controls), hence **BBJ not independent EAS replication**.

## What the coloc package actually requires

Official [coloc v6.0.1 data structures vignette](https://chr1swallace.github.io/coloc/articles/a02_data.html):

1. A **quantitative QTL** dataset using **original signed beta and varbeta=SE²** can be passed with a **direct, reliably supported quantitative phenotype sdY** without requiring an N for that specific input pathway. If the actual phenotype sdY is unavailable, coloc can estimate it from variance-of-beta, MAF and a scientifically defensible N. The latter should not be filled from a uniform top-line cohort N if original genotype/effective sample coverage varies without a documented model interpretation.
2. A **case-control GWAS** effect-size dataset with original beta and varbeta=SE² has a valid **type=cc** minimum-data pathway that **does not require N or case fraction s as separate scalar arguments** (in contrast to p-value-only inputs that do require N, s and MAF). This avoids mistakenly asserting that an arbitrary GWAS N_eff is a mathematical prerequisite for every Coloc.abf path. However the **original GWAS beta log-odds units**, per-SNP meta-analysis model, matching region, allele orientation, and sample overlap remain scientifically important; earlier technical GWAS source matches do not independently approve these model assumptions.
3. coloc.abf has **one causal variant per trait per region** as its base assumption; coloc.susie requires compatible signed LD, potentially different matrices for eQTL and GWAS ancestry, and separate QC.

Official [eQTL Catalogue normalisation methods](https://www.ebi.ac.uk/eqtl/Methods/) document cqn then inverse-normal transformation of gene expression phenotypes. That **suggests** near-standardised quantitative phenotypes, but it does NOT independently verify the actual sdY for each of the 36 specific study×gene phenotypes. Do not silently set all sdY=1 until the source-specific transformed expression values or authoritative release metadata establishes this.

## Original raw QTL effect/SE/MAF/AN diagnostic

Calculated for **every original canonical SNP×gene-assay record** passing the earlier nonpalindromic, canonical allele, beta/SE/MAF and source reference verification.

**Descriptive single-SNP residual-standard-deviation proxy** (NOT direct phenotype sdY):

- Formula: proxy = original SE(beta) × sqrt[ 2 × original molecular QTL MAF × (1-MAF) × original per-variant (AN/2) ].
- This follows a simplified additive linear regression with genotype variance near 2f(1-f), and **ignores phenotype and genotype PC covariates, genotype uncertainty, imputation and residual finite-sample effects**. It is a **diagnostic of scale only**. It is neither a measured phenotype sdY nor an independently approved input to coloc. No posterior inference or prior conversion is performed.
- All source variant N histograms are **identical to the independently stored original sample-number audit**. The residual proxy is **calculated from original beta SE and MAF/AN**, not synthesised from known RNA normalisation assumptions.

| Original source | Real QTL gene×tissue assays | Original strict nonpal SNP×assay records | Groups with varying AN/2 | Source AN/2 min–max | Fraction original SNP×assay records with N equal to dataset catalogue N | Median proxy across gene groups (descriptive median-of-medians) | Group-median proxy min–max |
|---|---:|---:|---:|---:|---:|---:|---:|
| QTD000131 GTEx Aorta | 10 | **44,693** | **10/10** | **366–387** | **98.86%** | **0.7613** | **0.7251–0.9625** |
| QTD000171 GTEx Cortex | 11 | **45,570** | **11/11** | **186–205** | **99.04%** | **0.8973** | **0.5690–1.0072** |
| QTD000434 ROSMAP DLPFC | 15 | **58,081** | **0/15** | **560–560** | **100%** | **0.7259** | **0.3884–0.9067** |
| **Total** | **36** | **148,344** | **21/36** | **varies** | **variable** | — | **0.3884–1.0072** |

The per-source group-median ranges are not effect sizes, p-values, disease associations or rankings. The proxy can differ from 1 due to covariate-adjusted residual variance **even if the underlying transformed expression phenotype truly has sdY ≈ 1**. Thus these values cannot be used to infer actual phenotype sdY is not 1. They do demonstrate that substituting a guessed single residual-variance number is methodologically inappropriate. Fifteen ROSMAP groups share genotyped AN/2=560, but this **does not confirm per-gene phenotype expression/covariate effective N**.

### All 36 per-group proxy medians (context/QC, not association strength)

| Dataset | Gene | Median proxy | Dataset | Gene | Median proxy |
|---|---|---:|---|---|---:|
| Aorta | SLK | 0.8243 | Aorta | GSTO2 | 0.8040 |
| Aorta | NT5C2 | 0.7467 | Aorta | NEURL | 0.8909 |
| Aorta | CALHM2 | 0.7670 | Aorta | GSTO1 | 0.7280 |
| Aorta | PDCD11 | 0.7312 | Aorta | SFR1 | 0.7251 |
| Aorta | C10orf32 | 0.9625 | Aorta | RP11-416N2.4 | 0.7556 |
| Cortex | GSTO2 | 0.9211 | Cortex | NT5C2 | 0.8720 |
| Cortex | CALHM2 | 0.9122 | Cortex | CYP17A1 | 1.0072 |
| Cortex | GSTO1 | 0.7820 | Cortex | ITPRIP | 0.5690 |
| Cortex | C10orf32 | 0.8973 | Cortex | WDR96 | 0.8095 |
| Cortex | RP11-127L20.3 | 0.9688 | Cortex | RP11-724N1.1 | 0.9869 |
| Cortex | RP11-416N2.4 | 0.8449 | ROSMAP | SLK | 0.5617 |
| ROSMAP | GSTO2 | 0.7875 | ROSMAP | NEURL | 0.5469 |
| ROSMAP | CCDC147 | 0.8748 | ROSMAP | CALHM2 | 0.8788 |
| ROSMAP | INA | 0.3884 | ROSMAP | GSTO1 | 0.7238 |
| ROSMAP | ITPRIP | 0.6101 | ROSMAP | CNNM2 | 0.5437 |
| ROSMAP | C10orf32 | 0.7259 | ROSMAP | USMG5 | 0.7044 |
| ROSMAP | WDR96 | 0.9067 | ROSMAP | AS3MT | 0.8107 |
| ROSMAP | RP11-127L20.3 | 0.7562 | ROSMAP | RP11-724N1.1 | 0.7738 |

### Actionable approval sequence, not an arbitrary statistical guess

- **Exposure quantitative path:** obtain official source study-specific QTL input molecular expression matrix (normalised qcnorm results) and sample_meta/covariate matching for the 3 datasets; measure direct per-gene transformed expression sdY, or document a faithful regression-compatible sample-size methodology with original genotype dosages, missingness and phenotype inclusion. No default sdY=1 merely because inverse-normal was used.
- **Outcome CC signed effect-size path:** validate the publisher GWAS beta units (log OR / correct additive case-control meta-analysis), the source variant beta/SE² correspondence, and allele matching across one scientifically defensible shared cis region; **do not invent GWAS per-SNP N_eff if beta/varbeta route is chosen**. Publication-discovery cohort counts remain valid descriptive metadata, BBJ overlap remains explicit. If a p-values-only path is used, it has different N/s requirements.
- **Formal Coloc:** apply exact gene×tissue variant locus boundaries, coverage, allele QC and prior sensitivity once both actual input datasets are scientifically signed off. Multiple causal signals demand SuSiE and ancestry-appropriate LD, not invented or reused EAS/EUR population matrices.
- **Readiness status:** QTL sdY not directly observed, QTL model-specific N unapproved, GWAS log-odds units unreviewed, 8 EUR ROSMAP allele mismatches and INDEL/palindromic harmonisation retained; **coloc.abf 0/36 eligible/executed, coloc.susie 0/36 eligible/executed, no PP.H4**.

Only three source real QTL Tabix region reads and one GitHub Actions aggregate artifact (source hashes, 36 group sample histograms and descriptive proxy statistics) were created. Scripts: scripts/is/audit_p3a_36_qtl_sdY_proxy_input_gate.py; tests/test_is_p3a_36_qtl_sdY_proxy_input_gate.py; workflow .github/workflows/is-p3a-36-qtl-sdy-proxy-n-model.yml. Scientific blocker tracking [Issue #307](https://github.com/kimtk94/IS_Analysis_V3/issues/307) remains OPEN.
