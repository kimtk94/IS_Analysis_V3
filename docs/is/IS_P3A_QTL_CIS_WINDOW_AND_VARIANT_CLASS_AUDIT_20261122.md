# IS P3A — official Ensembl 105 gene-cis reference and excluded-variant source audit

**Scope:** MASTER DEGREE / ischemic stroke, 2026-11-22. **36 authentic QTL assay groups; 20 genes; chr10 G0015**. This document is a **real-source technical QC audit**, not a completed Coloc or evidence of IS causality. Formal `coloc.abf` eligible/executed **0/36**; `coloc.susie` eligible/executed **0/36**, all HOLD.

## A. Real original eQTL Catalogue source: 34,735 excluded records explained

- [Final updated official remote Tabix source + variant-class CI #45027599911](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45027599911) **SUCCESS**, **19/19 new source tests PASS**, repo-wide **287 executed; 286 passed and one skipped**.
- Each of the three data sources was actually read again using one bounded GRCh38 chr10 PySam/HTSlib indexed query per source; 36 group keys preserved. The previous CI passed 645,751 region-wide association rows, 163,885 original candidate-gene association records and 111,178 QC-valid *nonpalindromic single-nucleotide gene-assay* records.
- The strict SNP-only screen excluded exactly 34,735 records. The **source's own `type` field identified 100% as `INDEL`**, and allele REF/ALT length confirmed all have multibase alleles. This is not evidence of corrupt genotypes or REF mismatches. Do **not** write these off as erroneous in the manuscript.

| QTL accession | Tissue | Excluded by SNP-only filter | Source `type=INDEL` | MULTIBASE REF/ALT | Noncanonical variant key | Non-ACGT base | REF=ALT |
|---|---|---:|---:|---:|---:|---:|---:|
| `QTD000131` | GTEx artery aorta | 4,958 | 4,958 | 4,958 | 0 | 0 | 0 |
| `QTD000171` | GTEx brain cortex | 4,558 | 4,558 | 4,558 | 0 | 0 | 0 |
| `QTD000434` | ROSMAP brain DLPFC | 25,219 | 25,219 | 25,219 | 0 | 0 | 0 |
| **Total** | Three source datasets | **34,735** | **34,735** | **34,735** | **0** | **0** | **0** |

**Caution:** True indel normalization, full REF length match against GRCh38, conversion to GRCh37, repeat/multiallelic handling and GWAS allele matching are **not yet verified**. The source's INDEL type is source metadata; treating these as benign source variants does not automatically approve them in the Coloc input. The existing nonpalindromic SNP assay counts exclude indels by design and are not a list of 111,178 independent genomic SNPs.

Source code: `scripts/is/audit_p3a_indexed_qtl_real_chr10.py` adds disjoint original variant exclusion reason classes and gene-assay observed position bounds. Six additional synthetic variant negative tests were added in `tests/test_is_p3a_indexed_qtl_real_chr10.py`.

## B. Authentic public Ensembl 105 gene expression phenotype-position annotation

[Official Zenodo eQTL Catalogue phenotype metadata record](https://zenodo.org/records/4715946), file `gene_counts_Ensembl_105_phenotype_metadata.tsv.gz`, published MD5 `b999b21dcffeecbfe0bf886a5a91997b`. It is only about 1.9 MB compressed. The actual GitHub Actions job downloads the file, **independently checks the original publisher MD5**, reads the official `gene_id / chromosome / phenotype_pos` fields, maps all 20 genes to GRCh38, and compares their **reference-based ±1,000,000 bp** expected cis spans to the previously audited GRCh38 **shared chr10:102,403,919–104,460,593** region.

The eQTL Catalogue 2023 paper reports ±1Mb cis association testing around the beginning of each gene (https://pmc.ncbi.nlm.nih.gov/articles/PMC10538656/). The current reference is **Ensembl release 105**, but each original QTL dataset's actual release-specific phenotype annotation identity has **not** been independently checked against that reference. These are **positional reference-coverage estimates**, not a claim of actual detected variant or complete source-specific nominal cis test coverage.

**Actual verified GitHub Actions [#45027498753](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45027498753): SUCCESS**, 12/12 new negative/positive tests, 281 repository-wide tests (280 passed, 1 skipped). The resulting original Ensembl105-based coverage for **all 20 genes**:

| Ensembl 105 gene | TSS/phenotype_pos GRCh38 | In prior common ~2.05Mb context? | Reference ±1Mb overlap |
|---|---:|---|---:|
| CALHM2 | 103,452,402 | **Entirely inside** | **100.00%** |
| PDCD11 | 103,396,626 | Left 7,293 bp not covered | 99.64% |
| ATP5MK | 103,396,492 | Left 7,427 bp not covered | 99.63% |
| NEURL1 | 103,493,705 | Right 33,112 bp not covered | 98.34% |
| NT5C2 | 103,277,605 | Left 126,314 bp not covered | 93.68% |
| INA | 103,277,138 | Left 126,781 bp not covered | 93.66% |
| ENSG00000273108 | 103,608,619 | Right 148,026 bp not covered | 92.60% |
| CNNM2 | 102,918,294 | Left 485,625 bp not covered | 75.72% |
| ENSG00000272912 | 102,914,585 | Left 489,334 bp not covered | 75.53% |
| SLK | 103,967,140 | Right 506,547 bp not covered | 74.67% |
| AS3MT | 102,869,470 | Left 534,449 bp not covered | 73.28% |
| BORCS7 | 102,854,259 | Left 549,660 bp not covered | 72.52% |
| CYP17A1 | 102,837,472 | Left 566,447 bp not covered | 71.68% |
| SFR1 | 104,122,058 | Right 661,465 bp not covered | 66.93% |
| CFAP43 | 104,232,364 | Right 771,771 bp not covered | 61.41% |
| GSTO1 | 104,235,356 | Right 774,763 bp not covered | 61.26% |
| GSTO2 | 104,268,873 | Right 808,280 bp not covered | 59.59% |
| ITPRIP-AS1 | 104,323,369 | Right 862,776 bp not covered | 56.86% |
| ITPRIP | 104,338,465 | Right 877,872 bp not covered | 56.11% |
| CFAP58 | 104,353,833 | Right 893,240 bp not covered | 55.34% |

**Critical coverage conclusion:** Only **1/20 genes (CALHM2)**, corresponding to **3/36 gene×tissue group keys**, has its *reference-predicted* ±1Mb span fully inside the old common query region. **19/20 genes and 33/36 groups** have at least one spatially missing reference-cis segment.

The **conservative union of the 20 Ensembl105 reference windows** would be **GRCh38 chr10:101,837,472–105,353,833**, or **3,516,362 bp** (inclusive), vs the original common query 2,056,675 bp. This requires a **566,447 bp left extension and an 893,240 bp right extension**. This is a **proposed retrieval interval**, not a proven dataset-specific tested cis interval or declaration that all genes have exactly 2,000,001 tested bases.

The directly saved 36-group audit contains actual `QTL_dataset, gene_id, predicted ±1Mb, expected positional coverage, observed candidate QTL records`, with every `coloc_abf_eligible=false`.

## C. Scientific next steps / acceptance gate

1. Confirm individual source study-specific phenotype-position/version against each of the three original QTL sources; check whether Ensembl105 positional reference is exact, and derive any necessary corrected cis bounds per gene.
2. If consistent, perform **one bounded, indexed union window fetch per dataset** around `chr10:101837472–105353833`, not the whole 11GB. Check all 36 gene×assay groups, original per-variant beta/SE/MAF/ALT, real spatial coverage and QTL sample ascertainment.
3. Normalize INDELs explicitly (reference FASTA, REF/ALT length, allelic representation and liftOver handling), count retained/rejected and do not introduce duplicate variant-IDs or strand direction errors.
4. Match **actual** GRCh37 publisher GIGASTROKE stroke GWAS per-SNP allele/beta/SE/EAF, verify precise cases/controls/case fraction and BBJ cohort overlap, with signed ancestral effect orientation before any `coloc.abf` run.
5. Only when all required full-cis and per-source N evidence is validated: formal `coloc.abf`, report PP.H0–PP.H4, H3/H4, priors and sensitivity. `coloc.susie` additionally requires appropriate ancestry-specific signed LD.

**Present scientific verdict: 36/36 HOLD; 0/36 ABF or SuSiE eligible/executed.** This cannot be described as proven stroke genes, Coloc H4 evidence or independence from BBJ.
