# IS P3A — Independent official ROSMAP QC report vs 2021 metadata vs Release7 N

**2026-11-23 KST · MASTER DEGREE · ischemic stroke only.**

## Source-level question and what changed

The previous [merged PR #319](https://github.com/kimtk94/IS_Analysis_V3/pull/319) independently verified 2021 ROSMAP source count matrix 581 pseudonymous columns with **569/569** joint RNA+genotype QC-pass samples retained; all 15 authentic ROSMAP genes are present. This leaves a **nine-sample numerical difference** against eQTL Catalogue Release7 QTD000434 nominal sample size **560**.

A separate official [EMBL-EBI eQTL Catalogue ROSMAP QC report](https://www.ebi.ac.uk/eqtl/static/QC_reports/ROSMAP_QC_report.html) explicitly publishes the following **historical source QC-stage numbers**:

| Independently named source / reported stage | Count | Evidence status |
| --- | ---: | --- |
| Official ROSMAP QC report, before quality control | **581** | Direct official QC report |
| Official ROSMAP QC report, after quality control | **576** | Direct official QC report |
| Report genotyped individuals in ROSMAP | **581** | Direct official QC report |
| Report matching genotype-expression | **576** | Direct official QC report |
| Report nonmatching genotypes | **0** | Direct official QC report |
| 2021 Zenodo ROSMAP.tsv, both RNA+genotype QC passed | **569** | Independent publisher MD5-verified file |
| Release7 QTD000434 ROSMAP brain_naive nominal sample N | **560** | Official eQTL Catalogue Release7 metadata |

The numerical differences **576 − 569 = 7**, **569 − 560 = 9**, and **576 − 560 = 16** are **descriptive differences between sources**, not a certified sequence of subject exclusions. Do NOT label five, seven, nine or sixteen *particular donors* as removed at any stage without source-comparable sample IDs, exact release version, genotype/phenotype intersection and covariate matrix.

### Why this matters

- The historical ROSMAP QC report supplies an **independent official intermediate sample-stage reference of 576**; it is not automatically the final Release7 qtlmap input.
- The report itself says matching genotypes 576 and nonmatching 0, and historical QC-after 576. This is not direct evidence of why Release7 qtlmap declares 560.
- Previous original 2021 count matrix had all 569 joint metadata-QC-pass columns. Missing 2021 count-matrix columns cannot explain nine.
- The official [eQTL Catalogue study page](https://www.ebi.ac.uk/eqtl/Studies/) also lists historical ROSMAP study **576** samples / donors; its release correspondence to the 2021 file or later R7 qtlmap has not been established.
- [eQTL Catalogue Release notes](https://www.ebi.ac.uk/eqtl/Release_notes/) document reprocessing and normalization changes between releases. The [official qtlmap pipeline](https://github.com/eQTL-Catalogue/qtlmap/blob/master/docs/usage.md) requires intersecting normalized expression, sample_meta and VCF sample IDs (with covariates), which are not recovered from these aggregate records.

## Reproducible verification, scope and privacy

- Script: \`scripts/is/audit_p3a_rosmap_official_qc_report_stages.py\`
- Workflow: \`.github/workflows/is-p3a-rosmap-official-qc-report-lineage.yml\`
- Tests: \`tests/test_is_p3a_rosmap_official_qc_report_stages.py\`.
- Fetch live original QC report from EMBL-EBI HTTPS, parse only those five specifically named reported values, SHA-256 hash raw report.
- Fetch 2021 publisher Zenodo record, enforce strict deposited **ROSMAP.tsv** MD5/byte size, count private in-memory joint QC PASS IDs, save only count and source hash.
- Independently check official Release7 **QTD000434 ROSMAP brain_naive N=560** and its source hash.
- Produce 36-group HOLD JSON and explicit fields distinguishing independent evidence from proven sequential filtering.
- Preserve raw donor IDs and any RNA matrix out of Git history and Actions artifact; no source-specific per-donor QC table, genotype data, original counts or phenotype values exported.
- No guessed gene-level QTL regression N, sdY=1, GWAS beta log(OR), new formal coloc.abf or coloc.susie. All 36 remain HOLD until authentic source inputs are established.

## Remaining study gates

1. Identify and independently verify **exact Release7 QTD000434 qtlmap** normalized gene expression matrix, sample_meta, genotype matching and covariate inclusion. Determine whether it is publicly reachable or controlled-access. Any availability determination must cite the exact publisher source inventory, not an absence of search results.
2. Derive defensible **15 ROSMAP gene-specific sdY** or a fully justified quantitative-QTL alternate N/MAF model from appropriate original source. Audit GTEx 21 groups on the same scientific standard.
3. Separately confirm GIGASTROKE EAS/EUR AIS original beta release units and replication independent-cohort inputs. Formal coloc remains 0/36 until scientific evidence gate clears.

Technical CI PASS is not scientific Coloc approval. See [Issue #307](https://github.com/kimtk94/IS_Analysis_V3/issues/307).
