# IS P3A ROSMAP — 2021 processed count matrix vs QC sample intersection audit

**Scope:** Ischemic stroke P3A 36 original gene×tissue assays (20 unique genes, one chr10 locus). **2026-11-23.** This branch extends the merged PR #318 source verification and does not change the formal Coloc gate.

## Source and experimental question

Previous original source audit found **581 distinct ROSMAP DLPFC pseudonymous metadata IDs**, **569 passing RNA and genotype QC**, while the independent **eQTL Catalogue Release 7 QTD000434** reports **N=560**. This nine-sample numerical discrepancy is currently **unattributed**. It is not legitimate to conclude that nine specific donors were removed because of an unstated filter.

The original 2021 eQTL Catalogue Zenodo record [4678936](https://zenodo.org/records/4678936) contains **ROSMAP_merged_counts_pseudo.tsv.gz** (25,484,860 bytes; publisher MD5 `d215dc801f687e15d5b938356d80d183`), in addition to the previously verified 53,064-byte **ROSMAP.tsv** metadata. The new workflow downloads both actual publisher files, checks their exact published MD5, validates the Release 7 official QTD000434 tissue-size record, and performs private in-memory set intersections:

1. 2021 processed count matrix pseudonymous sample **columns** present in the 2021 sample metadata.
2. The intersection with **2021 RNA+genotype QC PASS** samples (without exporting IDs).
3. The 2021 QC-pass samples **missing from the count-matrix header**, and matrix columns whose metadata fail QC.
4. Presence of the **15 authentic ROSMAP gene IDs** in the original processed count matrix, without exporting source expression values.
5. Explicit **36/36 scientific HOLD** output, regardless of whether the observed 2021 matrix sample intersection is 560, 569, or another count.

## Critical interpretation boundary

The 2021 Zenodo **processed gene counts** are *not* the exact **Release 7 qtlmap normalized expression phenotypes**. The 2021 sample intersection cannot alone establish the cause of the nine-donor difference, original qtlmap genotype-VCF/sample_meta/covariate intersection, any individual gene-specific regression N, or directly measured normalized quantitative phenotype `sdY`. Even an observed 2021 matrix intersection numerically equal to 560 remains descriptive only.

Per [official qtlmap documentation](https://github.com/eQTL-Catalogue/qtlmap/blob/master/docs/usage.md), the actual association workflow intersects the normalized phenotype input and sample_meta and uses genotype data. Release 7 association statistics use a later processing release; the public 2021 matrix must not be silently substituted for the association-stage phenotype matrix.

**No formal coloc.abf / coloc.susie, no PP.H4, no fabricated sample N, no sdY=1 and no assumption that AIS GWAS beta is exact log(OR).** Previously verified GWAS beta unit and independent replication source gaps remain open.

## Execution / acceptance criteria

- Code: `scripts/is/audit_p3a_rosmap_2021_count_matrix_intersection.py`
- Tests: `tests/test_is_p3a_rosmap_2021_count_matrix_intersection.py` (synthetic ID privacy, real source checksums, mismatched headers, malformed rows, wrong QTD/N, status HOLD)
- GitHub Actions: `.github/workflows/is-p3a-rosmap-2021-count-matrix-intersection.yml`
- Source: publisher-record-pinned MD5 and byte size; Release 7 original metadata SHA-256 preserved.
- Artifact: aggregate-only JSON of sample counts, target-gene presence and 36 scientific HOLD decisions; no original sample/donor identifiers, expression values or count matrices saved.
- After workflow success, attach **the actual Actions run ID and aggregate observations** to [Issue #307](https://github.com/kimtk94/IS_Analysis_V3/issues/307) before any claims about sample intersection.
- Next original-source evidence needed: Release 7 qtlmap-normalized expression, matching sample_meta/genotype/covariates (or a justified N/MAF-based quantitative Coloc input path), plus GIGASTROKE released beta effect-scale provenance.

Scientific meaning: **technical CI SUCCESS is not equivalent to formal Coloc approval**.


## Executed original-source result — 2026-11-23

[Original 2021 matrix audit workflow #45150798976](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45150798976) **SUCCESS** after correcting a brittle *whole Zenodo API JSON hash* gate. Whole record JSON fingerprint changed from the previously observed 818a5ebd… to **b8556d82b90be276c638c914f09ac8975855673d32faaecd19b3f0cee9060aab**; the deposited ROSMAP source *file checksums and sizes themselves* were independently rechecked and agreed with the publisher. The workflow now validates original record identity, 32-file inventory, and publisher MD5/length, while recording (not falsely requiring) the mutable record API JSON SHA.

- 2021 ROSMAP metadata: **581 unique source IDs; 569 RNA+genotype QC PASS; 12 other**.
- Real original processed count matrix: **581/581 metadata source IDs present as columns**; **569/569 QC PASS** included; **12/12 QC-not-PASS** also present; **zero of the 569 QC-pass source IDs absent**.
- Original 2021 processed count matrix feature rows: **58,870**; all **15/15** authentic original ROSMAP genes found in source matrix gene IDs.
- Original count matrix file: **25,484,860 compressed bytes; verified publisher MD5 `d215dc801f687e15d5b938356d80d183`; measured compressed SHA256 `17a13cc0a0f47e96537b2242c592ab4ce388b617510e53ff8a42e9e66e9b5f73`**.
- Release7 original dataset metadata QTD000434: **560**. The original 2021 count matrix **does not explain the nine-sample difference through absent matrix columns**; any downstream sample exclusion, genotype intersection, covariate completeness, phenotype selection or release-stage mismatch remains **unattributed**.
- Actions whole-repository test suite: **422 total, 421 pass, 1 unrelated Rscript skip**. The 12 new synthetic tests also passed.
- The aggregate artifact contains **no pseudonymous donor or sample identifier values, no original RNA expression values**, only counts, gene IDs and source hashes.
- **Coloc remains 0/36 approved or executed**. The 2021 count matrix must **not** be used to compute Release7 normalized phenotype sdY, and a 2021 matrix sample count does **not** certify a Release7 gene-specific regression N.

**Source-level next gate:** locate or independently reproduce the precise QTD000434 Release7 qtlmap normalized phenotype, sample_meta/genotype/covariate inclusion and gene-specific sdY evidence. Review GWAS beta unit provenance separately; do not call BBJ independent EAS replication of GIGASTROKE.
