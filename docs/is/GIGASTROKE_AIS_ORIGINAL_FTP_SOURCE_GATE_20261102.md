# IS P3A — GCST90104545 original GWAS summary-statistics source gate (2026-11-02)

## Narrow scope
The official GWAS Catalog v2 metadata accession and ancestry gate passed 5/5:
[GitHub Actions 41433562541](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41433562541).
The current phase tests **original full-summary-statistics provenance** for
East Asian *any ischemic stroke* GCST90104545, not just accession metadata.

## Primary source of truth
The GWAS Catalog [summary-statistics documentation](https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics/)
states that v2 study REST metadata is distinct from original full GWAS
summary-statistics data, supplied through accession-binned FTP directories.
Every original GWAS summary-statistics file can carry adjacent YAML with
file_type, data_file_name, data_file_md5sum, case_count and control_count;
availability and exact field coverage depend on the individual study.

Expected public folder:
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST90104545/

These URL and common YAML/header conventions are **documentation-derived expectations** until actual GitHub runner network evidence confirms them. Do not mark original data accessible based on generated URL alone.

## Bounded original-data CI

- Script: scripts/is/audit_gigastroke_ais_original_ftp.py
- Tests: tests/test_is_gigastroke_ais_original_ftp.py
- Workflow: .github/workflows/is-gigastroke-ais-original-ftp.yml

The CI checks:
1. Official FTP index availability and exact accession file listing.
2. Unique original metadata YAML and declared original filename.
3. Original official YAML MD5 and accession, with selective source metadata.
4. HTTP HEAD size when available, bounded compressed prefix and decompressed
   first-line column names, including chromosome, base position, effect and
   other alleles, p-value, effect estimate, and standard error.
5. If the preflight is sufficiently complete, bounded full primary
   gzip download, original-file MD5 and SHA256, complete gzip/TSV scan and
   aggregate number of records. Hard cap: 1.5 GB compressed, 20 million rows.
6. Export **only aggregate JSON**, never the actual raw GWAS data or
   individual sample identifiers. Keep temporary data on GitHub runner only.

The automated full file audit will NOT authorize downstream MR/coloc, even
when file checksum/columns pass. Full-file authentic source proof cannot
automatically resolve study-to-canonical transform, GRCh37 reference FASTA,
signed-effect alignment, genotype LD or cohort overlap.

## Expected publication-level counts (not adopted as original N)

GCST90104545 is AIS, **not** GCST90104544 AS. The previously verified
phenotype reference reports 19,032 EAS AIS cases and 237,242 controls.
Those values are *published expectations*, not directly adopted into coloc.
Use YAML and original source-file metadata as a further line of evidence,
and do not override older GIGASTROKE EAS5 archived provenance automatically.

## Hold conditions

- No FTP listing or original file: HOLD_SOURCE_UNAVAILABLE.
- No uniquely identified accession metadata YAML: HOLD.
- YAML missing official MD5 or mismatched original-file bytes: HOLD.
- Required signed-effect/SE columns absent: HOLD.
- More than 1.5 GB compressed or 20 million rows: bounded HOLD, not a failure
  of the scientific hypothesis.
- Original GRCh37 REF/ALT orientation and comparison to locally archived
  canonical GIGASTROKE EAS-AIS: HOLD until exact original file bytes and
  source-to-canonical transformation rules verified independently.
- MR, SuSiE/coloc, case fraction, and independent EAS replication remain HOLD.


## Observed completed official original file audit — 2026-11-02

[GitHub Actions #41434867406](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41434867406) completed SUCCESS and actually downloaded and scanned the **entire** 260,689,208-byte official public file. These are actual data-derived results, not merely fixture tests.

| Verified source fact | Actual log |
| --- | --- |
| Exact filename | GCST90104545_buildGRCh37.tsv.gz |
| Official original file data MD5 | 904e12e33ec3c2eeefe4dc93bec4995c |
| Independently computed compressed file SHA-256 | c521d5695414bdedcc84d2ea21d588e330c75a672f36a87c3fd092b4c55687a3 |
| Full original data rows | 6,725,216 |
| Compressed bytes | 260,689,208 |
| Filename companion YAML SHA-256 | 7383b859ce9bb4811c15bc634592ee1544788536797a5c663644e1aa648f1d3e |
| Original genome assembly according to YAML | GRCh37 |
| YAML file_type | pre-GWAS-SSF |
| YAML is_harmonised | false |
| Source file effect estimate and uncertainty | beta and standard_error columns present |
| Allele columns | effect_allele, other_allele, effect_allele_frequency |
| Additional numerical columns | chromosome, base_pair_location, p_value, odds_ratio, ci_lower, ci_upper |

This **row count exactly matches** the 2026-10-11 separate Drive full-cis study's GIGASTROKE EAS-AIS canonical GWAS count (6,725,216). It does **not** establish row-for-row data equality, canonical variant ID equivalence, signed beta normalization, or source-to-canonical transformation identity without an independent SNP-level original/canonical comparison.

YAML has an MD5 and genome assembly, but the originally selected YAML metadata keys did **not include file-specific case_count or control_count fields**, so the reported phenotype study N and case fraction remain reference-only. In particular, do not claim that SE numeric availability alone proves proper effective N for coloc.

### Added whole-file scientific diagnostics

The original whole-file stream was extended to compute a **quantitative** summary while reading all rows already verified by MD5: count finite beta / positive SE / valid p, EAF within [0,1], non-A/C/G/T or identical allele pairs, A/T and C/G palindromic SNPs, and a 1-in-101 deterministic two-sided Wald p-value consistency sample. An invalid/missing numeric result produces a separate QC HOLD; the checksum verification remains a distinct proof of original bytes.

Inspect the newer CI run's true logs before citing these quantitative diagnostics. Their implementation is not proof that the 6.7M-row values have already passed. Original original-file effect allele / original REF/ALT orientation, exact N, ancestral LD and downstream colocalization are **HOLD**.


## Actual whole-genome numeric GWAS audit — 2026-11-02

The updated source audit [Actions #41436048672](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41436048672) completed **SUCCESS** with all tests passing. It re-fetched and independently rechecked the same public MD5 and SHA-256 original source file; 6,725,216 rows processed.

| Measurement | Actual scanned count |
| --- | ---: |
| Finite beta, finite positive SE, valid p between 0 and 1 | 6,725,216 |
| Valid effect-allele frequency in [0,1] | 6,725,216 |
| A/T and C/G palindromic allele pairs (all SNP rows) | 1,028,002 |
| Non-palindromic or other allele pairs | 5,697,214 |
| Every 101st source row tested for Wald p consistency where p values nonzero | 66,586 |
| Wald log10(p) differences greater than 0.2 in comparable sampled rows | 0 |

A palindromic allele-pair flag is a **strand orientation risk**, not a discovery or proof of incorrectly mapped variants; reference FASTA and study allele-frequency semantics must be used before any ALT-normalized source harmonization. Exactly matching 6,725,216 source and historical canonical row counts is strong lineage evidence, **not** per-variant equality proof. Original file-specific study N, case fraction and ancestry-matched LD remain unapproved.

The observed numeric QC log is pinned without individual-level data in [research/is/evidence/gigastroke_eas_ais_numeric_qc_20261102.json](../../research/is/evidence/gigastroke_eas_ais_numeric_qc_20261102.json).

## Three-locus SNP-level source-to-Drive anchor check

A separate small, pre-specified set of positional anchors was taken from the connected [2026-10-11 four-gene QTL–GWAS report](https://drive.google.com/file/d/1sIko_lf94xHcENl-CUSYT_RmYbHjZN4N/view):

| Gene | GRCh37 position | Unordered base alleles from report | Prior EAS AIS p |
| --- | --- | --- | --- |
| FGF5 | 4:81,182,554 | T/C | 2.888e-8 |
| CALHM2 | 10:105,447,782 | A/G | 1.332e-6 |
| CALHM2 | 10:105,459,834 | T/C | 7.792e-7 |

The pipeline now checks whether the official public raw file contains **exactly one** row with each coordinate and allele set, with p-value within the original source-log tolerance (0.05 log10 p). This is a small independent source-to-historical-crosswalk spot check — **not a full mapping, GTEx QTL validation, independent replication or colocalization**.

The updated CI containing the spot-check must complete before promoting the spot-check result beyond HOLD. If p or alleles disagree, the appropriate interpretation is a source-lineage/variant-normalization blocker, not automatically a negative biological association.


## Actual 32/32 signed-ALT crosswalk validation — 2026-11-02

[GitHub Actions #41437991270](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41437991270) completed the actual raw-file verification and full original GWAS scan. All 19 CI regression tests PASS; historical original file SHA-256 reproduced exactly:
c521d5695414bdedcc84d2ea21d588e330c75a672f36a87c3fd092b4c55687a3.

The independently selected 32 **nonpalindromic** SNPs came from the connected original
[Drive four-gene full-cis SNP evidence table](https://drive.google.com/file/d/1eGY5PmwevFTCVcBTP9lqvXOprCAxZ7NU/view)
(Drive TSV SHA-256 2ddd5e80e84f2798200d2c926bb8833210dbc50fd4c2178edf8f2626bc097571, 6,363,722 bytes).
Selection is 4 genes × 2 historical allele-order categories
(FASTA_REF_ALT, OTHER_EFFECT_ORDER) × 4 positional quantiles.
The 32 positions are unique, with no palindromic A/T or C/G alleles.
**All 32/32** present exactly once in original public GCST90104545; all
have the correct unordered REF/ALT pair, same signed beta after orienting
the original raw effect allele to ALT, same ALT allele frequency and same
p-value (prespecified numeric tolerances).

The actual per-SNP check predicates and aggregate counters were captured in
[research/is/evidence/gigastroke_ais_original_drive_32_snp_live_result_20261102.json](../../research/is/evidence/gigastroke_ais_original_drive_32_snp_live_result_20261102.json).
The reference panel itself is stored in
[research/is/evidence/gigastroke_ais_drive_32_nonpal_snp_panel_20261102.json](../../research/is/evidence/gigastroke_ais_drive_32_nonpal_snp_panel_20261102.json).

These outcomes support **sampled allele/sign source-to-archive linkage**,
not whole-GWAS canonical equivalence, original molecular QTL model identity,
full GWAS case N, independent BBJ replication, or colocalization.
In particular 32 nonpalindromic examples do not automatically
disambiguate the previously flagged 1,028,002 palindromic original rows.
**Scientific status HOLD** until full reference, cohort and QTL/LD evidence gates.


## Actual complete four-gene EBI cis × official AIS GWAS source reconciliation

The independently generated Drive source from 2026-10-11 has **20,689 gene–SNP records** (4 genes; not 20,689 unique SNPs), including **12,223 gene–SNP GWAS intersection rows**. After deduplicating repeated variant positions across genes with *zero conflicting repeated statistics*, the exact GIGASTROKE EAS-AIS comparison universe is **6,273 unique GRCh37 SNPs**:

- chromosome 4: **3,227**; chromosome 10: **3,046**.
- nonpalindromic: **5,326**; palindromic A/T and C/G: **947**.
- input freeze: [chr4](../../research/is/evidence/gigastroke_ais_drive_fullcis_chr4_20261102.tsv) and [chr10](../../research/is/evidence/gigastroke_ais_drive_fullcis_chr10_20261102.tsv), derived from connected Drive TSV SHA-256 `2ddd5e80e84f2798200d2c926bb8833210dbc50fd4c2178edf8f2626bc097571`.
- actual live verification: [GitHub Actions #41439023099](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41439023099) **SUCCESS**, 23/23 unit tests PASS, original FTP file SHA-256 `c521d5695414bdedcc84d2ea21d588e330c75a672f36a87c3fd092b4c55687a3` independently rechecked on the 6,725,216-source-row original file.

**Result of full-locus original-vs-Drive GWAS numeric test:**

| Stratum | Verified | Failed |
| --- | ---: | ---: |
| Nonpalindromic SNP position, allele pair, p, ALT-signed beta, ALT-AF | 5,326/5,326 | 0 |
| Palindromic SNP numeric matching, orientation NOT independently resolved | 947/947 | 0 |
| **Combined source-to-Drive crosswalk** | **6,273/6,273** | **0** |

All comparisons obeyed preset numerical tolerances, not byte-for-byte floating-point equality. The 947 palindromic SNPs matched the archived numeric encoding, **but this does not establish independent DNA strand resolution**. Their effect/REF/ALT strand orientation remains REVIEW until reference-genome and frequency-aware confirmation.

Result artifact:
[research/is/evidence/gigastroke_ais_original_drive_fullcis_6273_live_result_20261102.json](../../research/is/evidence/gigastroke_ais_original_drive_fullcis_6273_live_result_20261102.json).

**Scientific boundaries:** This validates the *GWAS* source lineage for all 6,273 SNPs intersecting the **separate EBI-reprocessed GTEx QTD000161 cis tested universe of four genes**. It is not proof that (a) this is the original GTEx v8 FastQTL model used by historic 646 coloc entries, (b) the other 6.7 million source GWAS rows or the entire 2,225-gene universe underwent full canonical equivalence tests, (c) GWAS case fraction/effective N is file-specific and approved, (d) appropriate LD is available, or (e) any causal target has been established.

**Final scientific status remains HOLD** for full P3A / MR / coloc; **source-to-historical-GWAS numerical crosswalk PASS** within explicitly defined 6,273 SNP universe.
