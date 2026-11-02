# IS P3A — GIGASTROKE EAS stroke subtype source identity audit (2026-11-02)

## Scope: five GWAS Catalog accessions

The original [Mishra et al., Nature 2022](https://www.nature.com/articles/s41586-022-05165-3)
GIGASTROKE collaboration reports multiple stroke outcomes and ancestry-specific
analyses. [Published EAS subtype counts](https://www.mdpi.com/2308-3425/11/8/247)
and [PGS Catalog accession-level source sample metadata](https://www.pgscatalog.org/score/PGS002725/)
give the following *reference expectations*, not independently verified
GWAS file-level N or case fractions:

| Accession | Phenotype | Cases | Controls | Total |
|---|---|---:|---:|---:|
| GCST90104544 | any stroke (AS) | 27,413 | 237,242 | 264,655 |
| GCST90104545 | any ischemic stroke (AIS) | 19,032 | 237,242 | 256,274 |
| GCST90104546 | cardioembolic stroke (CE) | 926 | 237,242 | 238,168 |
| GCST90104547 | large-artery stroke (LAA) | 1,735 | 237,242 | 238,977 |
| GCST90104548 | small-vessel stroke (SVS/SVD) | 5,532 | 237,242 | 242,774 |

CRITICAL: GCST90104544 and GCST90104545 MUST NOT be interchangeable.
The difference of 8,381 cases between AS and AIS is phenotypically meaningful.
Do not take N from the multi-ancestry GIGASTROKE overall GWAS or from
the different BBJ 2019 / Sakaue 2021 file.

## Distinct existing source in connected Drive

[2026-10-11 full-cis BBJ/GIGASTROKE manifest](https://drive.google.com/file/d/1DmacR2aa0HJqIOKoKRS2htL4LEzviQ11/view)
and [study report](https://drive.google.com/file/d/1sIko_lf94xHcENl-CUSYT_RmYbHjZN4N/view)
document original server-derived canonical GIGASTROKE EAS-AIS data:
- File label: GCST90104545_AIS_GRCh37.canonical.tsv.gz.
- Study-wide canonical GIGASTROKE EAS-AIS source row count: 6,725,216.
- Four QTD000161 EBI reprocessed GTEx cis-QTL gene tests contain
  6,273 position-mapped EAS-AIS GWAS rows, of which 2,918 were historically
  other:effect allele order and 3,355 REF:ALT order, FASTA37-anchored in
  this independent 2026-10-11 audit.
- Those data validate a **distinct four-gene crosswalk** and do NOT
  independently authorize the complete original GIGASTROKE source ZIP/file
  hash, signed allele convention, original SE, population-matched LD,
  per-study N/case fraction or coloc.
- EBI reprocessed QTL is not automatically the original historical
  GTEx v8 coloc QTL source.

## New CI job

- Script: scripts/is/audit_gigastroke_eas_v2.py.
- Unit tests: tests/test_is_gigastroke_eas_v2.py (9 synthetic controls).
- Workflow: .github/workflows/is-gigastroke-eas-v2.yml.
- Query 5 official GWAS Catalog REST API v2 study accession records, then
  their 5 published study ancestry endpoints, using bounded HTTPS JSON.
- Verify **exact study accession** + explicitly typed EAS ancestry +
  phenotype label when official fields are available; missing data is HOLD,
  not an invented PASS.
- Print only source status, response hash and reference aggregate counts;
  never output raw variant/sample identifiers.
- No unauthenticated access bypass, no full GWAS download, no server SSH.

GWAS Catalog retired REST API v1 in August 2026:
[official migration guidance](https://www.ebi.ac.uk/gwas/docs/news/rest-api-v2-migration-guide/).
The current v2 metadata API is **not** an API for full GWAS summary-statistics.
A verified study metadata response does not imply the availability or
validation of full genome-wide beta/SE/allele source files.

## Approval gates

1. Accession & ancestry/phenotype pairing — determine from real v2 CI artifacts.
2. Original compressed GIGASTROKE EAS-AIS file SHA256, exact variant
   position/build and effect/other allele orientation — HOLD.
3. Confirm actual original GWAS SE or defensible p/beta SE derivation,
   exact case/control N and case fraction for each of five outcomes — HOLD.
4. Independently compare canonical EAS-AIS to original source across
   a bounded set of SNP coordinates, excluding palindromic ambiguity — HOLD.
5. Verify cohort overlap with BBJ; BBJ is included in GIGASTROKE EAS
   discovery and CANNOT serve as independent replication without
   BBJ-excluded data — HOLD.
6. Run per-locus corrected MR and colocalization only after approved
   source/QTL allele/LD intersection; no PP.H4 or causal gene claim.

**Scientific status: HOLD.**
