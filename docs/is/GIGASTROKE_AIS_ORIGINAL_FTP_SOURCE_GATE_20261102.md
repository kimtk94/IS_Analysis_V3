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
