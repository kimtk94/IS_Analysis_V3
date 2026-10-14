# IS public nominal-QTL table integrity audit — 2026-10-14

## Provenance
- Source: `docs/data/is/IS_2225_GENES_PUBLIC_NOMINAL_QTL_READINESS.tsv`
- Branch: `research/is-broad-discovery-20261009`
- Git blob SHA: `a99539a73431dc208d010f1a09256e9aad414624`
- Method: fetched complete GitHub text and inspected rows using an in-memory tab-delimited parser. This is **not** a server execution of the Python validator or a scientific replication.

## Observed
- 2,225 data rows; 2,225 distinct stable gene IDs.
- Source statuses:
  - `NOT_IN_EAS_30_GROUPS`: 1,356
  - `RAW_NOMINAL_QTL_RECORDS_AVAILABLE_NO_COLOC`: 149
  - `PUBLIC_NOMINAL_QTL_REGIONAL_PILOT_NOT_EXECUTED`: 629
  - `NO_NOMINAL_ROWS_IN_TARGETED_SLICES_NOT_NEGATIVE`: 91
- Checks: duplicate gene IDs, positive-QTL rows missing source IDs, positive-QTL rows missing minimum p, zero-QTL rows with p, and out-of-range numeric p. **0 violations**.
- This result does not establish correct allele mapping, ancestry-matched LD, independent GWAS signals, MR, or causal genes.

## Immediate analysis queue
1. Extract 149 nominal-QTL-ready genes without hard-filtering the 2,225-gene discovery denominator.
2. Resolve source-specific variant IDs, genome build, REF/ALT, effect alleles, frequency, and sample overlap.
3. Harmonize GWAS and QTL variants; separate EAS/JPT from EUR reference LD.
4. For each locus record `MATCHED`, `SWAPPED`, `PALINDROMIC_UNRESOLVED`, `MISSING`, or `MISMATCH`.
5. Only compute coloc with validated variant overlap and explicit model assumptions; distinguish single-causal-variant ABF from multi-signal SuSiE.
6. Retain `NOT_ASSESSED` for the 629 pilot-unexecuted and 91 targeted-slice-missing genes.

## Status
- GitHub TSV consistency: PASS (limited checks)
- Python CI: NOT VERIFIED
- Server analysis: NOT EXECUTED
- GWAS/QTL allele harmonization: NOT VERIFIED
- Colocalization: NOT VERIFIED
