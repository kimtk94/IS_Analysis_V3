# GIGASTROKE original EUR raw summary-stat integrity gate

Evidence date: 2026-11-16 KST. This extends the publisher metadata/GRCh37 phenotype mapping lock already recorded in `config/ckm_stroke/gigastroke_eur_official.tsv`. It does NOT imply harmonized original GWAS inputs or LDSC readiness.

## Verification contract

1. Exactly one of AIS, LAA, CE or SVD is selected. No uncontrolled URL or local archive path is accepted.
2. The selected accession, build, publisher MD5 and expected trait must match the local manifest.
3. The actual EBI publisher `*-meta.yaml` and `md5sum.txt` are separately retrieved and compared with manifest labels, raw filename, GRCh37, sample size and MD5.
4. Stream the publisher raw `.tsv.gz` directly into MD5 without storing it. Enforce a strict 400 MB maximum, HTTP Content-Length completeness, nonempty payload and gzip magic bytes. Report JSON only.
5. Only matching **entire raw bytes** can produce `PASS_RAW_SOURCE_INTEGRITY`; all errors are HOLD or FAIL. Scientific status remains HOLD in every case.

### CLI execution (from repository root)

```bash
python scripts/ckm_stroke/verify_gigastroke_raw.py --phenotype LAA --out artifacts/ckm_stroke/raw_integrity_LAA.json
```

All 4 studies need independent evidence, and actual run logs should be linked here as results become available. The manual GitHub Actions workflow `.github/workflows/ckm-stroke-gigastroke-raw-integrity.yml` accepts exactly one phenotype per dispatch. A push only runs small offline tests, not four repeated large downloads. Results are stored as short 30-day JSON artifacts; raw files must not be committed or uploaded as Actions artifacts.

## Initial status

- **Original official phenotype/build metadata:** PASS (AIS, LAA, CE, SVD; GRCh37).
- **Declared source MD5 ledger:** available; printed MD5 values are not independent file checksums.
- **Raw-byte integrity:** NOT_RUN until all source streams have completed and reported MD5.
- **Whole-genome source allele/variant QC, LDSC, MR and colocalization:** HOLD/NOT_RUN.

Technical limitation: MD5 compares published compressed bytes only. It cannot verify GWAS association semantics, participant/sample overlap, genome-wide SNP inclusion, cohort ancestry or scientific inference.

See the PR Draft #284 and issue #285 for work sequencing. Do not run full source downloads repeatedly while GitHub Actions queues are saturated.
