# CKM P1 BMI full original-source feasibility checkpoint

Validated in the Orca project worktree on 2026-11-16. This page is an evidence record and operational runbook, not scientific approval for LDSC or MR.

## Completed local source streaming evidence

- Corrected 2018 official GIANT/Yengo EUR BMI original: `https://giant-consortium.web.broadinstitute.org/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz`
- Actual compressed bytes transferred: **46,754,937 / 46,754,937** (HTTP Content-Length matched)
- Complete gzip decompression and parse: **PASS**
- Variant records scanned: **2,336,269** (all 22 autosomes, none rejected by structural/number checks)
- 10/10 required original header fields: PASS.
- Min / max per-variant N: **485,648 / 795,640**.
- P-values underflowing to 0 in the exported representation: **44 rows**. Treat as precision loss, not exact zero p; do not invent exact values.
- Actual source SHA256: `0d6ed0ea97870916b830ccae349df94ca6f3cc68c025c79e784729af7f7136a4`
- Actual source MD5 (computed but **not** publisher-verified): `ccca6c58bc45a53aec16619497da3f00`
- Full JSON evidence locally: `artifacts/ckm_stroke/bmi_giant_full_stream_20261116.json`. This directory is gitignored; raw GWAS never stored.

## What PASS means

`PASS_STREAM_STRUCTURE_ONLY` asserts that the **entire original compressed byte stream** was received and parsed without detected number/allele/coordinate violations in the implemented checks. It does **not** independently verify the publisher checksum (none secured in this audit), exact rsID/coordinate duplication, allele frequency against full GWAS participants, population LD, sample overlap, cohort stratification, LDSC intercept, or any genetic signal. The study's true use-of-file terms still apply.

## Remaining scientific gates

1. Find or obtain independently published checksum; until then report `SOURCE_SHA256_RECORDED` not `PUBLISHER_MD5_PASS`.
2. Audit duplicated rsIDs and locus alleles and quantify S-LDSC usable HapMap3 SNP intersection; make a documented rule for 44 `P=0` values (do not silently convert to `1e-300`).
3. Validate `Tested_Allele` effect orientation and distinguish `Freq_Tested_Allele_in_HRS` from study-wide frequency. Harmonize against EUR reference LD with per-variant N.
4. Investigate sample overlap of GIANT+UKB and GIGASTROKE exposures/outcomes, model LDSC intercept and appropriate uncertainty; do not interpret rg before heritability/quality diagnostics.
5. Input approval gate for **SBP remains HOLD** (NHLBI GRASP historical redirect or OpenGWAS JWT); **T2D remains HOLD** (T2DGGI EUR download terms/permission). Distinct gates even if BMI technical scan passes.

## Reproducible execution

```powershell
py -3 -m unittest discover -s tests -p test_ckm_stroke_bmi_full_stream.py -v
py -3 scripts/ckm_stroke/audit_bmi_full_stream.py --out artifacts/ckm_stroke/bmi_giant_full_stream_20261116.json
```

GitHub Actions: `.github/workflows/ckm-stroke-bmi-full-source-qc.yml`.
Every push to the relevant source/test/workflow files runs a small offline regression. The original 46.75MB stream runs **only** on a user-dispatched workflow or on a deliberate push whose final commit message includes `[bmi-full]`. It stores **only** a JSON artifact (30-day retention).

## Orca agent readiness

- Orca installed on Windows, linked to `%USERPROFILE%\orca\IS_Analysis_V3`.
- Installed official `@openai/codex` CLI version **0.169.0**.
- Windows npm `codex.ps1` execution-policy collision was resolved by reversible rename `codex.ps1.disabled`, leaving `codex.cmd` intact. Windows execution policy was **not** weakened.
- Codex login status: **not logged in**. Orca terminal displays **Sign in with ChatGPT** option. Account owner must authenticate once. No credentials should ever be placed in GitHub.
