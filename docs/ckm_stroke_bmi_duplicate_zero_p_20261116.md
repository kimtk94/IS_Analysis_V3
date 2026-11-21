# CKM × Stroke P1 — GIANT BMI 9 duplicates and 44 exported P=0 records

**Evidence date:** 2026-11-16. **Primary evidence:** original corrected GIANT/Yengo 2018 EUR BMI `.txt.gz` (46,754,937 bytes; full SHA-256 `0d6ed0ea97870916b830ccae349df94ca6f3cc68c025c79e784729af7f7136a4`). Script re-streams the original data twice when duplicates exist and verifies the complete compressed SHA-256 after **each pass**. No raw GWAS is retained or committed.

## Direct original-source observations

| Statistic | Actual count |
| --- | ---: |
| Total original BMI variant rows | 2,336,269 |
| Extra rows with a previously seen rsID | 9 |
| Distinct duplicate rsID groups | 9 |
| Extra rows with previously seen chr:pos | 9 |
| Same chr:pos under *different* rsID | 0 |
| Groups identical on **all ten raw fields** | 9/9 |
| Groups identical on ordered alleles, chr:pos, BETA, SE, P and N | 9/9 |
| Unique rsIDs after conservatively keeping one of each exact duplicate | **2,336,260** |
| Original P printed as zero | **44** |
| Full-source pinned-SHA verification | **2/2** |

The nine exact duplicate rsIDs were:

`rs11576885`, `rs13598`, `rs2160421`, `rs3780697`, `rs3819262`, `rs573071`, `rs627141`, `rs7591686`, `rs969206`.

### Deterministic deduplication decision

Each repeat is **byte-for-byte identical across all parsed fields**, including the original whitespace/rounding of BETA, SE, P, N and HRS allele frequency. There is no cross-rsID coordinate collision among these nine groups.

**PASS_TO_DEDUP_PREPARATION**: when generating a future derived LDSC input, keep the first original row number for each duplicated rsID and drop subsequent *identical* copies, log the 9 removed rows, and require post-transformation unique rsIDs. The final LDSC sumstats file has **not** been produced yet. Never silently retain the first row if discordant variants appear in a future original release.

### P=0 treatment decision — diagnostics, not source imputation

Original zero-P records: **44**. The source P remains the literal **0**; the audit **does not set it to 1e-300** or claim recovery of the original exact P. For *independent diagnostic purposes only*, compute signed `Z = BETA / SE` from the **rounded** source effect size/SE, and compute normal two-sided Wald `-log10(P)` using log-domain numerical calculations.

| Diagnostic | Value |
| --- | ---: |
| Approximate signed Z lower/upper | -46.875 / +46.8125 |
| Signed beta/SE Z positive / negative | 26 / 18 |
| Approximate normal-Wald `-log10(P)` range | **311.365 – 478.899** |
| Nonzero-P subsample (every 1,000th original row) | 2,336 variants |
| Median absolute difference between observed and beta/SE-inferred `-log10(P)` | 0.0101 |
| Maximum absolute difference in that same subsample | 1.4513 |

At extremely small P, the exported P value may underflow or be rounded. These Wald estimates are **not guaranteed identical to the original study's meta-analysis P**. Signed `BETA/SE` is appropriate as a separate statistical diagnostic only after confirming the effect allele, trait scale, coordinate build and N convention. Do not claim novel genetic signals from the transformation.

Upstream LDSC documentation:
- [Required fields for genetic correlation](https://github.com/bulik/ldsc/wiki/What-Data-Are-Necessary-to-Estimate-Genetic-Correlation%3F): SNP ID, alleles, per-SNP N and signed statistic.
- [Original munge_sumstats implementation](https://github.com/bulik/ldsc/blob/master/munge_sumstats.py): original conversion screens duplicated rsIDs, requires P in a usable range, and is sensitive to P=0.
- [Sumstats file schema](https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format): SNP, N, Z, A1, A2 with consistent A1 sign.

## Reproduction

Implementation: `scripts/ckm_stroke/audit_bmi_duplicate_zero_p.py` and tests `tests/test_ckm_stroke_bmi_duplicate_zero_p.py`. Local Orca fixture tests **8/8 PASS**. JSON evidence (locally, gitignored) `artifacts/ckm_stroke/bmi_duplicate_zero_p_20261116.json`.

```powershell
py -3 -m unittest discover -s tests -p test_ckm_stroke_bmi_duplicate_zero_p.py -v
py -3 -m scripts.ckm_stroke.audit_bmi_duplicate_zero_p --out artifacts/ckm_stroke/bmi_duplicate_zero_p_20261116.json
```

The new `.github/workflows/ckm-stroke-bmi-duplicate-pzero.yml` performs offline regression tests on relevant pushes; a deliberate commit containing `[bmi-dup-full]` or manual dispatch re-streams original GWAS twice and archives only a short JSON report.

## Remaining input gates / decision

- **Duplicate resolution: technical PASS** for exact duplicates only; future derived sumstats must assert 9 dropped and unique rsIDs before use.
- **P=0 numerical behavior: quantified**; no unconditional rewrite of original P; decide and preregister a conversion policy that preserves sign and considers source roundoff.
- **HapMap3 and EUR LD reference: technical evidence exists**, but the embedded Zenodo HM3 list differs in size from the MRCIEU pinned reference. Compare rsID sets and allele treatment before assuming interchangeability.
- **EUR reference identity to historical Alkes release and UKB sample overlap:** separate HOLD.
- **LDSC h2, intercept and stroke subtype rg:** **NOT_RUN**; original GIGASTROKE outcome harmonization and SBP/T2D original access still require validation.

**Overall scientific_status = HOLD**, regardless of successful source and unit tests.
