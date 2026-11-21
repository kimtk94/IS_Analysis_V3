# CKM–Stroke BMI: exact variant duplicates and GRCh37 HapMap3 overlap

Audit date: **2026-11-16**. Execution: Orca worktree on NGKR-SAL-TH, verified public GIANT BMI original plus MRCIEU/ldsc GRCh37 HapMap3. These metrics are technical input QC only; **scientific LDSC readiness remains HOLD**.

## Pinned source provenance

- Corrected GIANT/Yengo 2018 EUR BMI original: `https://giant-consortium.web.broadinstitute.org/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz`.
- Full BMI original SHA256 (already checked on Orca and independent GitHub runner): `0d6ed0ea97870916b830ccae349df94ca6f3cc68c025c79e784729af7f7136a4`, compressed size 46,754,937 bytes.
- Reference: MRCIEU/ldsc GRCh37 `w_hm3.snplist` file at repository commit `e9f9739c9c74ad5b4c653a6409a703cc005c62dc`. Raw URL: `https://raw.githubusercontent.com/MRCIEU/ldsc/e9f9739c9c74ad5b4c653a6409a703cc005c62dc/GRCh37/data/hm3/w_hm3.snplist`. SHA256 `4ed5a031f88edff7f94663841160d0fff750c45d4d55a630515b4d6de082c49e`, 17,170,716 bytes.
- MRCIEU reference has **1,217,311 SNP records** plus one header. Initial 1,204,346 unique marker labels included **12,965 repeated marker rows**. The sole internally conflicting label was the literal missing-identifier placeholder **`NA`**, not a real rsID. Excluding `NA` leaves **1,204,345** eligible reference labels. An earlier naive counter reported 8,958 conflicting-allele observations from repeated `NA` lines; these do **not** represent 8,958 independent SNPs.

## Orca whole-source observed results

| Metric | Observed |
| --- | ---: |
| BMI original full rows | 2,336,269 |
| BMI distinct rsIDs | 2,336,260 |
| BMI duplicate rsID rows | 9 |
| Duplicate chr:pos occurrences | 9 |
| Duplicate same chr:pos+unordered allele pair rows | 9 |
| Conservative reference unique rsIDs after 1 conflict exclusion | 1,204,345 |
| BMI–HapMap3 matching rsID rows (includes 2 duplicate rows) | 1,012,697 |
| BMI–HapMap3 **matching unique rsIDs** | **1,012,695** |
| % of conservative HapMap3 reference represented in BMI | **84.0868%** |
| % of BMI unique rsIDs in conservative HapMap3 reference | **43.3468%** |
| Allele pair compatibility for matched rows | 1,012,697 direct unordered-pair matches, no reference complement/mismatch flagged |

**Interpretation:** Strong common-SNP overlap suitable to continue *LDSC preparation*, but **not proof of usable LDSC regression SNP count or of genetic signal**. Nine BMI duplicates require an explicit deterministic deduplication rule (do not take the first record blindly without checking effect-stat agreement). Reference's conflicting rsID must be excluded. Allele pair concordance does not confirm phenotype effect orientation, sample alignment or within-cohort allele frequency. No full-cohort EAF exists here: `Freq_Tested_Allele_in_HRS` is a separate HRS reference frequency.

## Additional hard gates before LDSC

1. Choose a published and version-locked EUR LD-score + regression-weight reference, and independently verify its publisher/hash, 22-chromosome coverage and build. The Broad group historical `w_hm3.snplist.bz2` and `eur_w_ld_chr.tar.bz2` HTTP endpoints returned 404 on 2026-11-16. The MRCIEU GitHub HapMap3 list is accessible, but it alone is **not** the LD-score files or regression weights. A mirrored EUR LD-score set may be usable only after checking release lineage, license, reference sample count and checksum.
2. Inspect the 9 duplicated rsIDs/chr:pos (including P, beta, SE, N) before dropping; record a stable deduplication policy and number retained. Because no variant rows are saved in this audit, a dedicated bounded diagnostic must be run to compare effects if needed.
3. Exclude the 44 original BMI P=0 underflow records from a P-based pipeline or derive Z=BETA/SE under a documented, pre-registered procedure. Do not substitute arbitrary `1e-300` values or claim infinite significance.
4. Evaluate per-variant N, HapMap3 SNP list allele orientation, MHC exclusions, ambiguous palindromes, population LD and outlier alleles when constructing LDSC munged inputs.
5. Investigate participant overlap among BMI GIANT+UKB, GIGASTROKE EUR and additional CKM exposures; use intercept checks, test heritability quality and report uncertainty before stroke subtype `rg`.
6. SBP/EUR and T2D/EUR original release access still requires authenticated authorization/terms. Never bypass.

## Reproduction and CI

`scripts/ckm_stroke/audit_bmi_hm3_overlap.py` streams the original + pinned reference and writes only aggregate JSON to gitignored `artifacts/`. Full-source run:

```powershell
py -3 -m scripts.ckm_stroke.audit_bmi_hm3_overlap --out artifacts/ckm_stroke/bmi_hm3_overlap_20261116.json
```

Offline tests:

```powershell
py -3 -m unittest discover -s tests -p test_ckm_stroke_bmi_hm3_overlap.py -v
```

`.github/workflows/ckm-stroke-bmi-hm3-qc.yml` runs fixture tests on push; it streams the two original inputs **only** on explicit workflow_dispatch or on a commit containing `[hm3-full]`. No full GWAS/reference file is retained or uploaded. CI SUCCESS would demonstrate repeatability, not scientific LDSC success.
