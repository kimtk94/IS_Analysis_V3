# CKM × Stroke — BMI 1,007,899 rsID coordinates vs EUR LD score; GIGASTROKE EUR outcome schema

**Evidence date:** 2026-11-16 KST · Orca on NGKR-SAL-TH. Scientific analysis remains **HOLD**; no LDSC h² or genetic correlation was run.

## 1. Full primary-source coordinate join (real rows)

Reference inputs:
- GIANT Yengo 2018 corrected EUR BMI original: 2,336,269 rows in 46,754,937-byte gzip, pinned source SHA-256 `0d6ed0ea97870916b830ccae349df94ca6f3cc68c025c79e784729af7f7136a4`.
- Three-way strict rsID/allele panel: MRCIEU GRCh37 HM3 ∩ Zenodo HM3 ∩ canonical 22-chromosome EUR LD score = 1,172,984 unique rsIDs. The BMI subset contains exactly **1,007,899** unique non-palindromic SNPs.
- [Zenodo EUR LD score mirror](https://zenodo.org/records/8182036), 33,357,890 compressed bytes; published MD5 `e2f16343c4cfaa76caa7d0c03d26b489`. Canonical `{chr}.l2.ldscore.gz` and `{chr}.l2.M_5_50` for chromosomes 1–22. The archive's noncanonical `6_old` historical files are **excluded**.
- Candidate output created in a prior stage: `artifacts/ckm_stroke/bmi_3way_conservative.sumstats.gz` (9,660,069 bytes), ignored by Git. Candidate canonical plaintext SHA-256 `ed8f713827cd3281bea5f88a9cb42a3eabfecb1ee67c912ca7e86f151e4399bf`.

Exact Orca observed results:

| Metric | Observed |
| --- | ---: |
| Canonical EUR LD score rows across 22 autosomes | 1,290,028 |
| BMI candidate rsIDs | 1,007,899 |
| Candidate rsIDs matched in canonical EUR LD score | **1,007,899** |
| Exact CHR+BP coordinates concordant to GIANT BMI source | **1,007,899** |
| Candidate with different reference CHR/BP | **0** |
| Candidate absent from canonical EUR LD score | **0** |
| Repeated candidate rsIDs in LD score | **0** |
| Candidate markers within chr6 25–34 Mb MHC interval | **0** |
| Noncanonical tar members ignored | 4 |

The MHC count is **0 in this candidate subset** because the canonical reference SNP coverage used here excludes those SNPs. It does not establish that the original GIANT BMI summary statistics were MHC-free.

The script re-opens the original BMI and the pinned Zenodo compressed archive. Full source SHA-256 (BMI), full compressed Zenodo MD5, 22/22 canonical LD-score/M_5_50 members and the local candidate's canonical plaintext SHA-256 are all checked. No raw GWAS or LDscore archive is saved, extracted, or uploaded.

Technical result `PASS_SNP_ID_AND_GRCH37_COORDINATES_FOR_BMI_CANDIDATE`; **ready_for_ldsc=false**. Coordinate match alone does not prove effect-allele harmonization or absence of sample overlap.

## 2. GIGASTROKE EUR outcome official original schema (4 bounded gzip-header probes)

Catalog manifest `config/ckm_stroke/gigastroke_eur_official.tsv` was independently checked against each original GWAS Catalog publisher YAML and `md5sum.txt`, and **the actual first gzip header line** of each official `buildGRCh37.tsv.gz` was read. Original raw archives were **not** downloaded in full; a published md5sum sidecar is **not** a verified raw-source MD5.

| Trait | GWAS Catalog | Build | Publisher YAML + MD5 sidecar | Bounded source gzip header |
| --- | --- | --- | --- | --- |
| Any ischemic stroke (AIS) | `GCST90104540` | GRCh37 | PASS | PASS |
| Large-artery atherosclerosis (LAA) | `GCST90104542` | GRCh37 | PASS | PASS |
| Cardioembolic stroke (CE) | `GCST90104541` | GRCh37 | PASS | PASS |
| Small-vessel disease (SVD) | `GCST90104543` | GRCh37 | PASS | PASS |

All four original source headers have **exactly these eleven columns**:

```text
chromosome
base_pair_location
effect_allele_frequency
beta
standard_error
p_value
odds_ratio
ci_lower
ci_upper
effect_allele
other_allele
```

**No rsID column** exists in these raw files. Downstream outcome harmonization must join on **GRCh37 chr:pos** using unambiguous unique loci. Then compare `effect_allele` and `other_allele` with BMI A1/A2 (BMI A1 is Tested_Allele), explicitly handle swap/sign flip, reverse complement, palindromic markers, multiallelic loci, and per-outcome missingness. Never match on position alone and claim effect orientation. Report for each AIS/LAA/CE/SVD: positions queried, matched/unmatched, alleles same/swapped/strand ambiguous/mismatch, signed beta and SE QC, case-control N accounting, and excluded variants.

GIGASTROKE `effect_allele_frequency` cannot by itself establish BMI cohort-wide frequency. The `Freq_Tested_Allele_in_HRS` from the GIANT BMI source is an HRS-specific frequency rather than a global sample allele frequency.

Published official [GWAS Catalog summary statistics file documentation](https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics/) distinguishes original from harmonized GWAS files and states that the full sumstats FTP site is updated nightly.

## 3. Code, execution and CI

Source scripts:
- `scripts/ckm_stroke/audit_bmi_ldscore_coordinates.py`
- `scripts/ckm_stroke/probe_gigastroke_eur_header.py`
- `tests/test_ckm_stroke_bmi_ldscore_coordinates.py`
- `tests/test_ckm_stroke_probe_gigastroke_eur_header.py`

Reproduce in Orca (not raw source downloads by default in CI):

```powershell
py -3 -m unittest discover -s tests -p test_ckm_stroke_bmi_ldscore_coordinates.py -v
py -3 -m unittest discover -s tests -p test_ckm_stroke_probe_gigastroke_eur_header.py -v
py -3 -m scripts.ckm_stroke.audit_bmi_ldscore_coordinates --candidate artifacts/ckm_stroke/bmi_3way_conservative.sumstats.gz --provenance artifacts/ckm_stroke/bmi_ldsc_candidate_qc_20261116.json --out artifacts/ckm_stroke/bmi_ldscore_coordinate_qc_20261116.json
py -3 -m scripts.ckm_stroke.probe_gigastroke_eur_header --out artifacts/ckm_stroke/gigastroke_eur_outcome_header_qc_20261116.json
```

GitHub Actions workflow `.github/workflows/ckm-stroke-bmi-coordinates-gigastroke-headers.yml`:
- Runs 11 offline tests automatically when the code changes.
- BMI full original-source coordinate verification: explicit `[bmi-coordinate-full]` commit tag or manual dispatch. It generates the 1M candidate temporarily, validates all candidate positions, deletes the candidate, and uploads only compact JSON.
- GIGASTROKE publisher sidecar+bounded gzip headers: explicit `[gigastroke-header-full]` commit tag or manual dispatch. It uploads only a compact JSON report, no full outcomes.
- Both full jobs remain separate so source access failures do not become false scientific PASS statuses.

## 4. Pending scientific gates

1. Full raw-byte MD5 and full row/schema QC of all four GIGASTROKE outcome files (only header/source metadata checked here).
2. Explicit four-way BMI↔GIGASTROKE coordinate + effect-allele and effect-direction harmonization on complete outcome GWAS. No outcome `h2` / `rg` performed.
3. Confirm original GIGASTROKE sample ancestry and BMI vs GIGASTROKE UK Biobank overlap. Check phenotype N and LDSC intercept sensitivity.
4. Inspect outcome per-variant effect allele frequency, palindromic/multi-allelic and strand handling; add per-subtype missingness and unmatched flowchart evidence.
5. SBP/T2D primary GWAS availability and access/legal approval remain separate HOLD.

**Overall scientific_status=HOLD. Draft PR #284 remains unmerged.**
