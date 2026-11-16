# CKM × Stroke — Locke 2015 GIANT EUR BMI pre-UKB exposure sensitivity, full raw-source QC

**2026-11-17 KST · Scientific decision: HOLD (not an LDSC GO).** This work adds a historical *pre-UK Biobank* BMI exposure source as a **sensitivity candidate**. It does not replace the previously validated 2018 Yengo GIANT+UKB BMI primary source or claim that the 2015 GIANT cohorts are independent of the GIGASTROKE stroke GWAS.

## 1. First-party source evidence, 2015

- Official publisher [GIANT Consortium 2012–2015 summary statistics](https://giant-consortium.web.broadinstitute.org/GIANT_consortium_data_files), **Locke et al. (2015) BMI EUR**, filename `SNP_gwas_mc_merge_nogc.tbl.uniq.gz` (GWAS Catalog/OpenGWAS cross-index `ieu-a-835`, 2015 EUR, sample-size ceiling ~322,154, GRCh37/hg19).
- Exact original: `https://giant-consortium.web.broadinstitute.org/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz`.
- Exact official [2012–2015 GIANT data-field documentation](https://giant-consortium.web.broadinstitute.org/images/b/b9/Data_File_Description_2012-2015_GIANT_GWAS_Summary_Statistics.docx) was parsed and independently checked. It explicitly defines allele 1 as the **effect allele** for BETA, aligned to **hg19 + strand**; allele 2 is the other allele. It defines the provided allele frequency as **HapMap CEU reference frequency** (not the sampled European BMI cohort frequency), and N as the number of observations.
- Original gzip fields are **exactly** `SNP A1 A2 Freq1.Hapmap b se p N`. Original p-values reflect genomic-control adjustment, and the primary publisher says BETA/SE/p are rounded to limited significant digits for data disclosure protection; **BETA/SE need not reproduce the published adjusted P exactly**.
- Direct original HTTP byte-range first-64-KB/gzip-header probe: HTTP 206, official publisher, 35,869,809 total bytes; no third-party GWAS replacement.
- **Complete original stream in Orca:** 35,869,809 compressed bytes, computed **SHA-256 `ddca3b30985e3428f4762f8e8aadd385504627b9b825cff1ecdc4a5e2e358b44`**, gzip CRC/full EOF valid. This is a **directly computed complete-source baseline**, not an independent publisher-posted MD5 or SHA checksum; CI must independently re-download and reproduce this pinned hash.

## 2. Original QC anomalies — do not silently fix source files

Actual original has **2,554,637 association rows, 2,554,623 distinct rsIDs**, hence **14 extra rows associated with 13 duplicated rsIDs** despite the word `uniq` in the filename. Several duplicated IDs have reversed alleles and opposite BETA signs, or differing HapMap CEU frequency. All source duplicate rsIDs must be **quarantined in their entirety**, not arbitrarily deduplicated to first occurrence or minimum P.

Additional actual source observations:
- **53,635 unique rsID rows with fractional `N`**; example `rs10000027` `N=72852.4`. Preserve as numerical N rather than cast to int, and do not assert that a rounded N represents an actual number of independent individuals.
- **54,674 unique rsID rows with missing/invalid HapMap CEU reference frequency**, but valid BETA/SE/N. This reference field is **not required for effect-direction harmonization** and must not cause otherwise valid association records to be deleted.
- No variant-level effect was selected based on P-value or inferred by replacing missing values. No raw GWAS association table is persisted.

## 3. Exact overlap to pinned 2018 BMI panel (2018 panel is *only the reference for candidate SNPs*)

Prior 2018 corrected GIANT+UKB conservative EUR BMI panel: **1,007,899 non-palindromic unique rsIDs**, canonical plain-text SHA256 `ed8f713827cd3281bea5f88a9cb42a3eabfecb1ee67c912ca7e86f151e4399bf`. The new comparator does **not** import 2018 effect sizes into the 2015 sensitivity exposure, only its validated SNP/allele pairs for technical QC.

| Primary comparison | Exact result |
| --- | ---: |
| 2018 candidate SNPs | 1,007,899 |
| 2015 original rsIDs present among candidates | 1,007,860 |
| Candidate rsIDs absent from the 2015 original | **39** |
| Candidate IDs with ambiguous duplicate source rows removed | **2** (`rs3780697`, `rs627141`) |
| **Unambiguously harmonizable candidate SNPs** | **1,007,858 (99.9959%)** |
| 2015 A1/A2 orientation exactly same as 2018 BMI A1/A2 | **516,430** |
| 2015 effect allele swapped vs 2018 A1 (signed Z flip needed) | **491,428** |
| Complement-only, mismatching, rejected allele numeric candidates after duplicate exclusion | **0** |
| Published 2015 P ≤ 5×10⁻⁸ among these variants | 799 (informational only) |

Matched valid 2015 per-SNP N: **min 50,044.4; P10 216,666; median 233,687; P90 234,005; max 322,153**. This is a **continuous reported N field**, not a constant or an automatically approved LDSC integer N.

The source has **no explicit chromosome/BP columns**; this audit matches source 2015 rsIDs and allele pairs to an already validated GRCh37 2018 BMI rsID panel. It **does not independently certify every 2015 SNP's chromosome or position** or prove sample independence. Previous 44-point Ensembl GRCh37 orthogonal reference checks were performed on the **2018 panel**, not on new 2015 allele effects.

## 4. Study-design interpretation and gates

- The 2015 Locke study predates the 2018 GIANT+UK Biobank meta-analysis and is a materially useful **UKB-avoidance sensitivity candidate**, but "pre-UKB" **is not** "no overlap with GIGASTROKE". The original GIGASTROKE [2022 corrected Supplementary Table 1](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05492-5/MediaObjects/41586_2022_5492_MOESM2_ESM.xlsx) includes EUR cohorts such as ARIC, CHS, FHS, and UKBiobank. Potential **non-UKB shared cohorts with 2015 GIANT require explicit cohort-by-cohort comparison** rather than assumption.
- Changing BMI exposure source also changes GWAS power, study methods, participant selection, SNP coverage and meta-analysis inputs; it **cannot isolate sample-overlap bias by itself**.
- LDSC 2015 exposure needs per-variant N preservation, allele-oriented signed effects with full conversion/QC, correct HapMap3/LD-score chr:BP mapping and genome-build checks, regression model assumptions, and the GIGASTROKE case/control N_eff uncertainty addressed in `docs/ckm_stroke_gigastroke_cohort_overlap_neff_20261117.md`.
- An rsID/allele matched count is **not** an estimate of genetic correlation, heritability, MR causality, or shared individuals.
- Scientific status **`HOLD`**, `ready_for_ldsc=false`, `h2/rg/MR/coloc=NOT_RUN`. Draft PR #284 stays unmerged.

## 5. Independent Linux CI source re-fetch and exact cross-platform evidence

- [GitHub Actions #43983463946](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43983463946): **SUCCESS**, both offline tests (8/8) and entire-original-source validation; timestamp 2026-11-16 23:57:01 UTC (2026-11-17 08:57:01 KST).
- [Artifact #14283897455](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43983463946/artifacts/14283897455): JSON-only temporary 2018 candidate generation provenance and final 2015 source QC; no compressed original association files or full candidate tables uploaded.
- The exact Linux JSON artifact was copied to Orca and independently compared to the original Windows result: all **19 QC count keys**, all per-variant matched-N distribution statistics, original compressed source bytes **35,869,809**, SHA256 `ddca3b30985e3428f4762f8e8aadd385504627b9b825cff1ecdc4a5e2e358b44`, candidate total and match percentage were **identical**; comparison exited 0. Note this raw-source hash is a **pinned primary-source observation reproduced independently**, not publisher-published source checksum.
- New full-source data is **technically validated** as a pre-UKB-era alternate BMI exposure reference; `scientific_status=HOLD` and `ready_for_ldsc=false` remain unchanged.
- **Critical future sensitivity design constraint:** The 1,007,858 loci are a *QC intersection of 2015 original with the 2018-derived HM3 candidate*. They are **not** a standalone 2015-selected variant universe. For a methodologically independent 2015 sensitivity GWAS, derive its SNP set directly from the 2015 original and independently authenticated GRCh37 reference/LD-score and HM3 sources rather than selecting it conditional on the 2018 BMI panel. Compare 2018 and 2015 only after each independently established QC and harmonization universe is defined.
- Complete independence from GIGASTROKE is not established: 2015 GIANT may share non-UKB recruitment studies with stroke cohorts. Participant overlap remains `POSSIBLE_UNQUANTIFIED`, not zero.

## 6. Reproduction

- Main program: `scripts/ckm_stroke/audit_locke_2015_bmi_overlap.py`. Complete raw gzip SHA256 pinned to the Orca observation and rechecked in CI, 2,554,637 original rows and 14 duplicates pinned, strict hashed 2018 candidate consistency, per-row numerical QC, explicit allele swap/sign categories, SNP-specific N quartiles, aggregate-only output.
- Tests: `tests/test_ckm_stroke_locke_2015_bmi_overlap.py` — 8 offline test cases including non-integer N, missing reference frequency, duplicate-rsID quarantine, effects orientation and source drift.
- CI `.github/workflows/ckm-stroke-locke-2015-bmi-source-audit.yml`: automatic lightweight tests; one real ~35.9 MB official-source stream only on `[locke-2015-full]` commit marker or manual dispatch. No original source archive or candidate gzip uploaded; only evidence JSON.
- Windows Orca local ignored JSON: `artifacts/ckm_stroke/locke_2015_eur_bmi_fullsource_qc_20261117.json`.

```powershell
py -3 -m unittest discover -s tests -p test_ckm_stroke_locke_2015_bmi_overlap.py -v
py -3 -m scripts.ckm_stroke.audit_locke_2015_bmi_overlap --candidate artifacts/ckm_stroke/bmi_3way_conservative.sumstats.gz --provenance artifacts/ckm_stroke/bmi_ldsc_candidate_qc_20261116.json --out artifacts/ckm_stroke/locke_2015_eur_bmi_fullsource_qc_20261117.json
```

**Next scientific decision:** independently reproduce full source counts/SHA in Linux, then design a **sidecar-free per-SNP 2015 signed Z/N candidate** (original 2015 BETA/SE, oriented to original 2015 A1) with missingness/rounding safeguards. Restrict initial use to sensitivity/QC until a formal LD score h² / intercept / sample-overlap plan is approved; never treat the pre-UKB source as fully participant-independent without cohort evidence.
