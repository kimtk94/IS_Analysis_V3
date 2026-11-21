# CKM × Stroke: four original GIGASTROKE EUR files — full MD5 and BMI allele orientation

**Date:** 2026-11-16. **Scope:** technical GWAS source integrity + chr:pos/effect allele category QC. No LDSC h², rg, MR, causal inference, or outcome effect-size association analysis is asserted. The status is **scientific HOLD** by design, regardless of source integrity outcomes.

## Verified official publisher source contract

The repository manifest `config/ckm_stroke/gigastroke_eur_official.tsv` pins four official GRCh37 original EUR GWAS Catalog sources. The prior independent header probe checked their YAML, original publisher MD5 sidecar and 11-column gzip headers.

| Outcome | Accession | Full publisher gzip HEAD Content-Length | Original compressed MD5 recorded in official manifest |
| --- | --- | ---: | --- |
| AIS (ischemic stroke) | GCST90104540 | 283,482,604 bytes | `1bfe4ae8a5042fb24cb0a562b05b0d2f` |
| LAA (large artery stroke) | GCST90104542 | 230,156,887 bytes | `71ea9fd7c06791a831861053b6a1cfbe` |
| CE (cardioembolic stroke) | GCST90104541 | 263,057,516 bytes | `c252262586a6ea84b985f30dd7ab32be` |
| SVD (small vessel stroke) | GCST90104543 | 237,118,943 bytes | `f41bc5f55b3cab238e305e3e86363f69` |

The total is ~1.014 GB of compressed originals. Do not run automatic repeated full downloads on every commit. Official URLs are derived from pinned `GCST9010454x` catalog accession and never guessed from input. Separate original publisher YAML and `md5sum.txt` are checked again before each raw scan.

**Actual field schema, all four:** `chromosome base_pair_location effect_allele_frequency beta standard_error p_value odds_ratio ci_lower ci_upper effect_allele other_allele` (11 TSV columns). **No rsID column**, **no per-variant N column**; published `sample_total` is not automatically per-variant N or effective case-control N. Reference build: GRCh37.

## Exact candidate and harmonization contract

- Build BMI candidate coordinate map from original corrected GIANT BMI compressed source and previously validated derived local candidate (1,007,899 SNPs), checking original SHA-256 and candidate plaintext SHA-256.
- Preserve BMI A1 as original `Tested_Allele`; BMI A2 is `Other_Allele`. Do not use HRS reference allele frequency as study-wide EAF.
- Join on **unique GRCh37 `chromosome:base_pair_location`** only. No rsID is fabricated for outcome rows.
- Classify matching outcome `effect_allele, other_allele` relative to BMI A1/A2: `same` (outcome signed beta retained), `swapped` (flip signed beta), `complement_same` (retain signed beta), `complement_swapped` (flip signed beta), `allele_incompatible`, `non_ACGT_alleles`, `same_allele_invalid`, `palindromic_ambiguous`, `matched_numeric_invalid`. Strand complement is explicitly counted, not silently applied.
- Multiple outcome rows at the same BMI chr:pos **invalidate that coordinate** and are excluded rather than silently keeping the first or best P value.
- Check valid, finite original outcome BETA, SE (SE > 0), EAF [0,1], P [0,1] for the matched SNP set. Effect orientation examples verify signed BETA/SE calculations, but **do not emit complete harmonized per-SNP effects** until downstream audit/approval.
- Iterate the entire raw compressed outcome through a SHA-256/MD5 hashing reader without saving its bytes, and require exact publisher MD5, bytes, gzip decompression, source header, rows, official source host and original metadata. A published sidecar alone does not count as MD5 verified.
- Report per outcome: original rows, chromosomes, BMI match vs nonmatch, same/swapped/complement same/swapped, allele mismatch, invalid matching stats, duplicate-position exclusions, final uniquely alignable loci, unmatched BMI loci, and error evidence. Report **only aggregates** to JSON. Do not upload original/candidate GWAS gzip to GitHub.
- Source and allele-classification technical PASS does **not** equal an approved LDSC input. In particular outcome-specific allele-aware signed-Z/N preparation, case-control effective sample size, ancestry distribution, genome-build anomalies, UK Biobank sample overlap, LDSC intercept and heterogeneity require separate validation.

## Code and reproducibility

- Source: `scripts/ckm_stroke/audit_gigastroke_bmi_allele_harmonization.py`
- Offline test suite: `tests/test_ckm_stroke_gigastroke_bmi_allele_harmonization.py`. Covers opposite effects, reverse complements, numeric anomalies, duplicate positions, MD5 failures and scientific HOLD.
- CI: `.github/workflows/ckm-stroke-gigastroke-four-full-md5-harmonization.yml`, fixture checks on push; four original files only on deliberate `[gigastroke-full4]` push or manual workflow dispatch. Build candidate **once**, stream each of the 4 original gz sources **once**, then remove candidate and publish only two short JSON reports (30-day retention).
- Local Orca example for one subtype:

```powershell
py -3 -m unittest discover -s tests -p test_ckm_stroke_gigastroke_bmi_allele_harmonization.py -v
py -3 -m scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization --phenotype LAA --candidate artifacts/ckm_stroke/bmi_3way_conservative.sumstats.gz --provenance artifacts/ckm_stroke/bmi_ldsc_candidate_qc_20261116.json --out artifacts/ckm_stroke/liaison_LAA_full_md5_alleles.json
```

Future cross-platform comparison: both Windows Orca and Linux CI must report the **same** subtype-specific source MD5, original row count and allele orientation/exclusion histogram. When any source fails, keep the status `HOLD_ORIGINAL_SOURCE_OR_HARMONIZATION`, distinguish from a valid zero match, and do not present scientific PASS.

## Pending research decisions

1. Independent full original data integrity and allele/orientation counts for **AIS, LAA, CE, SVD** (see CI evidence; only mark PASS after actual logs).
2. Evaluate effect allele-frequency vs BMI HRS frequency noncomparability, palindromic exclusions, and eventual per-locus phase/strand.
3. Translate source `sample_total` into justified cohort-specific effective N for case-control LDSC; verify cohort overlap with BMI/UKB.
4. Audit outcome effect scale (case-control log odds), per-variant missingness and formal LDSC h²/intercept; only then propose pairwise genetic correlation.
5. SBP/T2D authenticated original sources remain independent HOLD and must not be bypassed.

**Scientific status: HOLD; `ready_for_ldsc=false`. PR #284 Draft remains unmerged.**


## 2026-11-21 pre-munge per-chromosome source intersection (distinct from LDSC regression SNP count)

New code adds `chr_attrition_pre_munge` to the original four-GWAS MD5-checked streaming audit, partitioning the independently pinned 1,007,899 BMI 2018/HM3/EUR LD-score candidate coordinates for each chromosome 1–22. Every chromosome records candidate positions, unique outcome positions, invalid duplicate positions, absent positions, valid aligned/numeric positions and incompatible/invalid positions. The sum must equal the original whole-genome counts, and the identical candidate chromosome profile must be observed in all four outcome comparisons.

The new independent summarizer `scripts/ckm_stroke/audit_ldsc_chr_attrition.py` verifies the existing whole-genome original publisher MD5 and allele-count reviewer first; it never trusts arbitrary counts from an unsigned JSON. Its own tests and the per-coordinate partition tests reject altered chr indices, missing chromosomes, incorrect subtotals, invalid source MD5 and subtype-specific panel mismatch.

**Important:** pre-munge overlap is *not* the SNP set used in LDSC's h² regression. These counts are upper-bound candidates conditional on the Yengo 2018-derived panel, not independent whole-genome outcome HM3-universe counts. No unapproved case-control N, LD-score weight merge, effect Z threshold, large-chi-square exclusion, MAF or GC handling is inferred. Final regression SNP matching and h² remain **NOT_RUN/HOLD**. The official LDSC documentation warns about analyses with fewer than 200,000 merged regression SNPs and favors 600,000+ for stable precision (not a substitute for input approval).


## Verified full original-input source audit: 2026-11-21 completed run

The original full-source CI [#44853148970](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44853148970) finished **SUCCESS** at 2026-11-21 02:10:47 UTC, Git commit [9465210c](https://github.com/kimtk94/IS_Analysis_V3/commit/9465210cd6bb7638b858e8ac9c83464b5b6181c7). It independently rebuilt the pinned Yengo 2018 conservative source panel from 2,336,269 GIANT original rows; independently streamed all four original EUR stroke GWAS gzip files with complete publisher MD5 checks; passed four-source global invariants and the new 22-chromosome partition invariants. Unit tests: **12 + 8 + 7 = 27/27 PASS**.

Source evidence artifact (four aggregate-only JSON files): [#14686313611](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44853148970/artifacts/14686313611). Its content was downloaded and independently summed by chromosome, matching the existing whole-genome counts for all four outcomes.

### Actual original-source per-chromosome pre-munge overlap (not h² regression SNPs)

The counts below are the original GWAS SNPs uniquely aligned by chr:position and allele with the **preselected BMI 2018 EUR/HM3/LD-score intersection**. They are **not** the final LDSC weight-LD merged SNP counts.

| Chr | BMI × HM3/EUR LD-score 후보 | AIS | LAA | CE | SVD |
|---:|---:|---:|---:|---:|---:|
| 1 | 80,163 | 80,070 | 75,734 | 78,986 | 76,497 |
| 2 | 86,907 | 86,812 | 82,614 | 85,750 | 83,464 |
| 3 | 70,820 | 70,757 | 67,499 | 70,141 | 68,206 |
| 4 | 63,805 | 63,745 | 60,544 | 63,132 | 61,303 |
| 5 | 64,621 | 64,578 | 61,531 | 64,011 | 62,402 |
| 6 | 62,832 | 62,794 | 59,964 | 62,285 | 60,670 |
| 7 | 56,037 | 55,997 | 53,394 | 55,528 | 53,957 |
| 8 | 56,550 | 56,512 | 53,855 | 56,009 | 54,562 |
| 9 | 47,013 | 46,960 | 44,588 | 46,416 | 44,959 |
| 10 | 55,000 | 54,962 | 52,605 | 54,467 | 52,845 |
| 11 | 52,284 | 52,222 | 49,708 | 51,732 | 50,212 |
| 12 | 49,529 | 49,471 | 46,977 | 48,934 | 47,495 |
| 13 | 39,098 | 39,061 | 37,081 | 38,676 | 37,591 |
| 14 | 33,835 | 33,778 | 32,033 | 33,427 | 32,405 |
| 15 | 30,439 | 30,414 | 28,668 | 30,062 | 28,966 |
| 16 | 30,409 | 30,384 | 28,656 | 30,018 | 29,030 |
| 17 | 26,049 | 26,027 | 24,552 | 25,727 | 24,784 |
| 18 | 30,619 | 30,597 | 28,900 | 30,306 | 29,369 |
| 19 | 16,906 | 16,882 | 15,792 | 16,668 | 15,967 |
| 20 | 26,375 | 26,346 | 25,010 | 26,085 | 25,313 |
| 21 | 14,314 | 14,286 | 13,689 | 14,165 | 13,762 |
| 22 | 14,294 | 14,277 | 13,535 | 14,137 | 13,597 |
| **합계** | **1,007,899** | **1,006,932** | **956,929** | **996,662** | **967,356** |

Detailed exclusions from the actual 22-chromosome report:

| Outcome | Candidate absent from outcome | Outcome duplicated candidate positions | Aligned candidate numeric/allele invalid | Pre-munge aligned share |
|---|---:|---:|---:|---:|
| AIS | 967 | 0 | 0 | 99.90406% |
| LAA | 50,970 | 0 | 0 | 94.94295% |
| CE | 11,237 | 0 | 0 | 98.88511% |
| SVD | 40,543 | 0 | 0 | 95.97747% |

**Newly reviewed allele-orientation warning:** the full-source aggregate JSON classifies **all matched SNPs for every subtype as 'same'** (AIS 1,006,932, LAA 956,929, CE 996,662, SVD 967,356), with **zero** swapped, complement-same or complement-swapped calls. The existing independent review sets `uniform_allele_orientation_needs_independent_review=true`. This is technically internally consistent with the two aligned sources; it is **not independent proof** that underlying GWAS effect signs are correctly harmonized. Check the original GIANT tested-allele convention and EBI harmonized outcome effect-allele convention against an independent variant/reference orientation source before any signed rg/MR/coloc or LDSC cross-trait interpretation. Do not silently treat uniformly aligned status as evidence of strand/exposure correctness.

**Strict scope and scientific conclusion:** This is a verified **pre-munge input intersection** only. It does **not** establish post-munging SNP counts, the precise LD-score-weighted regression SNP universe, outcome case/control meta-analysis effective N, original GC correction equivalence, nor h² accuracy. All four outcomes retain `scientific_status=HOLD`, `ready_for_ldsc=false`; **h², intercept, rg, MR and coloc were not run**.

### Next independently actionable gates

1. Compare GIANT tested-allele and EBI/GWAS Catalog effect-allele definitions at a small independently referenced, allele-strand-anchored sample. Separate plausible shared allele convention from a sign inversion or strand error.
2. For each subtype, inspect post-munge proxy SNP loss and LD-score-weight intersections **without inventing N_eff**; report bounds until the data model permits valid munging.
3. Resolve case/control meta-analysis effective-N assumptions, per-variant cohort participation, and public GWAS genomic-control provenance. Only then approve one-trait h² QC, followed by cross-trait rg.


## 2026-11-21 independent effect-allele order study challenge

The 2026-11-21 original GIGASTROKE four-file MD5 check found all shared Yengo-2018 candidate SNPs in the source same orientation. Uniform alignment to a harmonized reference does not prove biological effect allele orientation. An additional source audit independently recreates original Locke-2015 BMI SNP/A1/A2 using only Locke and pinned HM3/LD references, without using Yengo selection; retains 11,506 fractional N observations and verifies the exact independent candidate SHA. It recreates the separately pinned Yengo-2018 candidate and chromosome coordinates, then compares original source A1/A2 across shared rsIDs and four EBI 128KiB GWAS prefixes (strictly validated beta ~ log OR). This bounded sample is not a genome-wide independent effect direction certification, does not retest source whole-file MD5 in this job and must never be called sample or participant independence. Only aggregate orientation counts are emitted. Unit tests prove reversed, complement and nonmatching allele handling; no h2/rg/N_eff approved. Code: scripts/ckm_stroke/audit_locke_yengo_gigastroke_allele_order.py; CI: ckm-stroke-ldsc-singletrait-preregression.yml; new evidence artifact: gigastroke_independent_2015_allele_order_ci.json.


## Independent 2015 vs 2018 vs GIGASTROKE original-data result — 2026-11-21 verified

**Exact CI:** [GitHub Actions #44855617332](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44855617332) SUCCESS; new full-source code commit b189dd5205ba44970b1c88e56a280e755c9b76f6. Combined offline contract 28/28 PASS. Small JSON result artifact [#14687228228](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44855617332/artifacts/14687228228) (includes gigastroke_independent_2015_allele_order_ci.json).

* **Locke 2015:** SHA-256 verified independent 1,070,753 rows; 11,506 preserve fractional N; no Yengo SNP selection, original Locke source completely reverified.
* **Yengo 2018:** SHA-256 verified 1,007,899 candidates from independently rebuilt GIANT EUR source with pinned EUR LDscore/HM3 references.
* **Two independent BMI source effect A1/A2 orders:** exactly 1,007,858 rsIDs in both; 516,430 have the same A1/A2 order; **491,428 have swapped order (48.7596% of shared SNPs)**; 41 Yengo candidates absent from Locke independent candidate. Complement categories: zero for these direct pair-matched nonpalindromic candidates. Thus source-study allele order is demonstrably *not* invariant, even with independently matching alleles.
* **Four EBI source gzip prefixes:** 128KiB per trait, original GRCh37 association header, 13,201 *complete* rows and beta approximately log(OR) validated. **Only 1,858 rows** coincide with both independently constructed BMI panels, so source-prefix summary is not a random genome-wide sample or whole-genome effect-direction certification.

| GIGASTROKE EUR trait | Complete source prefix rows | Three-study shared prefix rsIDs | Source effect A1 vs Yengo A1 same | Source effect A1 vs Locke A1 swapped | Locke reversal proportion within three-way subset |
| --- | ---: | ---: | ---: | ---: | ---: |
| AIS | 3,413 | 424 | 424 | 190 | 44.81% |
| LAA | 3,256 | 493 | 493 | 226 | 45.84% |
| CE | 3,273 | 455 | 455 | 205 | 45.05% |
| SVD | 3,259 | 486 | 486 | 226 | 46.50% |
| **Total** | **13,201** | **1,858** | **1,858** | **847** | **45.59%** |

**Scientific interpretation:** The full original GIGASTROKE vs Yengo candidate audit's 100% same A1/A2 is not a general property of independently published BMI effect alleles. The original Locke A1 orientation differs at nearly half of jointly available BMI loci, and the direction difference recurs in sampled original GIGASTROKE associations. Correct cross-trait analyses MUST explicitly choose, document and verify a single study effect-allele orientation and sign-flip signed Z/BETA when necessary. Source A1 matching alone cannot certify causal correctness, whole-genome consistency or sample independence.

**Limitations / HOLD:** Whole-file MD5 for all four GIGASTROKE outcomes passed separately in #44853148970; this bounded experiment does NOT repeat those full source hashes, nor analyze all original stroke variants. The first 128KiB of each file are not claimed to be an unbiased sample. Full SNP-by-SNP allele harmonization, effect-frequency QC, GC-correction equivalence, post-munge LDSC SNP merge, outcome meta-analysis N_eff and single-trait h²/rg remain unapproved and not run. PR #284 stays Draft.


## 2026-11-21 full-GWAS independent Locke direction challenge — CI execution contract

This update extends the existing four-original-GIGASTROKE-MD5 streaming action. It rebuilds the fully SHA-pinned Locke 2015 independent EUR BMI SNP panel (1,070,753 SNPs) as well as the Yengo 2018 BMI panel (1,007,899 SNPs) and verifies 1,007,858 shared independent candidates including 491,428 reversed source effect-allele orders. Without downloading the original four ~1.014GB total stroke GWASs a second time, the existing single-pass scanner now classifies each genome-wide matched original coordinate against Locke 2015, preserving duplicate-coordinate exclusions. Full 22-chromosome category subtotals and independent original MD5/previous expected counts must agree.

Added independent full-source summary reviewer and offline tests: scripts/ckm_stroke/audit_gigastroke_independent_full_direction.py and tests/test_ckm_stroke_gigastroke_independent_full_direction.py. The action publishes only small aggregate provenance JSON, never signed per-SNP effects. Study-specific A1 orientation comparisons are NOT evidence of causal effect sign agreement. Across-study allele frequency harmonization is NOT established by GIANT HRS or pooled GWAS EAF. No N_eff is approved, no LDSC h2 or rg run, and scientific status remains HOLD. Source-dependent numerical results require an actual successful full-original GitHub Actions run and review.


## Whole-original Locke 2015 direction audit — verified 2026-11-21

**Linux CI:** [GitHub Actions #44856750909](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44856750909) **SUCCESS**, commit [7d498e58](https://github.com/kimtk94/IS_Analysis_V3/commit/7d498e5855d3f29350ce72d0cf7843f4f9842061), original-source evidence artifact [#14688324016](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44856750909/artifacts/14688324016). Offline fail-closed tests **36/36 PASS**. Locke 2015 (1,070,753) and Yengo 2018 (1,007,899) panels independently recreated using pinned source and EUR HM3/LD reference checksums. All four GIGASTROKE original EUR GWAS gzip files passed publisher MD5, original source row counts, allele direction review and 22-chromosome accounting. Aggregate technical status: PASS_FULL_ORIGINAL_MD5_LOCKE_2015_ALLELE_ORDER. Scientific status: HOLD.

| Stroke subtype | Full-original Yengo-aligned SNPs | Also selected by independent Locke 2015 | Same effect A1 order vs Locke | Reversed effect A1 order vs Locke | Reversed % among Locke-shared | Yengo matched but absent in Locke |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| AIS | 1,006,932 | 1,006,893 | 515,926 | **490,967** | **48.7606%** | 39 |
| LAA | 956,929 | 956,908 | 490,081 | **466,827** | **48.7849%** | 21 |
| CE | 996,662 | 996,635 | 510,535 | **486,100** | **48.7741%** | 27 |
| SVD | 967,356 | 967,335 | 495,463 | **471,872** | **48.7806%** | 21 |
| **Total phenotype-specific entries; SNPs repeated between phenotypes** | **3,927,879** | **3,927,771** | **2,012,005** | **1,915,766** | **48.7749%** | **108** |

**Interpretation:** Opposite *effect-allele coding order* is not evidence of opposite biological effect. Yengo 2018 and original stroke alleles were same-oriented on their full intersecting loci. Locke 2015 is from a separate original source and has reversed A1 order at ~49% of Locke/Yengo/stroke jointly present variants. Any downstream signed BETA/Z alignment MUST transform study-specific effect signs to a chosen reference allele when needed. The four stroke original gzip archives were independently streamed, publisher full-file MD5 checked, and only compact allele-count aggregates were uploaded.

**Limitations:** These are pre-munge candidate intersections conditional on the Yengo 2018/HM3/EUR-LDscore panel, **not actual LDSC regression SNP counts**. Exact individual-level source participant overlap and allele-frequency equivalence between GIANT HRS and GIGASTROKE cohorts are unknown. GIGASTROKE meta-analysis per-variant N_eff, per-SNP effect scale/GC status across the complete source and final LDSC regression merge remain unapproved. h2/intercept/rg/MR/coloc NOT_RUN, scientific HOLD; Draft PR #284 unmerged.

## 2026-11-21: Original EAF / HRS reference / P-vs-Wald full-source preregression gate

The established original EUR GIGASTROKE four-file streaming pipeline is extended to check **only originally matched, valid, unique source coordinates** after duplicate quarantine. A SHA-pinned GIANT 2018 original column named Freq_Tested_Allele_in_HRS is read within the existing source-verified BMI pass and matched by GRCh37 coordinate to the original outcome; its frequency is a **reference HRS population** and not the GIGASTROKE original outcome study cohort EAF. Following A1/A2 orientation (flip EAF to 1-EAF on swapped direction only), report absolute frequency-difference histogram (2%, 5%, 10%, 20% bands). There is no scientific requirement that the HRS and GIGASTROKE cohorts have matching EAFs, so large shifts are **descriptive / REVIEW**, not automatic ancestry-mismatch conclusions.

The original GIGASTROKE fields BETA, SE and P define a descriptive per-SNP P-vs-Wald input consistency check. Compare original -log10(P) with the two-sided normal-Wald -log10(P) implied by BETA/SE, with stable log-tail handling. Count zero/underflowed original P separately (no arbitrary replacement). Also summarize source EAF MAF distribution, EAF exact 0/1 boundaries, and BETA versus ln(OR). **GC or cov-LDSC preprocessing cannot be inferred from these differences alone**; only source provenance and original supplementary methods can resolve it. Do not silently rescale P, SE or beta.

Source code: scripts/ckm_stroke/audit_gigastroke_eaf_gc.py and extensions to scripts/ckm_stroke/audit_gigastroke_bmi_allele_harmonization.py. Fail-closed fixture checks: tests/test_ckm_stroke_gigastroke_eaf_gc.py. GitHub Actions pipeline: .github/workflows/ckm-stroke-gigastroke-four-full-md5-harmonization.yml with [gigastroke-full4] opt-in. The old pinned original MD5, full allele direction and chromosome ledgers are rechecked. Only bounded aggregate JSON is uploaded, never per-SNP effects/frequencies.

**Approval boundary:** A successful descriptive audit is not an approved N_eff, not a validated GC correction coefficient, not a final post-munge LDSC regression SNP count and not a completed h2 or rg. Keep scientific HOLD. Original LDSC munge behavior: https://github.com/bulik/ldsc/blob/master/munge_sumstats.py .


## GC provenance and original GIGASTROKE methods (checked 2026-11-21)

The original Nature 2022 GIGASTROKE study (Mishra et al., https://doi.org/10.1038/s41586-022-05165-3) explicitly reports using covariate-adjusted LDSC (cov-LDSC) on ancestry-specific GWAS meta-analyses WITHOUT genomic-control correction to test for genomic inflation and estimate SNP heritability. This published analysis description does NOT prove that the corresponding public GWAS Catalog EBI file has identical normalization to the internal discovery input, nor does it document per-SNP effect or meta-analysis effective N.

The new full-source EAF/P-Wald diagnostic distinguishes source EAF and the independent GIANT original column Freq_Tested_Allele_in_HRS. HRS frequency is a cohort-specific descriptive reference, NOT the corresponding original GIGASTROKE cohort EAF. Only compare after verified effect allele coding orientation; observed frequency discrepancies do NOT establish genotype strand errors or ancestry admixture without other supporting evidence.

Reported P versus beta/SE-derived Wald P differences are descriptive. Separate underflowed original P=0 and DO NOT impute arbitrary numerical P. Original P, BETA and SE are never silently modified by an inferred genomic-control correction or an estimated inflation coefficient; marked divergences are REVIEW, not a final GC diagnostic.

Outstanding gates remain the original EBI publicly provided effect scaling / GC provenance, valid SNP-specific case-control meta-analysis effective N, ancestry-specific sample overlap, independently validated LDSC regression-weight reference SNP merge, final post-munge counts and h2/rg. Science remains HOLD, no LDSC/h2/rg/MR performed by this QC workflow.


## 2026-11-21 original GIGASTROKE EAF / P-Wald full-source results and secondary audit

**Verified original-source workflow**: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44860586028 SUCCESS; immutable source commit 7c1ca53a57f3e44b1cdb8e2cf19753b943cf45ba. Compact full original MD5 + EAF/GC evidence artifact ID 14689563883. EAF/Wald technical status PASS_ORIGINAL_MD5_FOUR_EAF_AND_WALD_DESCRIPTIVE_ONLY. Scientific status HOLD, no LDSC or downstream causal interpretation. No per-SNP GWAS effects released.

| Phenotype | Matched pre-munge SNPs | Original MAF below 0.05 | abs EAF vs GIANT HRS reference >= 0.05 | HRS >= 0.10 | HRS >= 0.20 | abs -log10(P) vs beta/SE Wald delta 0.1 to 1 | Delta >= 1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| AIS | 1,006,932 | 71,858 | 1,835 | 95 | 1 | 14 | 0 |
| LAA | 956,929 | 29,737 | 301 | 86 | 1 | 0 | 0 |
| CE | 996,662 | 64,994 | 332 | 88 | 1 | 14 | 0 |
| SVD | 967,356 | 39,757 | 308 | 88 | 1 | 0 | 0 |

**Primary review targets**: each outcome has *one* >0.20 absolute original EAF versus GIANT HRS reference frequency discrepancy on **chromosome 11**. This is repeated CHROMOSOME identity, not proof of the same SNP. P-vs-Wald discrepancies in the 0.1–1 log10 band: AIS chr1=1, chr4=9, chr9=1, chr12=3; CE chr4=14; none in LAA/SVD. No phenotype has a P-vs-Wald discrepancy >=1 log10. Original BETA = log(OR) to strict 1e-10 accuracy for all valid matched SNPs.

**Follow-up workflow**: .github/workflows/ckm-gigastroke-eaf-wald-artifact-review.yml, triggered by changes to its own source/tests or via workflow_dispatch. It verifies that the input source run completed successfully on its exact pinned commit, downloads only its trusted preexisting compact JSON artifact and reruns the independent publisher-MD5/source-ledger reviewer; no extra original GWAS download or SNP-level effect upload. Publishes JSON and Markdown chromosome-level anomaly summaries. Code scripts/ckm_stroke/review_gigastroke_eaf_wald_anomalies.py; tests/test_ckm_stroke_gigastroke_eaf_wald_anomalies.py.

**Scientific restriction**: GIANT HRS is a separate frequency reference, not GIGASTROKE outcome cohort EAF. Frequency differences cannot prove strand errors or ancestry mismatch by themselves. P-Wald differences are not evidence of a GC correction factor. Aggregate source matches are *not* final post-munge LDSC regression SNP counts; outcome per-SNP N_eff remains unknown, and h2, intercept, rg, MR and coloc are NOT RUN. Do not modify original effect sizes/P-values based on this descriptive audit.


## 2026-11-21 targeted original source chr11 SNP identity and chr4 P/Wald rounding QC

**Source evidence:** [GitHub Actions #44863365730](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44863365730) SUCCESS, commit [4860eff034949d82802aeeb0e9c32294f03d0669](https://github.com/kimtk94/IS_Analysis_V3/commit/4860eff034949d82802aeeb0e9c32294f03d0669), [artifact #14690824408](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44863365730/artifacts/14690824408). Seven offline tests PASS. All four original GIGASTROKE EUR public gzip files were independently streamed and original publisher full-file MD5, full source SHA256 and previously proven chr4/chr11 EAF/Wald histogram counts reproduced. The output retains only selected source SNP identities and numeric-precision diagnostics, no per-variant beta or SE.

### One shared chr11 SNP, not merely shared chromosome

The sole >=0.20 GIANT HRS-reference frequency departure in all four source outcomes maps to **rs198472**, GRCh37 **chr11:61,527,266**, allele set **A/G**, with the candidate rsID mapped by pinned GIANT Yengo2018 original source. All four represent the same coordinate, rsID and compatible allele set.

| Stroke outcome | GIGASTROKE EUR EAF reoriented to GIANT Tested_Allele | GIANT HRS reference frequency | Absolute difference |
| --- | ---: | ---: | ---: |
| AIS | 0.7749 | 0.5341 | 0.2408 |
| LAA | 0.7534 | 0.5341 | 0.2193 |
| CE | 0.7575 | 0.5341 | 0.2234 |
| SVD | 0.7588 | 0.5341 | 0.2247 |

**Interpretation restriction:** the GIANT HRS column is a frequency reference from a different population/cohort, not the GIGASTROKE cohort-specific EAF; therefore a 20pp difference does not alone establish a strand error, sample mix-up or ancestry failure. The locus lies near the chr11q12.2 FADS1/FADS2 LC-PUFA/lipid biology region; the published hg19 interval is approximately chr11:61,540,615–61,664,170 (https://pmc.ncbi.nlm.nih.gov/articles/PMC5874031/). rs198472 is approximately 13kb proximal to this interval, but physical proximity is not proof of a particular causal gene.

### chr4 moderate P-versus-Wald discrepancy

The source EAF/Wald full-scan classified AIS chr4 **9** variants and CE chr4 **14** variants with absolute discrepancy 0.1–1 in -log10(P) between original reported P and two-sided Wald P implied by source BETA/SE. On direct original-source extraction, **all 23** are compatible with rounding intervals of the originally printed BETA/SE/P decimal strings. This does NOT demonstrate that rounding was the actual explanation, and it is NOT a genomic-control correction estimate. No original source effects were modified.

Exactly **three** variants occur in both AIS and CE, matching coordinate + rsID: rs10032150 (chr4:111,646,618); rs17042088 (chr4:111,654,814); rs17042171 (chr4:111,708,287). All lie in the 4q25 region near PITX2, a known atrial-fibrillation and cardioembolic-stroke risk locus (https://www.nature.com/articles/nature06007; https://pmc.ncbi.nlm.nih.gov/articles/PMC5386395/). CE additionally has rs2200733 (chr4:111,710,169), a documented AF risk SNP. AIS affected SNPs were grouped near chr4:111.65–111.72 Mb (5) and chr4:155.44–155.54 Mb (4); the latter lies near the fibrinogen region, with GRCh37 FGB chr4:155,484,132–155,493,959 and FGA chr4:155,504,278–155,511,894 (ClinGen).

**Still HOLD:** exact original GIGASTROKE case/control meta-analysis SNP-specific effective N, sample overlap and final post-munge LDSC regression-weight reference merge are not verified. No h2/intercept/rg/MR/coloc was calculated. This is technical source-QC validation only. Draft PR #284 remains unmerged.


## 2026-11-21 next gate — rs198472 GRCh37 1000 Genomes EUR allele frequency versus original GIGASTROKE

The MD5 and SHA256-verified four-outcome chr11 outlier rs198472 (GRCh37 chr11:61,527,266, A/G) has reported source GIGASTROKE EAF reoriented to the GIANT 2018 BMI tested allele of AIS=0.7749, LAA=0.7534, CE=0.7575 and SVD=0.7588. GIANT 2018 original HRS tested-allele reference frequency = 0.5341. This discrepancy by itself DOES NOT establish cohort or strand error.

To obtain an INDEPENDENT European population reference, fetch rs198472 via the official Ensembl GRCh37 REST variation endpoint with population allele frequencies for 1000GENOMES:phase_3:EUR. The exact original GIANT Yengo2018 tested allele is recovered from the source gzip by streaming to the published SHA256 completion sentinel. Never assume a sign or a frequency from the sorted A/G allele set. Require the Ensembl variant ID, GRCh37 chromosome coordinate, alleles, and a biallelic EUR frequency partition summing to one. No EUR estimate is accepted until these checks pass. A reliable request may be unavailable due to network/API issues: report technical HOLD if so rather than inventing allele frequency.

**Methods sources:** https://www.internationalgenome.org/allele-frequency-calculator-documentation/ , https://grch37.rest.ensembl.org/documentation/info/variation_id , https://www.ensembl.org/info/genome/variation/species/populations.html .

**Case/control and N_eff evidence:** The published source reports AIS 62,100 cases / 1,234,808 controls, LAA 6,399 / 1,234,808, CE 10,804 / 1,234,808, and SVD 6,811 / 1,234,808. Illustrative 4*Ncase*Ncontrol/(Ncase+Ncontrol) is computed for transparent arithmetic only and must NEVER be inserted as original per-SNP effective N in a multi-study inverse-variance meta-analysis. Public original EUR files have exactly eleven columns and no per-variant case count, control count or N_eff, whereas the meta-analysis combines many heterogeneous cohorts and variant availability differs across contributing studies. Total N is not equivalent to effective N for highly imbalanced binary outcomes. Exact contributor N_eff, GC transformation provenance and LDSC regression-weight SNP merge remain HOLD. Reference: https://www.nature.com/articles/s41586-022-05165-3 .

Execution code: scripts/ckm_stroke/verify_rs198472_independent_EUR_and_neff.py, tests/test_ckm_stroke_rs198472_independent_EUR_neff.py and .github/workflows/ckm-rs198472-1000g-eur-neff.yml. Never publish raw per-variant effects, and never merge Draft PR #284 without approval.


## 2026-11-21 rs198472 independent 1000 Genomes Phase 3 EUR original-source resolution — VERIFIED

**Original-independent CI:** [GitHub Actions #44865055885](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44865055885), SUCCESS, 9/9 fail-closed offline tests PASS; [immutable evidence artifact #14691656333](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44865055885/artifacts/14691656333). Research code commit 6f897a2df785e19e015aad0da6b64eee7941236d. The original GIGASTROKE four-subtype chr11 identity source was pinned to original verified workflow #44863365730. The Yengo 2018 GIANT original gzip was independently streamed with its previously pinned source SHA256, confirming rs198472, GRCh37 chr11:61,527,266, tested allele **A**, other allele G and HRS tested-A frequency 0.5341. Ensembl **GRCh37** original variation REST returned rs198472 A/G mapping and independent `1000GENOMES:phase_3:EUR` A frequency **0.7554672** (G frequency **0.2445328**), and passed a strict EUR biallelic sum-to-one gate. Neither GRCh38 nor a different effect allele was used.

| Cohort reference / phenotype | Allele A frequency | Absolute deviation vs 1000G EUR A |
| --- | ---: | ---: |
| 1000 Genomes Phase 3 EUR | **0.7554672** | 0 |
| GIANT Yengo 2018 HRS | 0.5341 | **0.2213672** |
| GIGASTROKE AIS | 0.7749 | **0.0194328** |
| GIGASTROKE LAA | 0.7534 | **0.0020672** |
| GIGASTROKE CE | 0.7575 | **0.0020328** |
| GIGASTROKE SVD | 0.7588 | **0.0033328** |

**Decision:** The previously reported ~20 percentage point difference between GIGASTROKE EUR and the distinct GIANT HRS reference is **NOT reproduced against independent 1000 Genomes Phase 3 EUR**. The first-pass 'large between-study EAF discrepancy' must therefore NOT be classified as a GIGASTROKE strand/allele-coding error. HRS-specific source frequency is an open provenance/QC question, not evidence that the stroke frequency is wrong. 1000G EUR Phase 3 includes its own European reference samples and is not identical to the stroke study cohorts; AIS still differs by ~1.94 percentage points, not a perfect match. No original per-SNP effects or genotypes were exported.

**GIGASTROKE effective N decision:** the source meta-analysis comprises EUR cohorts and shows only 11 columns (no per-variant N/Ncase/Ncontrol). The phenotype total N values and illustrative harmonic 4*Ncase*Ncontrol/(Ncase+Ncontrol) calculations are:
- AIS cases=62,100 controls=1,234,808 total=1,296,908; illustrative N_eff=236,505.83
- LAA cases=6,399 controls=1,234,808 total=1,241,207; illustrative N_eff=25,464.04
- CE cases=10,804 controls=1,234,808 total=1,245,612; illustrative N_eff=42,841.16
- SVD cases=6,811 controls=1,234,808 total=1,241,619; illustrative N_eff=27,094.55

These are population-wide arithmetic illustrations, **NOT effective N per SNP**, and do not account for variant-specific missingness, imputation or per-contributing-study case/control balance. They are NOT approved for LDSC munging or regression. The original Mishra 2022 GIGASTROKE study used cross-cohort meta-analysis and cov-LDSC; exact original phenotype-specific allele coverage and cov-LDSC settings cannot be assumed from the published N totals alone. **Scientific HOLD remains**, `ready_for_ldsc=false`, no h2/intercept/rg/MR/coloc executed.

**Next practical scientific gate:** identify a defensible per-variant N_eff source, or explicitly design a prespecified sensitivity-only constant-N analysis (label unvalidated and separate from primary) after confirming methodology and a source-identical population reference, then perform actual outcome munging and independent 22 chromosome EUR LD-score/weight SNP merge. Prefer a stroke GWAS with documented per-variant N_eff when the original EUR meta-analysis does not expose it. Draft PR #284 remains unmerged.


## 2026-11-21 22-chromosome EUR LDSC weight-reference intersection and cohortwise-N source gate

The following new source-backed test is separated from all h2/intercept/rg/MR execution by a hard fail-closed scientific HOLD: scripts/ckm_stroke/audit_ldsc_weight_reference_merge.py and tests/test_ckm_stroke_ldsc_weight_reference_merge.py, running under .github/workflows/ckm-eur-ldsc-ref-weight-neff-gate.yml. It verifies all 22 chromosome LDSC .l2.ldscore.gz file SNP identifiers, the Zenodo reference w_hm3.snplist, and each chromosome .l2.M_5_50 normalizer, while separately auditing whether optional .l2.M is present from the original pinned archive (Zenodo 8182036, MD5 e2f16343c4cfaa76caa7d0c03d26b489, SHA256 9537f00eb0d163a935aaa2cf04b358b7cf21852279b9c7925802526f6060b069). It independently reconstructs the pinned 1,007,899 GIANT Yengo 2018 pre-munge candidate SNPs and compares chromosome-level rsID coverage with the publisher MD5-verified GIGASTROKE original genome-wide SNP attrition. The exact final **post-LDSC-munge**, weight-regression SNP set remains UNESTABLISHED because outcomes lack a justified SNP-specific N_eff and per-SNP outcome post-munge reference matching was not performed. A shared directory for LDSC ref-LD and weight-LD is **not** an independent reference-weight validation. The reference files are a technical input check only.

Original corrected GIGASTROKE 44-EUR-cohort source supplement (publisher SHA256 7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4), previously verified in workflow #44855617332 (original source SHA of that run b189dd5205ba44970b1c88e56a280e755c9b76f6), provides case/control counts for the meta-analytic cohorts. The source distinguishes a single pooled harmonic 4*Ncase*Ncontrol/Ntotal from the sum of harmonic effective sizes across cohorts with published control totals. These do **not** represent original per-SNP meta-analysis Neff.

| Subtype | Pooled harmonic Neff illustration | Cohortwise known-control sum Neff illustration | Cohortwise / pooled ratio | Original cases lacking cohort control totals |
| --- | ---: | ---: | ---: | ---: |
| AIS | 236,505.833 | 195,570.717 | 82.7% | 103 |
| LAA | 25,464.041 | 23,149.571 | 90.9% | 0 |
| CE | 42,841.160 | 39,605.585 | 92.4% | 28 |
| SVD | 27,094.551 | 23,758.723 | 87.7% | 39 |

Both source examples ignore variable per-SNP participating cohorts and meta-analysis weighting, and some cases lack known cohort-specific control counts. As emphasized by Grotzinger et al. (2023, https://pmc.ncbi.nlm.nih.gov/articles/PMC10066905/), pooled ascertainment can bias case-control GWAS meta-analysis liability-scale h2. Hence do not silently use these illustrations as primary LDSC N; explicit publication/author SNP-level participation evidence would be needed for a primary meta-analysis h2. Alternative well-documented public stroke GWAS case/control SNP-specific N should be assessed before substantive h2/rg. Published 2022 GIGASTROKE cov-LDSC measurements remain external literature, not output of this project.

**Execution status**: the source-backed CI run must report its own 22chr actual results before numerical reference coverage can be promoted. HOLD for outcome SNP-wise N_eff, post-munge SNP set, liability h2, intercept and rg regardless of passing reference QC. The Draft PR #284 stays unmerged.


### Exact original LDSC archive files: mandatory M_5_50 present, optional M absent

During source-backed CI #44865845950 and #44866026394, Zenodo 8182036 full compressed-byte MD5 and SHA256 both matched their pinned public checksum. The tar archive contains **22/22 EUR .l2.ldscore.gz**, **22/22 .l2.M_5_50**, **w_hm3.snplist**, and **0/22 .l2.M** (49 total members). The initial code incorrectly demanded optional .l2.M under a common-SNP normalizer workflow and appropriately stopped instead of silently passing. This has now been corrected explicitly: all 22 default common-SNP M_5_50 files are strictly mandatory, and the absence of all 22 M files is preserved as source provenance. An LDSC --not-M-5-50 workflow is **not supported** with this source archive and must not be silently attempted.

The correction is a file-inventory policy fix, not a change to GWAS association scores, reference variant IDs, primary N_eff assumptions or scientific HOLD. It does not establish LDSC post-munge regression SNP numbers. The source archive is not independent of the reference originally used to preselect Yengo 2018 variants.


### VERIFIED 22-chromosome LD reference completion — 2026-11-21

**[Original-source CI #44866183271](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44866183271) SUCCESS**; independent aggregate artifact [#14692221160](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44866183271/artifacts/14692221160). The source archive was 33,357,890 bytes (published MD5 and full original SHA256 PASS), with 49 archive members. Exactly 22 LD-score gzip chromosomes, 22 M_5_50 normalization files, original w_hm3.snplist and **no optional M** files. The default common-variant M_5_50 mode is technically complete; --not-M-5-50 remains unsupported for this specific published archive.

All **1,007,899** pinned GIANT Yengo2018 pre-munge candidate rsIDs mapped uniquely to one of the 22 canonical EUR LD-score chromosomes and to the source-bundled HM3 panel, and their 22-chromosome counts matched the previous independently publisher-MD5-verified four GIGASTROKE original reference partition. The prior four original outcome intersections remain AIS=1,006,932; LAA=956,929; CE=996,662; SVD=967,356. These are **pre-munge** matches, not an LDSC merge. Prior Yengo candidate construction intentionally used this same LD-score archive as a filter, so these findings provide independent implementation/aggregate integrity verification, not statistical independence from the underlying reference source.

The corrected official original 44-cohort EUR outcome table yields pooled versus known-control cohortwise illustrative harmonic effective-N differences of 17.3083% (AIS), 9.0892% (LAA), 7.5525% (CE) and 12.3118% (SVD). Case counts with missing contributing cohort control totals: AIS 103; LAA 0; CE 28; SVD 39. Both scenarios remain **NOT APPROVED for SNP-specific N_eff**. The input statistical/meta-analysis weighting uncertainty is not solved by any reference-file success. No new LDSC h2, intercept, liability transformation, rg, MR or coloc values have been calculated. Scientific and manuscript readiness remain HOLD.

**Next validation:** compute GIGASTROKE four-outcome *post-munge* variant-to-weight-score intersections only under an explicit, externally justified source N_eff policy; otherwise use a separate stroke GWAS with documented per-variant contributing sample count, or pre-register a clearly labeled constant-N sensitivity run that is never mistaken for the primary analysis. Independently confirm allele-frequency, INFO/MAF and genomic-control handling against source methods.


## 2026-11-21 official alternative stroke GWAS audit — per-variant N_eff availability

Research decision is organized into PRIMARY vs REPLICATION vs SUBTYPE COMPARISON; it is NOT appropriate to call another GWAS an N_eff solution based on its phenotype total N alone.

1. **GIGASTROKE EUR 2022 (primary, four outcomes)**: Original GWAS Catalog 11-column GRCh37 files have no per-SNP effective sample size or per-SNP case/control sample count. Publisher-MD5 four-outcome CI #44860586028 and follow-up EUR reference verification #44866183271 passed, but per-variant N_eff and final LDSC regression SNP set remain HOLD. No primary LDSC h2/rg has run.
2. **MEGASTROKE EUR 2018 (possible external direction/subtype comparison)**: Official https://www.megastroke.org/ lists only MarkerName, Allele1, Allele2, Freq1, Effect, StdErr and P-value. The EUR meta-analysis contains 40,585 stroke cases and 406,111 controls; listed stroke traits are AS, AIS, LAS, CES, SVS. The publication's variant inclusion filter n_cases > 50% is a threshold, NOT a released per-variant n_cases column. The site has a data-use policy and project description/download step. Earlier MEGASTROKE contributing cohorts may overlap GIGASTROKE; the independence of study participants is not established. Not an automatic solution to N_eff.
3. **FinnGen R13 (potential phenotypically defined replication)**: Official https://finngen.gitbook.io/documentation/data-download.md says Release 13 and access through an online form, with manifest https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv. Official https://finngen.gitbook.io/documentation/data-description.md documents GRCh38 variants, 13 GWAS columns (#chrom, pos, ref, alt, rsids, nearest_genes, pval, mlogp, beta, sebeta, af_alt, af_alt_cases, af_alt_controls). The public schema does NOT provide SNP-wise n_cases, n_controls or N_eff. Allele frequencies among cases/controls must never be confused with sample sizes. Stroke endpoints in the manifest must be reviewed for case definition and subtype match, not blindly called ischemic stroke. GRCh38-to-GRCh37 conversion and population differences must be controlled.

Source-verified executable: scripts/ckm_stroke/audit_alternative_stroke_gwas_neff_sources.py; fail-closed tests: tests/test_ckm_stroke_alternative_gwas_neff_sources.py; GitHub Action: .github/workflows/ckm-stroke-alternative-gwas-neff-source-audit.yml. It reads public consortium descriptions and, when technically available, the FinnGen R13 manifest only, and publishes a small metadata audit without downloading controlled GWAS or publishing individual-level data. Manifest access failure is recorded separately; it never manufactures candidate endpoint coverage.

**Primary scientific HOLD resolution requirement:** an outcome with independently verified allele/effect convention, EUR-appropriate LD reference, harmonized build, well-defined ischemic-stroke subtype, and a defensible genomewide effective sample-size model (preferably per variant or based on documented participation and weights). Case/control harmonic N calculated from published total counts and cohort-level incomplete counts are exploratory sensitivity illustrations; cannot substitute automatically for per-SNP meta-analysis N. Even if an alternative GWAS is available for replication, its availability alone does not remove GIGASTROKE's primary N_eff problem. No new h2, intercept, rg, MR or coloc is asserted.


### 2026-11-21 CI-completed alternative EUR GWAS and FinnGen R13 endpoint manifest triage

**[Live official schema/manifest CI #44868507970](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44868507970) SUCCESS** and **10/10 fail-closed tests PASS**. [Artifact #14692489529](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44868507970/artifacts/14692489529) contains the exact source provenance, full official document/manifest SHA256, bounded no-N-evidence report, and candidate codes. External sources checked: official MEGASTROKE public schema, official FinnGen R13 public schema/access instructions, and the official R13 GCS manifest. The previous publisher MD5 GIGASTROKE source audit #44860586028 was revalidated by exact repository, run SUCCESS and source commit hash.

**Official source findings:**

| Data source | Officially described genomewide variant-wise N/Neff field | Outcome use | Gate |
| --- | --- | --- | --- |
| GIGASTROKE EUR 2022 (original four outcomes) | **No** (11 original fields) | PRIMARY MR/coloc associations and LDSC candidate pending Neff | HOLD for new LDSC |
| MEGASTROKE EUR 2018 | **No** (7 described fields: MarkerName/Allele1/Allele2/Freq1/Effect/StdErr/P-value; n_cases>50% is a variant inclusion filter) | EUR AIS/LAS/CES/SVS subtype-direction comparison, note original cohort overlap possibility | HOLD for per-SNP LDSC Neff |
| FinnGen R13 | **No** (13 documented GRCh38 fields; allele frequencies are NOT sample counts) | Candidate Finnish phenotype-defined replication AFTER endpoint definition and release access | HOLD for per-SNP LDSC Neff and genome build harmonization |

**FinnGen R13 manifest actually retrieved and parsed:** original manifest SHA256 c9caafa9b98ee5ef050a705766bd1c451fc2154fb8e81b069022012a2f75ae20, **2,754 endpoints**; **15** initial keyword matches (an intentionally broad scan, NOT 15 ischemic stroke phenotypes). Three cases merit phenotype-definition review: **I9_STR**, **I9_STR_EMBOLIC**, **I9_STR_EMBOLIC_WITH_FLIMMER**. Four identified bleeding codes are excluded from the ischemic stroke candidate set: **I9_ICH**, **I9_SAH**, **I9_STR_SAH**, **I9_SAHANEUR**. Other keyword hits include non-stroke phenotypes such as splenic infarction, benign cerebral meningeal tumors and congenital cerebrovascular malformations; false-positive exclusions are essential.

The actual FinnGen manifest fields in this test did not provide published phenotype-description labels, so NONE of the three candidate endpoints is yet validated as the clinical GIGASTROKE AIS or CE subtype. No LAA- or SVD-specific phenotype was verified. Any FinnGen replication will need the R13 endpoint dictionary/ICD code definitions, sufficient case counts, and explicit cohort-overlap assessment. R13 effect alleles use GRCh38 ALT and require validated build/allele harmonization to original GRCh37 outcomes, not direct positional joining.

**Scientifically decisive negative result:** no alternative source has yet demonstrated a per-SNP case/control/effective N for approved primary LDSC h2/rg. The source checks **do not authorize an arbitrary constant N or pooled 4*Ncase*Ncontrol/Ntotal**. However, MR/coloc pipelines using properly harmonized effect sizes and uncertainties may be assessed separately subject to their own QC; LDSC Neff HOLD does not itself disprove those designs. No new h2/intercept/rg, causal inference or experimental result is asserted. Draft PR #284 remains unmerged.


## 2026-11-21 FinnGen R13 original clinical definitions + released GWAS eligibility — SOURCE VERIFIED

[**GitHub Actions #44869603525**](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44869603525) SUCCESS; source-backed [artifact #14693172602](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44869603525/artifacts/14693172602); 10/10 offline positive/negative tests PASS. Exact official original input sources:

- Public DF13 endpoint case definitions (original 2025-08-14): https://www.finngen.fi/sites/default/files/inline-files/FINNGEN_ENDPOINTS_DF13_Final_2025-08-14_public.xlsx ; SHA256 0004624f727b96b4a56c6fac75ff182eb8bc3daf78cdf690bdad276b31fca6c7.
- Public DF13 control definitions: https://www.finngen.fi/sites/default/files/inline-files/Controls_FINNGEN_ENDPOINTS_DF13_Final_2025-08-14_public.xlsx ; SHA256 7cbbf9072296e34f6bde44b2958c7520c7caedeade924d941e5bbc4fa89cc29b.
- Official R13 public GWAS manifest (2,754 GWAS phenotypes): https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv ; SHA256 c9caafa9b98ee5ef050a705766bd1c451fc2154fb8e81b069022012a2f75ae20. Actual headers: phenocode, phenotype, category, num_cases, num_controls, path_bucket, path_https. Endpoint case/control counts are **release-level phenotype sample metadata**, not per-SNP effective N.

| Endpoint | R13 manifest cases | R13 manifest controls | Original DF13 case ICD-10 / inclusion | Original DF13 control definition | Released R13 GWAS? | Role |
| --- | ---: | ---: | --- | --- | --- | --- |
| I9_STR_EXH | NA | NA | I63 and I64, excludes I636; includes I9_STR_EMBOLIC | Control exclude I9_TIA | **NO** (definition exists; not listed among GWAS files) | Clinically preferable hemorrhage-excluded stroke, BUT NO R13 GWAS result; I64 is *unspecified stroke*, so even this definition cannot be called imaging-confirmed AIS |
| I9_STR | 36,461 | 446,552 | **I61**, I63, I64 | Control exclude I9_STR_EXH and I9_TIA | YES | Broad stroke phenotype **including intracerebral hemorrhage**; NOT directly eligible for strict GIGASTROKE AIS replication |
| I9_STR_EMBOLIC | 1,919 | 450,257 | I631 (cerebral infarction due to embolism) | Control exclude I9_CEREBVASC | YES | Candidate embolic-stroke direction replication only; NOT proven TOAST CE equivalent |
| I9_STR_EMBOLIC_WITH_FLIMMER | 1,041 | 432,403 | Includes I9_STR_EMBOLIC AND requires I9_AF | **Control conditions: not I9_AF** | YES | AF-associated embolic-stroke candidate; compare cautiously because the control group explicitly excludes AF and differs from base embolic-stroke controls |
| I9_ICH | 5,576 | 446,545 | I61 | Control exclude I9_STR_EXH and I9_TIA | YES | EXCLUDE, hemorrhagic stroke |
| I9_SAH | 4,365 | 446,966 | I60 | Control exclude I9_STR_EXH and I9_TIA | YES | EXCLUDE, subarachnoid hemorrhage |

**Correction to earlier shortlist:** The earlier keyword-driven manifest shortlist listed I9_STR, I9_STR_EMBOLIC and I9_STR_EMBOLIC_WITH_FLIMMER only. Official DF13 clinical source reveals **I9_STR_EXH**, but **its public R13 GWAS is not listed** in the pinned manifest. It is a definition/available GWAS mismatch, not a missing record to invent or impute. The I9_STR definition is mixed because it directly includes I61 intracerebral hemorrhage; do not use it as clean ischemic stroke. R13 has no demonstrated published LAA/SVD-specific GWAS in these verified endpoints. These are clinical definitions, not ground-truth TOAST stroke subtypes, and sample overlap with the GIGASTROKE 2022 consortium remains unverified. Do not claim independent replication without overlap accounting.

**GIGASTROKE original corrected 44-cohort N_eff scenarios are descriptive only** (pinned source SHA256 7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4):

| Stroke subtype | Pooled case/control harmonic illustration | Sum of cohortwise known-control harmonic illustrations | Relative difference | Cases with missing cohort control count |
| --- | ---: | ---: | ---: | ---: |
| AIS | 236,505.833 | 195,570.717 | 17.3083% | 103 |
| LAA | 25,464.041 | 23,149.571 | 9.0892% | 0 |
| CE | 42,841.160 | 39,605.585 | 7.5525% | 28 |
| SVD | 27,094.551 | 23,758.723 | 12.3118% | 39 |

**Do not mislabel these as lower/upper bounds** for true SNP-level effective sample sizes; both omit variant-specific contributing studies/weights and some cohort control denominators are unobserved. The 44-cohort source comparison demonstrates sensitivity to sample-size definition but does not license assigning either number to the original 11-column stroke GWAS for primary LDSC h2 or rg. No new LDSC h2/intercept/rg or MR/coloc was executed.

**Next research gate:** If FinnGen CE-direction replication is desired, request/obtain the official accessible R13 summary-stat files (not only the metadata), verify SHA/build/ALT allele and proper phenotyping, assess statistical power with 1,919/1,041 cases and possible overlap with GIGASTROKE, and avoid equating AF-enriched case/control ascertainment to TOAST CE. For primary AIS replication, identify another independently source-verifiable stroke GWAS that excludes hemorrhage and provides variantwise N / credible N method, instead of substituting FinnGen I9_STR. Run original GIGASTROKE MR/coloc only under its separate statistical/effect-harmonization QC: LDSC N_eff HOLD is not a blanket prohibition on all methods. Draft PR #284 remains unmerged.
