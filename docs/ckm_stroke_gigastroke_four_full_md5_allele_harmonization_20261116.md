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
