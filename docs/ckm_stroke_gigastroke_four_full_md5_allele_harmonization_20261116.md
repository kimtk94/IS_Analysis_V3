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

New code adds \`chr_attrition_pre_munge\` to the original four-GWAS MD5-checked streaming audit, partitioning the independently pinned 1,007,899 BMI 2018/HM3/EUR LD-score candidate coordinates for each chromosome 1–22. Every chromosome records candidate positions, unique outcome positions, invalid duplicate positions, absent positions, valid aligned/numeric positions and incompatible/invalid positions. The sum must equal the original whole-genome counts, and the identical candidate chromosome profile must be observed in all four outcome comparisons.

The new independent summarizer \`scripts/ckm_stroke/audit_ldsc_chr_attrition.py\` verifies the existing whole-genome original publisher MD5 and allele-count reviewer first; it never trusts arbitrary counts from an unsigned JSON. Its own tests and the per-coordinate partition tests reject altered chr indices, missing chromosomes, incorrect subtotals, invalid source MD5 and subtype-specific panel mismatch.

**Important:** pre-munge overlap is *not* the SNP set used in LDSC's h² regression. These counts are upper-bound candidates conditional on the Yengo 2018-derived panel, not independent whole-genome outcome HM3-universe counts. No unapproved case-control N, LD-score weight merge, effect Z threshold, large-chi-square exclusion, MAF or GC handling is inferred. Final regression SNP matching and h² remain **NOT_RUN/HOLD**. The official LDSC documentation warns about analyses with fewer than 200,000 merged regression SNPs and favors 600,000+ for stable precision (not a substitute for input approval).
