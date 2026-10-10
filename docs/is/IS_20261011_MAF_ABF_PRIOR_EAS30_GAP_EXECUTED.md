# IS Research — BBJ SNP Frequency, ABF Prior Fragility and EAS Discovery Scope

**Date:** 2026-10-11 KST  
**Branch:** `research/is-broad-discovery-20261009`  
**Evidence level:** **VERIFIED computational input/provenance** vs **EXPLORATORY statistical interpretation**; **CAUSAL GENE NOT ESTABLISHED**.

## Executive conclusions

The original eQTL coloc window is now reproducible at SNP, gene-tissue, reference-allele, and source-GWAS levels, but it does **not** validate specific causal genes. The main next bottleneck is wider molecular-QTL coverage, not inflating the H4 prior.

All original 2,225 positional genes, 646 archived BBJ/GTEx ABF assays, and 8,956 unique SNPs are retained. The 646 assays evaluate only 43 of 869 East-Asian-positioned genes across four of 30 EAS regions. The remaining **826 EAS genes/26 wholly unscreened regions** are not negatives.

## 1. Original GWAS / LD input selection — frozen subset

| Source gate | Count |
|---|---:|
| Original per-locus Stage3 BBJ SuSiE v3 SNP rows | **10,065** |
| Unique source-verified BBJ/GTEx coloc SNPs | **8,956** |
| Of archived 8,956, native BBJ ALT effect/beta/EAF/chain concordant | **8,956** |
| Of archived 8,956, present in source 1KG EAS pvar, post-QC pvar and Stage3 | **8,956** |
| Of archived 8,956, BBJ minor allele frequency below 0.01 | **0** |
| Palindromic A/T or C/G | **1,243** |
| Palindromic with BBJ **or** source EAS ALT AF in 0.4–0.6 | **182** |
| Absolute BBJ-EAS ALT AF difference ≥0.10 | **581** |
| Absolute BBJ-EAS ALT AF difference ≥0.20 | **1** |

The retention fraction **8,956/10,065 ≈ 89%** applies only to Stage3-to-coloc comparison. It is **not** an estimate of how many original native GWAS variants survived genotype MAF QC; the Stage3 input already underwent selection. The 1KG EAS `EAS_AF` field is an original VCF population INFO annotation. It is **not** GTEx cohort EAF or a guarantee of molecular-QTL LD matching.

### Extreme population ALT frequency discordance

- **GRCh37 `10:104952499:C:T`**, BBJ locus `BBJ_IS_L002`: native BBJ `AF_Allele2=0.372646741166674`.
- Original 1KG EAS PVAR `EAS_AF=0.9792`. Independent read-only **PLINK 2 `--freq`** over the archived 504-person genotype panel yielded **`ALT_FREQS=0.979167`, `OBS_CT=1008`**.
- The same original BBJ Stage3 row carries `info=0.828592511538718`; alleles and variant ID were checked as C/T at the same GRCh37 coordinate.
- This is a **source-supported frequency discordance**; do not label it an allele strand flip or a causal result without further locus and sample-level provenance validation. Manual GWAS QC and cohort-composition review recommended.

## 2. Complete five-prior ABF sensitivity

These are identical 646 archived gene-tissue SNP sets and fixed `p1=p2=1e-4`; only `p12` changes.

| Shared-signal prior `p12` | H4 > H3 assays | PP.H4 ≥ 0.5 | PP.H4 ≥ 0.8 | PP.H1 > 0.5 |
|---|---:|---:|---:|---:|
| 1e-6 | 1 | 0 | 0 | 561 |
| 3e-6 | 7 | 2 | 0 | 558 |
| **1e-5 (original baseline)** | **40** | **6** | **0** | **551** |
| 3e-5 | 435 | 13 | 4 | 533 |
| 1e-4 | 565 | 86 | 12 | 203 |

High-prior (`1e-4`) H4≥0.8 occurs in **12 assays belonging to 7 distinct genes across 3 original BBJ loci**. At the original baseline **0 assays** pass H4≥0.8. The single largest original H4 for these historical genes is `CALHM2=0.79768355`, still below 0.8. Higher-prior H4 is not replication and must **not** be promoted to shared causal variants.

The H1-dominated fraction (551/646 at baseline) shows that in many *tested molecular-QTL comparisons* GWAS-only association dominates under the specified single-signal model. It does **not** prove those candidate genes lack eQTLs: power, sample size, tissue selection, variant ascertainment, ancestry or multi-signal architecture remain alternatives.

### Reproducible Figure

[**H3 vs H4, same 646 assays at baseline vs high p12**](../figures/is/IS_646_ABF_H3_H4_BASELINE_VS_HIGHPRIOR.svg)

- Source script: `scripts/is/render_is_646_h3h4_sensitivity.R`
- The dots are assays, not 646 independent replications; the H4=0.8 horizontal criterion is descriptive, not a formally calibrated causal-gene threshold.

## 3. Full candidate universe: EAS screen gap versus non-EAS pipelines

| Positional gene category | Number | Correct interpretation |
|---|---:|---|
| **All unique Ensembl gene IDs** | **2,225** | Broad candidate universe, no causal claim |
| Belonging to one of 30 `IS_XDATA_G....` East Asian regions | **869** | EAS-positioned genes |
| EAS genes with any archived BBJ/GTEx ABF assay | **43** | 646 tests in just 4 EAS loci |
| **EAS genes without archived ABF tests** | **826** | **Not tested; not negative** |
| Outside the EAS positional regions | **1,356** | Needs other outcome/geography-specific pipelines |
| Genes linked to EUR GWAS regions (including EAS overlap) | **1,556** | Do not add directly to 869; overlap = 200 |

**Thirty** distinct EAS regions are present; only **4** have any previously archived molecular-QTL tests, leaving **26** entirely unassayed in this ABF branch.

Coverage-gap-first scheduling (not GWAS association ranking or biological importance ranking):

| Region | Positional genes | Previously tested genes | Unscreened genes |
|---|---:|---:|---:|
| `IS_XDATA_G0022` | 112 | 16 | **96** |
| `IS_XDATA_G0025` | 92 | 0 | **92** |
| `IS_XDATA_G0015` | 71 | 17 | **54** |
| `IS_XDATA_G0018` | 53 | 0 | **53** |
| `IS_XDATA_G0030` | 53 | 0 | **53** |
| `IS_XDATA_G0026` | 50 | 0 | **50** |
| `IS_XDATA_G0028` | 50 | 0 | **50** |

Note: 869 EAS genes and 1,556 EUR genes are not disjoint; the total of their membership counts is 2,425 because **200 gene IDs occur in both**. This is separate from the 2,425 positional *gene-to-region relationships* previously counted. Do not conflate all relationships with EAS-plus-EUR intersections.

### Prioritized scientific work, without selection of causal genes

1. Build an **all-30-EAS-region audit** of outcome GWAS sites, exact genome builds, allele direction, LD reference population, and native nominal molecular-QTL data availability. Begin with large coverage gaps (G0022, G0025, G0015, G0018, G0030) only as *workflow scheduling*. Do **not** use number of positional genes as a biological evidence score.
2. Find eQTL/sQTL data that include **non-significant cis tests** for all overlapping SNPs and target genes; published significant-only matrices cannot support unbiased whole-region coloc input. Differentiate nominal/sparse QTL, expression below detection, and inaccessible genotype.
3. Preserve baseline `p12=1e-5` as the primary historical comparison and stricter/higher priors solely as sensitivity. Use paired original GWAS/eQTL alleles and sample sizes, source QC, and separately ancestry-matched molecular QTL LD for any SuSiE multi-signal follow-up.
4. Publish the 2,225-gene readiness matrix with **MISSING/NOT_ASSESSED** labels, not low biological scores. Keep the original 646 read-only and differentiate H3/H4 from adult-control single-cell expression localizations.
5. Evaluate exceptional source AF discordance individually, not by automatic flips or post hoc SNP exclusion designed to raise H4.

## 4. Executed code/data provenance

All files are **new derivative outputs**, not modifications of archived inputs:

- `scripts/is/audit_is_8956_maf_palindrome_reference.py`  
  `/srv/is-analysis/results/is/audits/is_8956_maf_palindrome_selection_20261011_v1/`
- `scripts/is/summarize_is_646_prior_h3h4_2225_genes.py`  
  `/srv/is-analysis/results/is/audits/is_646_prior_h3h4_2225_full_20261011_v1/`
- `scripts/is/build_is_eas30_qtl_gap_priority.py`  
  `/srv/is-analysis/results/is/audits/is_eas30_qtl_gap_priority_20261011_v2/`
- `scripts/is/render_is_646_h3h4_sensitivity.R`  
  `docs/figures/is/IS_646_ABF_H3_H4_BASELINE_VS_HIGHPRIOR.svg`

The first EAS scope output directory `v1` contains an incomplete intermediate attempt (manifest key error); use **v2** only. All new derivatives include source hashes, no source file overwrite.

**Offline CI-level regression:** the previous **37** provenance tests plus **9** new guard tests (new EAS source AF, 646 prior-grid integrity, EAS 30 region mapping): **46/46 PASS** executed on `taehoonserver` 2026-10-11. GitHub workflow source updated. Hosted GitHub Actions job status not independently attested.

## Scientific status

**VERIFIED:** original 8,956 SNP source allele coding and provenance, Stage3 and EAS reference presence, genotype-backed ALT frequency of one outlier, 646 archived ABF posterior values and prior sensitivity, source-linked counts for 30 EAS regions/2,225 distinct positional genes.

**EXPLORATORY:** 12 H4≥0.8 high-prior assays, healthy control brain cell-type localization, locus-first workflow sequencing, alternative multiomic mechanism predictions.

**NOT VERIFIED / BLOCKED:** molecular-QTL genotype-backed effect orientations and tissue-matched LD, European/Japanese cross-cohort harmonization beyond tested records, full unscreened QTL nominal files in 26 regions, outcome GWAS independent replication, ischemic-stroke case-control expression effects, formal causal gene assignment. Keep ALDH2 mechanistic/conditional inference **BLOCKED** until its independent branch gates clear.
