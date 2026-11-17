# CKM × ischemic stroke — first-party study-cohort overlap risk: Locke BMI 2015 vs GIGASTROKE 2022

**Date:** 2026-11-17 KST. **Study question:** Can Locke et al. (2015) European BMI be treated as *sample-independent* from the four GIGASTROKE 2022 European ischemic stroke GWAS datasets? **Answer based on published source data: NO — shared recruiting cohorts are documented, but actual overlapping participants are not quantified.** This is a provenance / risk-screening analysis, *not* an estimate of the LDSC genetic covariance intercept or sample overlap fraction. `scientific_status=HOLD`, `ready_for_ldsc=false`.

## 1. Verified primary original supplements

The Locke 2015 Nature paper [Genetic studies of body mass index yield new insights for obesity biology](https://doi.org/10.1038/nature14177) reports 125 studies in the full GWAS+Metabochip study. We obtained its **actual Supplementary Table 1**, not an inferred list from authors/affiliations, from the institutional EMBL-EBI BioStudies Europe PMC public archival package [S-EPMC4382211](https://www.ebi.ac.uk/biostudies/studies/S-EPMC4382211):

- Official [`NIHMS668049-supplement-2.xlsx`](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-EPMC/211/S-EPMC4382211/Files/NIHMS668049-supplement-2.xlsx) — **1,626,884 compressed bytes**, complete SHA-256 **`dadd561bba67dbda9cbe33331a2fea195e638f0951ace75706149d137d2171f8`**; 30 worksheets. `S1-Study QC` is the original cohort/study design, ethnicity and sample size table, with sections for prior GWAS, updated GWAS, Metabochip and additional studies. `S2-Study QC` lists genotype/imputation/association methods. Primary PMC article: [PMC4382211](https://pmc.ncbi.nlm.nih.gov/articles/PMC4382211/).
- GIGASTROKE Mishra et al. (2022) [Nature](https://doi.org/10.1038/s41586-022-05165-3), publisher **corrected 2022 Supplementary Table 1** [`41586_2022_5492_MOESM2_ESM.xlsx`](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05492-5/MediaObjects/41586_2022_5492_MOESM2_ESM.xlsx) — **756,910 compressed bytes**, complete SHA-256 **`7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4`**; 44 named EUR discovery cohort rows 5–48 of `Suppl.Table1`, EUR total row 49. Columns K AIS, L CES, M SVS, N LAS, P controls. Total row AIS 62,100 / CE 10,804 / SVD 6,811 / LAA 6,399 / controls 1,234,808. Source was rechecked against official GWAS Catalog accessions (AIS GCST90104540, LAA GCST90104542, CE GCST90104541, SVD GCST90104543) earlier.
- Both primary XLSX files are stored **only in ignored Orca `artifacts/ckm_stroke/references/`**, not tracked in Git. The script reads both using standard-library bounded ZIP/OOXML, validates full SHA-256 and unique expected sheet profiles, and publishes an **aggregate cohort-label evidence JSON only**. No genotype, individual records or raw GWAS association rows are imported.

Two source archival pathways were investigated: direct PMC supplement access and the old Nature extension returned an HTML block, whereas EMBL-EBI BioStudies exposed the actual same-numbered Locke supplementary table file. The BioStudies file was byte-checked and sheet-verified; HTML/error text was never treated as XLSX.

## 2. Directly corresponding recruitment cohort labels — 10 of 44 GIGASTROKE EUR entries

The strict cohort-label category means the original publishers independently name the same specific study or wave. It **does not prove that a particular individual was used in both analyses**. Counts below are respective original author-reported sample sizes; they are **not cross-study overlapping participant counts**.

| Study family | Locke Supplement 1 (label; N as printed) | GIGASTROKE corrected Supplement 1 (source label; AIS cases / controls) | Label interpretation |
|---|---|---|---|
| AGES | AGES; 3,207 | AGES; 99 / 2,882 | Same named study |
| ARIC | ARIC; 8,108 | ARIC_EA*; 523 / 8,241 | Same named study, European ancestry |
| CHS | CHS; 3,228 | CHS_EA*; 520 / 2,615 | Same named study, European ancestry |
| FHS | FRAM (Framingham Heart Study); 8,094 | FHS_EA*; 214 / 8,062 | Same named study family |
| deCODE | deCODE; 26,799 | DECODE; 4,483 / 254,000 | Same source label; later project-scale recruitment may differ substantially |
| Rotterdam I | RS-I; 5,744 | RS1_EA*; 603 / 5,033 | Same named cohort wave |
| Rotterdam III | RSIII; 2,006 | RS3_EA**; 48 / 2,901 | Same named cohort wave |
| WGHS | WGHS; 22,888 | WGHS_EA*; 587 / 22,582 | Same named study |
| SHIP | SHIP; 4,092 | SHIP*; 41 / 1,915 | Same named study |
| Health ABC | Health ABC; 1,655 | HABC*; 144 / 74 | Same named study |

**Critical:** GIGASTROKE CHS and ARIC original cases/controls are defined for ischemic stroke outcomes, whereas Locke original N is a BMI phenotype-specific sample count. The table must not be used to infer that the smaller count is the actual shared sample: the participants, phenotype QC, genotyping and inclusion dates are not linked by identifier. The rows may overlap partially, completely, or not at all in individuals. The evidence supports **cohort-based sample-overlap risk**, not quantified individual sample overlap.

## 3. Cohort-family/umbrella links — 6 further source entries (not proof of same participant sample)

| GIGASTROKE source label | Locke named source(s) | Interpretation |
|---|---|---|
| PROSPERlogistic_EA* | PROSPER/PHASE | Same named trial/study family; different analysis subsets possible |
| HUNT** | HUNT 2 | Different study waves / participant mix unknown |
| RS2_EA* | RS-I | Rotterdam family, but RS2 ≠ RS-I; **not a direct same-cohort match** |
| EPICCVD* | EPIC-Obesity Study (also EPIC-Norfolk T2D present in Locke) | EPIC umbrella; source subcohort distinction crucial |
| EstBB** | EGCUT | Estonian biobank institutional source family; participant identity unresolved |
| UTARTU-EstBB** | EGCUT-370 | Estonian Genome Center University of Tartu source family; participant identity unresolved |

**Do not merge family-only links into the direct cohort count.** Linked studies can have disjoint samples or partially overlapping waves, especially prospective follow-up recruitment and multi-site consortia.

The remaining **28 of 44 GIGASTROKE source labels** are not matched by this strict predefined list. The absence of a matching label is **not evidence of zero participant overlap**; unlisted aliases, nested consortia and overlapping sample recruitment may exist.

Notable categorical case: `UKBiobank` is present in GIGASTROKE corrected Table 1 (AIS cases **2,211**, controls **369,419**) but not named as a Locke 2015 input cohort in the original published S1 table. Locke 2015 is an excellent **UK Biobank avoidance sensitivity source**, yet it still cannot be treated as independent of GIGASTROKE overall.

## 4. Scientific LDSC sample-overlap and effective-N policy

These original first-party tables establish **several potential shared study cohorts outside UKB**, which invalidates a blanket assertion that "pre-UK Biobank BMI = zero overlapping participants." The actual shared-person count \(N_s\) and covariance of BMI with ischemic stroke among the shared cases/controls cannot be identified from the two public cohort summaries.

The [LDSC official FAQ](https://github.com/bulik/ldsc/wiki/FAQ) and [Bulik-Sullivan et al. cross-trait LDSC derivation](https://pmc.ncbi.nlm.nih.gov/articles/PMC4797329/) make clear that constraining the cross-trait genetic covariance intercept requires an independently justified sample overlap value or secure proof of zero overlap. **Accordingly, a default constrained zero covariance intercept is NOT authorized**. An unconstrained intercept *diagnostic* may eventually be planned after the original effect-size scale, LDSC variant panel, N values and trait h² are validated, but do not interpret an estimated intercept as a direct headcount of overlapping subjects.

Outcome per-SNP N is absent in the four GIGASTROKE originals. Phenotype-wide \(4N_{case}N_{control}/(N_{case}+N_{control})\) and the sum of cohortwise harmonic values differ; see [cohort N_eff audit](ckm_stroke_gigastroke_cohort_overlap_neff_20261117.md). The two calculations are illustrations only. The Locke 2015 file has a noninteger per-variant reported `N` for 11,506 independently selected EUR HM3 variants; N and P-value effects require format and GC correction review. Do not silently copy the total BMI study size to all variants or round fractional N as if exact participant counts.

GIGASTROKE's published subtype h² precision differs strongly (original Cov-LDSC observed-scale h²/SE about AIS 9.8, CE 6.17, SVD 3.0, LAA 1.67). Therefore confirm h²/SE and intercept diagnostics before any rg interpretation, especially LAA. Cross-trait analysis must not be promoted from technical SNP harmonization alone.

## 5. Reproducible computation and evidence

- Script: `scripts/ckm_stroke/audit_locke_gigastroke_cohort_overlap.py` (first-party source SHA, original 30-worksheet Locke header, original 44-row GIG EUR table, exact/near-match source label tiers, required aliases, EUR totals, no shared-person estimate, `scientific_status=HOLD`).
- Tests: `tests/test_ckm_stroke_locke_gigastroke_cohort_overlap.py`, **8 offline fail-closed cases** including publisher source checksum, missing cohorts, alternate waves, EUR ancestry label, LAA/CE swapped columns, changed source case totals, repeated Locke study label, inaccessible source.
- GitHub workflow: `.github/workflows/ckm-stroke-locke-gigastroke-cohort-overlap.yml`. Fast fixture tests on push, independent original XLSX check triggered only by `[cohort-overlap-full]` or workflow dispatch. Raw publisher XLSX is read into memory only; uploaded artifact is aggregate JSON. 
- Orca ignored report: `artifacts/ckm_stroke/locke_gigastroke_cohort_overlap_20261117.json`; technical source and alias result `PASS_FIRSTPARTY_TWO_SUPPLEMENT_COHORT_NAME_RECONCILIATION`.

**Status:** Cohort-label provenance/technical audit **PASS**. Exact person-level overlap `UNKNOWN`; zero covariance intercept `NOT_APPROVED`; final LDSC h², rg, MR, coloc `NOT_RUN`; scientific `HOLD`. Draft PR #284 stays Draft, unmerged.

### Next permissible work

1. Design a **single-trait Locke 2015 BMI and subtype outcome h²/intercept diagnostic protocol** that requires original effect scale, per-variant N and GWAS study ancestry/build QC approval before runtime execution.
2. For sample-overlap sensitivities, compare the existing 2018 GIANT+UKB exposure with independently selected 2015 Locke exposure **using identical methods and an explicitly recorded intersection/regression SNP universe**, not arbitrary differences in SNP counts.
3. Only consider any cross-trait genetic covariance / rg after h² precision, intercept behavior and cohort-assumption limitations are approved and transparently reported; never treat cohort-name matches as an exact \(N_s\) estimate.
