# CKM × Ischemic Stroke — Locke 2015 independent EUR BMI SNP universe (not selected from 2018)

**Date:** 2026-11-17 KST. **Research classification:** alternate historical GIANT BMI exposure for sample-overlap sensitivity, not the primary exposure and not yet an approved LDSC input. Scientific status remains **HOLD**.

## 1. Scientific question and independence boundary

Earlier technical comparison of Locke 2015 to Yengo 2018 only tested the **1,007,899 SNPs chosen from the 2018 BMI GWAS** and found 1,007,858 overlaps. That restricted comparison cannot represent an *independent 2015-selected exposure universe*. This new pipeline deliberately selects variants from the **original Locke 2015 BMI EUR source itself**, using *only externally verified reference data*, **without reading any 2018 BMI candidate variants, association effects, p-values, or N**.

Source:
- [GIANT official Locke et al. 2015 European BMI gzip](https://giant-consortium.web.broadinstitute.org/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz), 35,869,809 bytes, complete original source SHA-256 `ddca3b30985e3428f4762f8e8aadd385504627b9b825cff1ecdc4a5e2e358b44` (reproduced by independent prior Linux CI [#43983463946](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43983463946)).
- [Official GIANT 2012–2015 file schema](https://giant-consortium.web.broadinstitute.org/images/b/b9/Data_File_Description_2012-2015_GIANT_GWAS_Summary_Statistics.docx): effect allele `A1` is hg19 plus strand, source N is per-SNP reported observations, `Freq1.Hapmap` is **HapMap CEU** reference frequency (not cohort-wide exposure EAF); P is genomically controlled and association coefficients limited to 2 significant digits.
- [Pinned MRCIEU GRCh37 HM3, exact git commit](https://raw.githubusercontent.com/MRCIEU/ldsc/e9f9739c9c74ad5b4c653a6409a703cc005c62dc/GRCh37/data/hm3/w_hm3.snplist), original SHA-256 `4ed5a031f88edff7f94663841160d0fff750c45d4d55a630515b4d6de082c49e`.
- [Zenodo official record 8182036 EUR LD score mirror](https://zenodo.org/records/8182036): 33,357,890 compressed bytes, **publisher posted MD5 `e2f16343c4cfaa76caa7d0c03d26b489`**, with bundled `w_hm3.snplist`, canonical 1–22 `.l2.ldscore.gz` and `.l2.M_5_50`. Historical `6_old` / noncanonical members are ignored.

## 2. 2015-only independent SNP selection (actual Orca source run)

**Exact technical selection**: eligible rsID must exist in MRCIEU HM3 ∩ Zenodo HM3 ∩ 22 canonical EUR LD score SNP sets, with an **identical unordered A1/A2 pair across both HM3 sources**. The actual eligible reference contains **1,172,984 unique rsIDs**. Then match these IDs against **Locke 2015 original**; drop *all* records for every source-duplicated rsID (13 distinct duplicate IDs, 14 extra rows); exclude ambiguous palindromic alleles, pair-incompatible SNPs, bad beta/se/N/P. Preserve Locke 2015 original **A1 effect direction and unrounded fractional source N**; compute signed `Z_FROM_BETA_SE` from the original *rounded* BETA and SE. Do not impute original P from that Wald Z.

| Source/selection QC | Actual independent source result |
| --- | ---: |
| Original Locke 2015 association rows | 2,554,637 |
| Eligible three-reference rsID panel | 1,172,984 |
| Original rows outside reference triple | 1,483,425 |
| Excluded original duplicate rsID records | 27 (14 additional rows across 13 IDs) |
| Excluded reference allele-pair discordance | 431 |
| Excluded palindromic among matched source rows | 1 |
| **2015 independently selected candidate SNPs** | **1,070,753** |
| 2015 candidate SNPs with fractional source N | **11,506** |
| 2015 source P ≤ 5×10⁻⁸ among selected | 839 |

Deterministic intermediate plain-text SHA-256 (after canonical row serialization): `a1cc5ec06cc81d2bb0a583eb76c0f7e81784cd5dc3f622270de61023f95de4bf`. Output in *ignored* `artifacts/ckm_stroke/locke_2015_independent.candidate.tsv.gz`. Header is **`SNP N_REPORTED Z_FROM_BETA_SE A1 A2`**, deliberately **not** `SNP N Z A1 A2` to prevent accidentally treating fractional reported N and a rounded-coefficient Z as approved LDSC munged summary statistics.

Important: No 2018 BMI file participates in source selection. The same MRCIEU+Zenodo reference *families* used by 2018 pipeline are shared methodologically; reference reuse is not 2018 effect selection.

## 3. Whole-panel EUR LD score coordinate assignment and independent GRCh37 sentinel audit

The original Locke 2015 file has **no chromosome/BP columns**. For the 1,070,753 independent 2015 candidates, chromosome+base-pair coordinates were assigned from 22 **canonical EUR LD score reference members** and checked for duplicate rsIDs or coordinate collisions. This assignment is **reference-derived, not an independent 2015 original coordinate**.

| Verification | Actual result |
| --- | ---: |
| Canonical EUR LD score rows across autosomes | 1,290,028 |
| 2015 candidate rsIDs present in EUR LD score | **1,070,753 / 1,070,753** |
| 2015 candidate rsIDs missing from canonical reference | **0** |
| Repeated 2015 candidate rsIDs in reference | **0** |
| Duplicate 2015 candidate reference coordinates | **0** |
| Complete 22 canonical LDscore and 22 M_5_50 | **PASS** |
| Ensembl GRCh37 forward-strand independent 1-bp sentinels | **44/44 matched allele pair**, 2 per chromosome |
| Genome-wide 2015 variant effect sign against Ensembl | **NOT verified** (sampled SNP-pair test only) |

44 SNPs were chosen by smallest SHA256(`locke-2015-EUR-LDscore-Ensembl44-original-only-v1|rsID`) per chromosome, independently of GWAS p values. The public Ensembl API checked each source 2015 A1/A2 pair against forward GRCh37 reference nucleotides. This verifies a stratified **sample**, not all 1.07M original forward-strand allele pairs or the actual effect BETA sign.

### Comparison to 2018 after independent 2015 generation (post hoc only)

After selection and QC were complete, the independent derived 2015 and previously pinned derived 2018 panels were **compared**, never mixed during selection:

| Set | Unique SNP IDs |
| --- | ---: |
| 2015 independently selected | 1,070,753 |
| 2018 conservative GIANT+UKB selected | 1,007,899 |
| **Common rsIDs** | **1,007,858** |
| **Only in 2015** | **62,895** |
| Only in 2018 | 41 |

The extra 2015 variants do **not** themselves imply a better-powered GWAS, additional causal biology, or independence from the GIGASTROKE cohorts. They arise from differing source SNP availability and QC.

## 4. Reproducibility, QC contract, and interpretation

- Independent selection script: `scripts/ckm_stroke/prepare_locke_2015_independent_panel.py` with 7 offline synthetic tests in `tests/test_ckm_stroke_locke_2015_independent_panel.py`. Checks full original pinned SHA-256, 13 known ambiguous duplicate IDs, all 2,554,637 rows and 14 extra duplicate records, HM3/Zenodo reference hashes, strict allele match and preserved fractional per-variant N, candidate plaintext SHA256.
- Independent coordinate + 44 Ensembl reference script: `scripts/ckm_stroke/audit_locke_2015_reference_coordinates.py` and 7 offline tests `tests/test_ckm_stroke_locke_2015_reference_coordinates.py`. Confirms all candidate ID/chr:BP unique, 22 reference file completeness and publisher archive MD5, 44 real external genome bases.
- GitHub CI: `.github/workflows/ckm-stroke-locke-2015-independent-panel.yml`. Offline tests on changes; source-heavy jobs opt-in by `[locke-2015-independent-full]` commit tag or manual dispatch. Temporary association table deleted before artifact publication; **only compact aggregate JSON and 44 public sentinel points** uploaded. No source GWAS, inferred full effect table or participant identifiers committed.
- Local evidence: `artifacts/ckm_stroke/locke_2015_independent_panel_20261117.json`, `artifacts/ckm_stroke/locke_2015_independent_reference_20261117.json` (ignored).

```powershell
py -3 -m unittest discover -s tests -p test_ckm_stroke_locke_2015_independent_panel.py -v
py -3 -m unittest discover -s tests -p test_ckm_stroke_locke_2015_reference_coordinates.py -v
py -3 -m scripts.ckm_stroke.prepare_locke_2015_independent_panel --emit-candidate artifacts/ckm_stroke/locke_2015_independent.candidate.tsv.gz --out artifacts/ckm_stroke/locke_2015_independent_panel_20261117.json
py -3 -m scripts.ckm_stroke.audit_locke_2015_reference_coordinates --candidate artifacts/ckm_stroke/locke_2015_independent.candidate.tsv.gz --provenance artifacts/ckm_stroke/locke_2015_independent_panel_20261117.json --out artifacts/ckm_stroke/locke_2015_independent_reference_20261117.json
```

**Scientific review state:** `HOLD`, `ready_for_ldsc=false`; no new LDSC h², intercept, rg, MR or coloc computed. The N reported by Locke 2015 is sometimes fractional; full LDSC formatting, N setting, signed association scale, adequate LDSC h² precision, and cohort overlap with GIGASTROKE remain unapproved. In particular "before UK Biobank era" **does not mean** "independent of all GIGASTROKE participants". Continue cohort-wise overlap reconciliation (e.g. shared non-UKB GIANT and GIGASTROKE cohorts) and prespecified h²/intercept QC before considering rg. No Draft PR merge.
