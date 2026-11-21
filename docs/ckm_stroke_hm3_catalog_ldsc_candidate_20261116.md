# CKM × Stroke — MRCIEU vs Zenodo HapMap3 + conservative GIANT BMI candidate

**Evidence date: 2026-11-16** · Execution: Orca on NGKR-SAL-TH. Original public datasets streamed and independently pinned by existing source-integrity SHA-256 / Zenodo MD5. Scientific interpretation remains **HOLD** throughout.

## 1. Exact reference identity and marker comparison

| Data | Source | Pinned provenance |
| --- | --- | --- |
| MRCIEU HM3 | [MRCIEU/ldsc GRCh37 commit e9f9739c](https://github.com/MRCIEU/ldsc/blob/e9f9739c9c74ad5b4c653a6409a703cc005c62dc/GRCh37/data/hm3/w_hm3.snplist) | 17,170,716 bytes, SHA-256 `4ed5a031f88edff7f94663841160d0fff750c45d4d55a630515b4d6de082c49e` |
| Zenodo bundled HM3 | [Zenodo 8182036](https://zenodo.org/records/8182036), `eur_w_ld_chr/w_hm3.snplist` | 17,264,312 bytes, member SHA-256 `ec73fca0b696e8beba465b51e52911676fcade375bcc9475d99c0ec30509d3ed` |
| Canonical 1–22 EUR LD score | Same Zenodo `eur_w_ld_chr.tar.gz` | 33,357,890 archive bytes, posted MD5 `e2f16343c4cfaa76caa7d0c03d26b489`, verified gzip/22 scores/22 M_5_50 |
| GIANT original BMI | [GIANT Yengo 2018 corrected gzip](https://giant-consortium.web.broadinstitute.org/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz) | 46,754,937 bytes, source SHA-256 `0d6ed0ea97870916b830ccae349df94ca6f3cc68c025c79e784729af7f7136a4` |

MRCIEU original includes **1,204,345** usable unique SNP IDs after excluding the conflicting `NA` placeholder. Zenodo HM3 has **1,217,311** unique IDs, no duplicate/NA label.

| Reference comparison | Observed |
| --- | ---: |
| MRCIEU ∩ Zenodo rsIDs | 1,204,285 |
| MRCIEU only | 60 |
| Zenodo only | 13,026 |
| Common IDs, unordered allele pair identical | **1,204,285 / 1,204,285** |
| Strand-complement candidates among common pairs | 0 |
| Incompatible allele pairs among common pairs | 0 |
| Canonical EUR LDscore 1–22 union of unique rsIDs | 1,290,028 |
| MRCIEU ∩ canonical LDscore | 1,172,984 |
| Zenodo HM3 ∩ canonical LDscore | 1,184,461 |
| **MRCIEU ∩ Zenodo HM3 ∩ LDscore** | **1,172,984** |
| Exact allele agreement on triple intersection | 1,172,984 / 1,172,984 |

**Interpretation:** The catalogs are **not interchangeable** (60 and 13,026 exclusive IDs), but the shared SNP allele sets agree. Use their strict triple intersection as the *technical candidate panel*. This does not imply matching GRCh37 BP/strand or outcome GWAS input suitability. The LDscore mirror is pinned to Zenodo, not independently proven byte-identical to the old Broad `.bz2` bundle.

Implementation: `scripts/ckm_stroke/audit_hm3_zenodo_catalog_diff.py`. Local report: gitignored `artifacts/ckm_stroke/hm3_zenodo_catalog_diff_20261116.json`. Uses full MRCIEU SHA and Zenodo archive MD5, parses all 22 canonical LDscore IDs and allele-pair sets without extracting/retaining raw catalogs. Nine offline fixtures PASS.

## 2. Conservative BMI LDSC candidate input

From the **2,336,269** published GIANT BMI rows:
- Exact duplicate **9** extra records removed only after verifying all original fields. No discordant duplicate group.
- **1,328,361** rows failed strict 3-way reference SNP identity/allele set selection.
- **1,007,899** unique BMI SNPs retained in the **candidate**, not a final scientific analysis dataset.
- Candidate includes **14** of the original **44** `P=0` underflow records, represented as signed `Z=BETA/SE`; no P-value reconstruction or `P=1e-300` imputation.
- No palindromic alleles passed the conservative common reference panel, and no candidate allele mismatch detected. The `Freq_Tested_Allele_in_HRS` source frequency was **not** used as if it were the whole-cohort effect-allele frequency.
- Candidate columns: **`SNP N Z A1 A2`**, where A1 is original GIANT `Tested_Allele`, A2 original `Other_Allele`, and N source variant-level sample size.
- Candidate plaintext SHA-256: `ed8f713827cd3281bea5f88a9cb42a3eabfecb1ee67c912ca7e86f151e4399bf`.
- Orca compressed candidate size: **9,660,069 bytes**. Local ignored path: `artifacts/ckm_stroke/bmi_3way_conservative.sumstats.gz`.
- **Independent local re-read:** 1,007,899 unique rsIDs; correct header, finite signed Z (range **−46.875 to +46.75**), valid nonpalindromic alleles, per-variant positive integer N, full plaintext SHA-256 matched. No raw/candidate file uploaded to GitHub.

This is a **derived technical candidate** only. `ready_for_ldsc=false`, `scientific_status=HOLD`, `run_ldsc=false`. A technical output does not clear the UK Biobank sample overlap, cohort heterogeneity, allele-direction verification with all outcomes, HRS reference frequency limitation, or LDSC regression/intercept reliability.

## 3. Reproduction and enforcement

```powershell
# From the Orca research worktree.
py -3 -m unittest discover -s tests -p test_ckm_stroke_hm3_zenodo_catalog_diff.py -v
py -3 -m unittest discover -s tests -p test_ckm_stroke_prepare_bmi_ldsc_candidate.py -v
py -3 -m unittest discover -s tests -p test_ckm_stroke_verify_bmi_ldsc_candidate.py -v
py -3 -m scripts.ckm_stroke.audit_hm3_zenodo_catalog_diff --out artifacts/ckm_stroke/hm3_zenodo_catalog_diff_20261116.json
py -3 -m scripts.ckm_stroke.prepare_bmi_ldsc_candidate --out artifacts/ckm_stroke/bmi_ldsc_candidate_qc_20261116.json --emit-candidate artifacts/ckm_stroke/bmi_3way_conservative.sumstats.gz
py -3 -m scripts.ckm_stroke.verify_bmi_ldsc_candidate --candidate artifacts/ckm_stroke/bmi_3way_conservative.sumstats.gz --provenance artifacts/ckm_stroke/bmi_ldsc_candidate_qc_20261116.json --out artifacts/ckm_stroke/bmi_candidate_independent_check_20261116.json
```

CI: `.github/workflows/ckm-stroke-hm3-triple-candidate.yml` runs offline 3-suite regression on relevant pushes; the full original source analysis occurs **only** under explicit `[hm3-triple-full]` commit marker or manual dispatch. Full CI generates a temporary candidate, checks its contents, and deletes the gzip. Only small JSON metrics are uploaded. Use its candidate *plaintext* SHA-256, not compressed gzip bytes, for cross-platform reproducibility.

## 4. Remaining scientific gate

- Verify conversion and harmonization against complete, officially sourced EUR Stroke subtype outcome GWAS headers and reference alleles.
- Reconcile genomic build/BP for common rsIDs, non-unique rsIDs in outcome files and ambiguous allele orientations.
- Investigate sample overlap of GIANT/UKB and GIGASTROKE EUR, intercept and regression reference lineage.
- Perform sensitivity check for standardized BMI phenotype scale, per-variant N and effect-size rounding; run LDSC `h2` diagnostics only **after** the data/governance and scientific input gate is signed off.
- SBP and T2D original access remain separate HOLD; do not bypass data use terms.

**Draft PR #284 remains unmerged**, and no LDSC `h2/r_g` run is authorized by this technical QA.
