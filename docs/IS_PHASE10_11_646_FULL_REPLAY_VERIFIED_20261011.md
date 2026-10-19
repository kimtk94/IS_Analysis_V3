# IS Phase10/11 — FULL original 646-source SNP coloc reproduction audit

**FULL_646_NUMERICAL_REPLAY independently validated against Drive-origin source artifacts, complete historical raw ABF master and server Git R script.**

## Evidence and scope

- 646/646 original SNP-level gene–tissue analyses from **four original BBJ GWAS regions** reproduce their archived five H0–H4 posterior probabilities at tolerance <=1e-8.
- Original 646 results span 43 gene IDs, 16 tissue datasets, and are **not** 646 independent loci nor the expanded 80 provisional candidate windows.
- Original run settings: coloc 5.2.3, p1=p2=1e-4, p12=1e-5, GWAS N=174686, cases=22664; GTEx QTL first-scalar N.
- Maximum absolute deviation logged by original R run: 5.55e-16; archived posterior table roundtrip difference 0.
- All 646 SNP denominators, unique SNP counts, first-scalar QTL N and QTL N ranges match historical source index/master.
- All 646 have SNP-varying QTL N; 646 input hashes recorded, and 30 local cached full input TSV hashes match the corresponding Colab source manifest.
- Both complete output status and PASS-only TSV bytes identical; terminal run log reports 646 PASS and 0 failures/missing.
- Colab run full 646 file execution was independently verified **as a computation**; this was not independent cohort molecular association evidence.

## Full original source H4 distribution

| Historical ABF posterior | Number of gene–tissue tests |
|---|---:|
| PP.H4 >= 0.50 | 6 |
| PP.H4 >= 0.75 | 2 |
| PP.H4 >= 0.80 | 0 |

## Interpretation and remaining scientific gates

- The historical computation is reproducible. It does NOT justify interpretation of any gene as proven to mediate stroke risk.
- The single-causal-signal coloc ABF prior and effect of tissue/gene testing multiplicity require explicit sensitivity analysis.
- All historical gene–tissue tests used QTL first-scalar N despite per-SNP N differences. Reproduction is not a robustness test of that assumption.
- EAS504 GWAS-side reference LD cannot substitute for ancestry/cohort-appropriate GTEx QTL LD in joint multi-signal colocalization.
- GRCh37/GRCh38 identity mapping, GWAS/eQTL effect-allele orientation and sample-size sensitivity remain separate validation gates.
- R3_7 human adult control brain donor expression supports only healthy cell localization; FGF5/C4orf22 missing features stay NOT_ASSESSABLE, not biological zeros.
- The expanded 80 provisional GWAS windows / 2,225 gene IDs remain a separate broader discovery branch.

## Drive provenance

- Verified original output folder: https://drive.google.com/drive/folders/1AKS8A2YlE_UrsINmUMdDtVsQpz_vz_Rn
- Colab source/index/master, full 646 status and source SHA256 manifests are cross-linked by SHA256 in IS_646_FULL_SOURCE_AUDIT.json.
- Full Drive source artifacts and all raw SNP files are outside Git; versioned server snapshots remain untouched.

## Original highest-H4 tests (descriptive, not gene causality)

| Gene | GTEx tissue | PP.H3 | PP.H4 |
|---|---|---:|---:|
| CALHM2 | Brain Cerebellar Hemisphere | 0.0377 | 0.7977 |
| FGF5 | Brain Cerebellar Hemisphere | 0.1408 | 0.7787 |
| NEURL1 | Brain Cerebellar Hemisphere | 0.1713 | 0.6776 |
| C4orf22 | Brain Cerebellar Hemisphere | 0.1125 | 0.6599 |
| C4orf22 | Brain Cerebellum | 0.0802 | 0.5597 |
| INA | Brain Anterior cingulate cortex BA24 | 0.0897 | 0.5364 |

The six rows above are selected from the original **646 gene–tissue tests** and are not multiplicity-adjusted. The top PP.H4 of the historical screen is 0.7977, so zero tests reach 0.80.

## Direct Drive output references (verified folder readback)

- [Original complete outputs folder](https://drive.google.com/drive/folders/1AKS8A2YlE_UrsINmUMdDtVsQpz_vz_Rn)
- [Full run manifest](https://drive.google.com/file/d/1ZPqOizan7H3AKXbSaTKwKibRENxWXbVy/view)
- [Complete 646-row coloc posterior status](https://drive.google.com/file/d/17WtpLgmzksUIKoaX8RF2UGYDyg2i1lqe/view)
- [646 input SHA256 checksums](https://drive.google.com/file/d/1h0ZrgtE074Q4XWqvfVgZ_98eZPxqXE3b/view)
- [R runtime log](https://drive.google.com/file/d/1EmqkTUU6GxYgrdBl5PoJQawfDS5Qu3hM/view)

## Release gates

- **G0 historical numeric reproducibility: COMPLETE** (all 646 source matched, full source identifiers and original executable R script pinned).
- **G1 effect-allele / GRCh37–38 / ancestry harmonization: PENDING**; successful exact rerun does not prove cohort-specific allele validity.
- **G2 original coloc model robustness: CONDITIONAL / PENDING**, including SNP-varying eQTL N, tissue selection, prior calibration, multi-signal model and appropriate GTEx QTL LD.
- **G3 independent molecular and disease-cell regulatory validation: PENDING**. Healthy temporal-lobe single-cell expression is only descriptive.
- Expanded 80 provisional GWAS components and 2,225 positional genes are a **different candidate denominator** and have not been colocalized by this legacy four-locus screen.

- Current code: scripts/is/audit_is_646_full_colab_source.py
- Standalone source proof JSON and server source artifacts (not in Git): results/is/audits/legacy_coloc_646_full_colab_20261011_v1/validation_v3/

**No causal gene-mechanism claim has been promoted by this milestone.**
