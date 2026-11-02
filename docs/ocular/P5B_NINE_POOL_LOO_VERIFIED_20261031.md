# Ocular P5b — VERIFIED nine leave-one-pool-out edgeR results (2026-10-31)

## Evidence and stage
- GitHub Actions **[run #41085982185](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41085982185)**: job **SUCCESS**, including downloading all 18 original public 10x matrices, creating ephemeral 57 matched-donor BeamA pseudobulk counts, fitting the original robust paired edgeR model, executing all **9** leave-one-pool-out QL models, deleting private donor counts and uploading aggregate-only results.
- **[Artifact #12993824564](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41085982185/artifacts/12993824564)**: `p5b_pool_sensitivity.json`, `p5b_leave_one_pool_out.tsv`, `p5b_top_gene_leave_pool.tsv`, corrected anonymous `p5b_MDS_annotated.png`, `p5b_pool_robustness.png`, primary gene-level results. No donor IDs or individual donor count matrices are in the artifact.
- Dataset GSE313736: 57 paired biological human cultured trabecular meshwork donors, BeamA Chen2025 fine classification score >=0.75, >=30 cells per donor per arm, same 18 original matrices as P4/P5. The first reference donor-paired QL fit uses `robust=TRUE`; all nine sensitivity fits use `robust=FALSE` for compute/runtime. Thus model-level FDR counts are **method sensitive**; interpret primarily signed logFC stability and non-identical FDR denominators carefully.

## Actual 9-pool LOO results
| Excluded pool | Remaining donor pairs | Genes tested | BH-FDR<0.05 | Correlation with primary log2FC | Primary significant common-gene FDR retention |
|---|---:|---:|---:|---:|---:|
| A | 49 | 15,944 | 13,323 | 0.997265 | 96.27% |
| B | 50 | 15,997 | 13,391 | 0.998206 | 96.30% |
| C | 51 | 16,008 | 13,580 | 0.999333 | 97.64% |
| D | 50 | 16,063 | 13,536 | 0.998875 | 97.34% |
| E | 50 | 16,032 | 13,532 | 0.999066 | 97.72% |
| F | 51 | 16,026 | 13,506 | 0.999300 | 97.69% |
| G | 50 | 16,003 | 13,507 | 0.998993 | 97.64% |
| H | 51 | 16,212 | 13,822 | 0.998737 | 98.25% |
| I | 54 | 16,149 | 13,892 | 0.999427 | 99.13% |

- **9/9 successful QL refits**, minimum gene-level log2FC Pearson r=**0.99726488**, median r=**0.99899329**. The top **100** primary genes all preserved effect direction in all nine refits (100/100 each).
- FDR retention is defined among **the primary significant genes present in the re-filtered common gene set**, not all 13,717; it ranges from 96.267% to 99.125%. Original full primary gene set: 16,050 tested, 13,717 FDR<0.05.
- **Pool G**: despite 92.743% original DEX BeamA barcode recall vs >=95% for other libraries, G omission yields 50 pairs, 13,507 FDR-significant, r=**0.998993**, retention **97.6385%**. Therefore G does not dominate global signed expression effects.
- Pool I reduced sample count: 54 pairs remain; r=0.999427, FDR retention=99.125%, but **original source attrition is high (7 described donors→4 final metadata)** and is not cured by LOO.

## Corrected MDS figure and limitation
- The actual `p5b_MDS_annotated.png` shows **all 114 anonymous matched pseudobulk samples** and 57 paired line segments, unlike the earlier invalid blank MDS. This is a real run output and was visually checked after download.
- First MDS dimension: *unadjusted, single-factor* linear R² for treatment=**0.95557721**, for pool=**0.00978584**. Second dimension treatment R²=**0.01824173**, pool R²=**0.08924798**. These are separate regressions, not a multivariable decomposition, and are not evidence that treatment causes 95.6% of gene-expression variability.
- Vehicle and DEX separate almost perfectly on MDS dimension 1. This is compatible with true strong DEX response **and also** an unmeasured treatment-correlated library/batch difference. The design sequenced DEX and Vehicle in different source 10x libraries; 57 donors are **not 57 separately prepared treatment-paired libraries**.
- The chart `p5b_pool_robustness.png` uses a [0,1] y-axis and appears visually saturated at the upper border; use the numeric TSV rather than that figure for subtle r differences.

## Scientific novelty and next gates
- Prior 2022–2026 experiments already report HIPK2, PTK2B, FKBP5, FAM107A, ZBTB16, ATP10A, EDIL3 and MAOA–MIEF1–DRP1 steroid response. Reproducing their DEX expression direction is **not new mechanism discovery**. See `P5B_LITERATURE_GENE_CROSSWALK_20261031.md`.
- **[P5d conservative 9 paired experimental-library pool unit edgeR](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41090827791)** also PASS: 21,081 genes tested, 13,600 BH-FDR<.05, logFC correlation with donor model 0.99624485 on 16,049 common genes; 80.77% of common primary significant genes remain significant. This strengthens within-study numerical robustness but cannot separate treatment from library prep.
- **[P5e classifier selection](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41091146021)** PASS: global BeamA fraction DEX higher by +2.406 percentage points, but median paired 59-donor fraction −0.421 percentage points (29 positive/30 negative), so cell-count aggregation and post-treatment annotation selection are potential bias.
- **[P5f original donor attrition](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41091777826)** PASS: 74 paper donors→70 official GSM descriptions→59 final metadata→57 high-confidence BeamA analysis, reasons for removal not verified. 316,799 originally reported singlets vs final 169,398 metadata cells.
- Independent replication PRJNA1243143 public ENA has 46 runs but lacks interpretable per-sample strain/time/treatment map: **HOLD**. No cultured-donor clinical steroid-induced ocular pressure outcome: **STOP** clinical adverse effect prediction.

## GO/HOLD verdict
- **GO** reproducible in-vitro molecular statistical association and *within-cohort* signed direction stability, regardless of any single pool.
- **HOLD** robust causal attribution to DEX versus condition-correlated 10x library effects, independent biological replication, truly novel molecular mechanism and dose-time harmonization.
- **STOP** steroid-induced clinical IOP genetic risk prediction from these cultured cell-only data.

## Remaining work, preregistered
1. Acquire verified control/DEX pair identity and processed counts for independent human TM RNA-seq from PRJNA1243143, or another validated external accession/clinical dataset. No unverified run-to-condition guessing.
2. Re-test donor-specific responses under treatment-independent broader cell-state definition and alternative BeamA classifier confidence thresholds to quantify selection effect; do not choose a threshold based on significant gene yields.
3. Evaluate the distribution of paired pool effects under 9 experimental-library replicates, including confidence intervals and gene-set benchmarks; explicitly separate 57-donor inference vs 9-library-unit experimental replication.
4. Make figure-level prior-paper comparison with full originally published gene sets (not only 5–10 positive controls), cross-check reference genome and ontology mapping, and prioritize novel, independently replicating candidates.
5. Protect all donor IDs, genetic traits and per-donor expression matrices. Public artifacts contain gene-level aggregate effects and anonymous global diagnostics only.
