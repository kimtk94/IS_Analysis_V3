# GSE301517 independent third-cohort TM RNA audit — actual verified results (2026-10-31)

## Exact reproducible evidence

- Official primary GEO: [GSE301517](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301517), part of Singh et al. 2026 (chromatin topology, IOP/POAG).
- **Final 3-study signed comparison Actions [#41186736618](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41186736618) SUCCESS**, [artifact #13037335888](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41186736618/artifacts/13037335888) actually retrieved, opened, examined.
- First source attempt [#41186409050](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41186409050) failed **transient HTTP 403** on the GEO CPM; bounded official-source retry [#41186499319](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41186499319) confirmed source. The initial third-cohort comparison [#41186622446](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41186622446) succeeded; latest extends it to the official published **100nM 48h** comparator. History preserved.
- Reproducible source: `scripts/ocular/ocular_gse301517_third_cohort_preflight.py` and `scripts/ocular/ocular_gse301517_crossstudy_cpm.py`. Existing reusable on-demand workflow `.github/workflows/ocular-third-tm-gse301517.yml`. No new P-number workflow proliferated.
- Source official GEO SOFT SHA256: `a73281e9a59cffa3f37205f9281f3ad59baa50f6a6ef399b4c4153ba73c3187e`.
- Official submitted **CPM gzip SHA256: `717f4c9861de4d6b169439ef56a781f50d93ec918c456a9d0e61873e4a6e5f09`**, compressed **7,905,333 bytes**, complete gzip CRC PASS; 85,298 distinct Ensembl IDs, 24 TSV columns, 14 source sample-expression columns. NOT 85,298 testable protein-coding genes.
- **21 GSM** genuine official SOFT records: 7 separate Resource reference samples (TM62 3, TM71 2, TM83 2) excluded from DEX contrast; **7 Control + 7 DEX** technical sample profiles: TM62 4+4, TM70 2+2, TM72 1+1. **Three independently named TM cell strains, NOT fourteen independent human donors**. No independent donor identity validation.
- Official GEO treatment protocol: 80% confluency, 2h 2% FBS preconditioning, **0.5 µM (500nM) dexamethasone in 2% FBS for 7 days**, ethanol-matched vehicle, medium/treatment replacement every other day. Existing GSE313736 500nM BeamA **48h** and Mehrotra 100nM bulk **48h** are not a randomized 7d-vs-48h or 100-vs-500nM factorial experiment.

## Numeric method and actual results

Methods: match `ENSG` IDs (no fuzzy gene name union) to frozen original GSE313736 **16,050** primary edgeR genes; **15,876** IDs matched, all had finite numeric CPM values across all 14 external samples, zero missing. In each of TM62, TM70 and TM72, average available replicate CPM by condition; compute `log2((mean Dex CPM + 1)/(mean Vehicle CPM + 1))` **within each strain**, then **equal-weight the three strains**. Pseudocount = 1 CPM, fixed. GSE301517 data are normalized CPM: **no new edgeR/DESeq2 test, no independent gene-wise p-value/FDR**, no false 14-donor degrees of freedom. The prior 500nM 48h result is a count-based donor-blocked edgeR log2FC and cannot be interpreted as directly effect-size interchangeable with a CPM + pseudocount ratio.

| Frozen candidate | 500nM 48h edgeR log2FC | 100nM 48h published log2FC | 500nM 7d mean of 3 strain CPM log2 ratios | 7d strains matching the 500nM 48h sign |
|---|---:|---:|---:|---:|
| **MYOC** | **+4.02569** | **+1.85959** | **+4.34878** | 3/3 |
| **TGFBI** | **−1.86921** | **−0.84295** | **−1.52948** | 3/3 |
| **MMP2** | **−0.58907** | **−0.35730** | **−0.75592** | 3/3 |
| **ACTA2** | **+0.46950** | **+0.37375** | **−0.39881** | 1/3 |
| FKBP5 (primary-only P9 reference) | +3.98401 | NOT USED/NOT VERIFIED in P9 two-study selected-marker output | +3.81151 | 3/3 |
| EGR3 | +0.74018 | −0.82498 | −1.04088 | 0/3 |
| LPAR1 | +0.74562 | −0.65992 | −0.26517 | 1/3 |
| SDC1 | −0.55354 | +0.63147 | −0.19110 | 2/3 |
| TFPI | +0.61464 | −0.72295 | −0.10354 | 1/3 |

- **Three-study same sign: 3/4 core** (MYOC, TGFBI, MMP2), **ACTA2 exception**; 7d vs primary alone **10/16** original P9 panel signed, with 3/4 core matching across all three strain directions.
- All four P9 secondary ECM genes (IGFBP3, LUM, NTM, VCAM1) are **mean signed agreement** in 7d vs 500nM 48h, but VCAM1 only **2/3 strains** same; thus do not claim four perfectly strain-stable.
- **EMT 194 shared genes**: 7d vs original 500nM 48h **156/194 (80.41%)** same sign. 7d vs independent 100nM 48h **145/194 (74.74%)** same. **ALL THREE** studies same sign **124/194 (63.92%)**. Gene-level comparison is descriptive; **not** a calibrated third-study pathway-replication p-value.
- Across all matched gene symbols in source, 7d mean arm CPM≥1 for **13,269** common genes; restricting both original 48h edgeR and 7d CPM ratio to abs(effect)≥.25 and 7d CPM≥1 gives **2,350/2,745 (85.61%)** same sign. These filters are exploratory and select on both outcomes; do **not** use this fraction as an unselected genome-wide replication rate.
- 7d intensity per strain and tiny 1-CPM pseudocount can distort zero/low expression signals. Strong genes MYOC/TGFBI/MMP2 were each sign-consistent in all three strains; low expression targets require separate assay and pseudocount sensitivity. Equal-weight strains avoids TM62 technical-replicate pseudoreplication but **does not prove each strain corresponds to a unique biological donor**.

## Updated thesis interpretation

**Good evidence:** glucocorticoid exposure associated expression of MYOC↑, TGFBI↓ and MMP2↓ in three **different public TM experimental cohorts** with different doses/durations/culture conditions, and broad EMT/ECM transcript program signed concordance.

**Important negative context:** ACTA2 48h↑/7d↓, and EGR3/LPAR1/TFPI strong differences by cohort. This is biologically interesting but **cannot be labeled time-dependent reversal or dose-dependent reversal**, because long vs short studies differ in TM strain identity, cell subpopulation, treatment duration, library chemistry and analysis model. No within-donor randomized time×dose contrast.

**Novelty risk:** both independent source papers already report DEX, ECM and transcription biology; a meta-analysis or cross-study gene atlas must explicitly demonstrate novel uncertainty resolution (cell-state/experimental-unit selection bias, reproducibility modeling and useful candidate tiering) rather than claim basic MYOC/EMT detection as discovery.

**GO:** computational 3-study descriptive signed replication, source/matrix/QC, prioritized 3-core candidate hypothesis, P9 16-panel re-tiering.
**HOLD:** prospective wet-lab same-donor dose/time function, protein and MMP2 enzyme activity, cross-strain origin independence, genuine effect-size comparability, annotation cutoff sensitivity and source author's full 63/17 DEG crosswalk.
**STOP:** causal 7-day vs 48-hour difference, validated steroid-induced ΔIOP predictor, any microbial causal mediation from unrelated eye-surface cross-sectional data.

## Next immediate work — no more unnecessary workflows

1. **High priority:** GSE313736 BeamA confidence 0.5/0.75 and alternative cell-state annotation sensitivity, report donor N and 3-core vs ACTA2/SDC1/EGR3 direction. This is needed for manuscript-grade causal framing.
2. **Check exact publication overlap:** authors' reported 63↑/17↓ original GSE313736 high-effect genes against our 83, with actual source list. Do not assume counts match gene names.
3. **Optional:** a new low-cost, frozen-source cross-study figure heatmap across core 4, secondary ECM, and the 194 EMT genes, clearly marked 48h edgeR vs 7d CPM pseudocount log2.
4. **Only conditional:** GSE301516 ATAC and GWAS/actual native glaucoma TM once provenance and genomic harmonization are resolved; no new clinical glaucoma classification claim.
