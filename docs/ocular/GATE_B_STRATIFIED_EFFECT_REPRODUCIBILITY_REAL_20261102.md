# GSE313736 Gate B: effect/abundance-stratified sensitivity from actual gene-level edgeR data

**Execution date:** 2026-11-02. **User decision:** do not contact original authors; author80 overlap stays a non-blocking limitations item. Proceed with available public transcriptomic studies and transparent uncertainty.

## Executed primary evidence

- [GitHub Actions #41389639386](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41389639386) — **SUCCESS**; two independent CI jobs PASS (stratification fixture test + **real archived gene-level edgeR statistics**).
- [Actual archive #13116138830](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41389639386/artifacts/13116138830) includes new effect-stratified **PNG Figure**, full 16-row strata TSV, P9 subset gene effects, source SHA256 manifest and Markdown report. Real source files were fetched on the runner from existing study Actions only. **No 18-matrix original 10x raw file redownload; no donor IDs, patient metadata or per-donor UMI matrices uploaded.**
- Reproducible: `scripts/ocular/ocular_gate_b_stratified_effect_reproducibility.py` and `.github/workflows/ocular-gateb-stratified-real-aggregate.yml`.
- Comparison A (same biological 57 donor pairs): BeamA score ≥0.75 from [#41085982185](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41085982185) versus score ≥0.50 from [#41199532259](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41199532259).
- Comparison B (fixed **same 55** donor pairs): BeamA versus scPred TM_Fibro from [#41206289766](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41206289766). Original combined workflow FAILED on a strict 16-gene comparator, while both real count-level models PASS; the corrected downstream reviewer [#41207583296](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41207583296) SUCCESS.

## Frozen numerical checks and stratified biological interpretation

| Metric | 57 donors: BeamA 0.75 vs 0.50 | 55 donors: BeamA vs scPred TM_Fibro |
|---|---:|---:|
| Genes testable in both | 16,050 | 13,975 |
| Overall log2FC Pearson r | **0.99911597** | **0.94990493** |
| Shared-gene same sign | **15,851 (98.76%)** | **12,634 (90.40%)** |
| BH FDR <0.05 in both | **13,563** | **10,255** |
| Median absolute log2FC difference | **0.01124618** | **0.10189687** |

**Effect magnitude strata below are based only on the baseline model, not a preregistered independent partition**. Denominator includes only Ensembl genes testable in both relevant comparisons.

| Baseline | 57 donor classification score: same sign | Fixed-55 annotation: same sign |
|---|---:|---:|
| abs(log2FC) <0.25 | **96.659% (5,957 tested)** | **78.258% (5,349 tested)** |
| 0.25 ≤ abs(log2FC) <1 | **100% (8,728)** | **97.666% (7,583)** |
| 1 ≤ abs(log2FC) <3 | **100% (1,282)** | **99.898% (984)** |
| abs(log2FC) ≥3 | **100% (83)** | **100% (59)** |

**Baseline expression quartiles:**

| Expression quartile | 57 donor score sign agreement | Fixed-55 annotation sign agreement |
|---|---:|---:|
| Q1 (lowest logCPM) | 97.907% | **87.064%** |
| Q2 | 98.679% | 89.236% |
| Q3 | 99.153% | 91.929% |
| Q4 (highest logCPM) | 99.302% | **93.389%** |

This stratification provides an important restriction on the interpretation of aggregate high Pearson correlation. With a different cell labeling approach (same 55 donors) the **smallest baseline transcript effects (<0.25) exhibit 21.742% direction disagreement**, versus zero direction reversals among the 59 common-tested genes with baseline abs(log2FC)≥3. The latter group is defined by baseline effect size and is subject to effect-based selection; it is **not independent replication** nor a newly discovered 59-gene author set. Within-score threshold comparison was near-identical but both versions share original classification algorithm, donors and source libraries.

MYOC↑, TGFBI↓, MMP2↓ remained signed/FDR significant in both models for both sensitivity contrasts. ERFE was **not testable** under the scPred TM_Fibro expression filter: do not turn its lack of a second estimate into evidence of absence or a nonsignificant negative result.

### Key study limitations and novelty ceiling

1. The original treatment and vehicle were prepared in **condition-specific 10x libraries**. Neither donor blocking nor two within-cohort label/filter sensitivity analyses identify a treatment-specific causal effect separately from technical processing.
2. The original analysis found **13,717/16,050 (~85.46%) BH FDR<0.05**, a large signal fraction that supports reporting **both adjusted statistical significance and minimum effect size**, while avoiding claims that 85% of genes have biologically meaningful independent treatment effects.
3. Near-zero log2FC sign agreement is intrinsically sensitive to small numeric variation; no binwise hypothesis test or independent p-value is inferred here. Stratification is descriptive and exploratory.
4. This analysis improves transparent evaluation of **selection/annotation robustness**, not evidence of altered protein abundance, MMP2 enzymatic activity, TM outflow, or clinical steroid-induced glaucoma susceptibility.
5. Original GSE313736 authors' 63↑/17↓ gene identities remain unavailable. **No author contact will be requested or sent under this project plan**. The author80-to-our83 overlap is not a thesis submission blocker; it must instead be clearly declared unverified and no result can imply author-level exact replication.

## Publication handling

- Main results: 57-donor and 55-donor global concordance + three preselected candidate signed effects, with uncertainty caveat.
- **Supplementary Figure (new):** two panels comparing sign concordance by baseline abs(log2FC) category and baseline logCPM quartile, with actual tested gene denominator annotated.
- Supplementary Table: 16 strata rows with counts, both-FDR count, median |Δlog2FC|; annotated P9 candidate table.
- **Next priority without author outreach:** examine batch/selection confounding boundary quantitatively using *existing non-identifying aggregates*; if those cannot identify causality, preserve HOLD rather than inventing corrections. Strengthen independent external comparison based on verified per-strain and per-timepoint source mapping (Mehrotra/GSE301517) only when compatible.

**Outcome:** computational sensitivity/stratified figure **GO**; between-library causality **HOLD**; exact original authors' DEG80 overlap **HOLD / explicitly non-blocking**; clinical prediction **STOP**.
