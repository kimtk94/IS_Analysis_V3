# P5r — three parallel public aggregate audit gates: actual results, 2026-10-31

## Execution provenance
- **[GitHub Actions #41172214095](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41172214095): 3/3 MATRIX JOBS SUCCESS**. Independently reviewed mode-specific uploaded ZIP contents and job logs.
- `thresholds`: artifact [#13031144426](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41172214095/artifacts/13031144426).
- `timecourse`: artifact [#13031403190](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41172214095/artifacts/13031403190).
- `consensus`: artifact [#13031313532](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41172214095/artifacts/13031313532).
- Auditing code `scripts/ocular/ocular_p5r_threeway_aggregate_audits.py`. Matrix workflow `.github/workflows/ocular-p5r-threeway-aggregate.yml`. Temporary PR trigger reverted to `workflow_dispatch` after actual run.
- Inputs are prior **public gene-level aggregated files only**: P5p run 41161006241, P5m run 41160429359, P5o run 41160840453. No donor counts, private identity, or clinical IOP outcomes accessed.

## A — post-hoc candidate cutpoint sensitivity
24 actual combinations: mean primary logCPM cutoffs 0,1,2,3,4,5 and minimum absolute magnitude applied to *all three* effect estimates (primary BeamA, external bulk TM, conservative 9-pool model) 0.25,0.5,0.75,1.0; pooled library BH-FDR<.05 and same direction as primary are retained.

| Mean logCPM cutoff | Count at all three abs(log2FC)>=0.5 | Count at all three abs(log2FC)>=0.75 |
|---|---:|---:|
| >=0 | 13 | 1 (ERFE) |
| >=1 | 11 | 0 |
| >=2 | **9** | **0** |
| >=3 | **8** | 0 |
| >=4 | **7** | 0 |
| >=5 | **5** | 0 |

At >=3, ADAMTS15 (mean logCPM 2.50) is excluded. At >=4, ABHD3 (3.52) is additionally excluded. At >=5, EGR3 (4.29) and SDC1 (4.43) are also excluded. **This exposes post-hoc threshold sensitivity, not poor experimental replication of those genes.** At >=2/effect>=.5 the exact original nine P5p symbols are reproduced with their Ensembl IDs; no new gene is independently discovered.

## B — significance-aware official external 100nM TM time course for 12 targets
Official published 100nM bulk TM gene contrasts, *not* a matched-dose experiment with GSE313736 500nM BeamA:

| Gene | 1h external log2FC / FDR | 6h external log2FC / FDR | 48h external log2FC / FDR | Interpretation |
|---|---|---|---|---|
| **EGR3** | **-0.357 / 0.00930** | -0.504 / 0.9996 | **-0.825 / 0.00993** | 1h and 48h signed decreases significant; 6h not significant |
| **ERFE** | +0.170 / 0.99999 | **+0.864 / 0.04026** | **+0.804 / 0.01242** | 6h and 48h signed increases significant |
| CHAC1 | -0.299 / 0.99999 | -0.305 / 0.9996 | **+0.575 / 0.01174** | Early sign differs, but early estimates nonsignificant |
| PARM1 | -0.089 / 0.99999 | +0.375 / 0.9996 | **-0.709 / 0.000614** | Early effects nonsignificant |

Among exactly 12 candidates: 1/12 FDR significant at 1h, 1/12 at 6h, **2/12 any early significance**, **8/12 at least one early raw sign flip** versus external 48h. A timepoint-specific FDR result is not a formal dose×time or cross-time contrast; published summary effects without covariance do **not** establish a temporal reversal. In particular, EGR3's 6h FDR does not establish lack of effect, and sign flips of nonsignificant effects are not valid reversal evidence.

## C — P5o×P5p source-by-source agreement checks

All 96 unique Ensembl IDs join, all 96 external/primary 48h effect signs oppose. Original descriptive QC is reproduced *exactly*: **33/96** all-nine source-pool signed concordance; **96/96** all-nine equal-pool excluded means signed; **75/96** pool-bootstrap descriptive 95% intervals exclude zero.

Among the 12 prespecified P5p targets:
- **10/12** same sign in all 9 original source pools. **DDIT4L and RAB32: 8/9** (not 9/9).
- **12/12** keep signed pooled mean across all nine omitted-pool calculations.
- **12/12** have original pool-bootstrap descriptive 95% intervals excluding zero.
- Selected per-donor direction counts from 57 include: STK17B 56/57, ABHD3 55/57, LPAR1 54/57, RAB32 54/57, DDIT4L 45/57, SDC1 47/57, ERFE 50/57, CHAC1 50/57, PARM1 48/57. This is **within-source descriptive consistency**, not independent biological or clinical replication.

## Decision, limitations and next gate

**GO**: the aggregate-only computational crosschecks, external timepoint FDR-aware summaries and P5p candidate cutpoint grid are reproducibly completed. These are validation of previous numerical summaries, **not three new independent cohorts or new evidence of clinical relevance**.

**P5q formal 12-gene leave-one-pool-out edgeR QL remains PENDING** until [actual workflow #41170641379](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41170641379) completes and its 216 gene-level rows / 24 by-model summaries are manually reviewed. P5r's descriptive pool bootstrap cannot substitute for P5q inferential FDR.

**HOLD**: 100nM vs 500nM dose effect, treatment vs library-preparation batch disentanglement, post-treatment BeamA selection effects, unique TM mechanism and any protein effect. To resolve dose attribution, require independently cultured donor-matched Vehicle/100/500nM ×1/6/48h with randomized sample/library batches and pre-specified donor-blocked dose×time interaction plus qPCR/protein QC. Low-expression ERFE/CHAC1/PARM1 need detection-limit/assay reliability checks.

**STOP** clinical steroid-induced delta-IOP prediction: no treated patient-level response phenotype or independent patient cohort.

## Artifact security check
P5r uploaded files are mode-specific public gene-level aggregate TSV/JSON only. They do not contain original UMI matrices, donor IDs, raw read files, or per-donor trajectories.
