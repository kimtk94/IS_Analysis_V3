# GSE301517 7-day TM experiment — individual strain contribution and strict cross-source direction

**Executed:** 2026-11-02 KST. **Study scope:** published human TM transcript data only, not original GSE313736 authors' unavailable DEG80 list; no investigator outreach required or planned.

## Actual source and CI

- [GitHub Actions #41423487666](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41423487666) **2/2 SUCCESS**, self-test and actual previously verified publisher-derived gene-level aggregate audit.
- [Actual result artifact #13130831706](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41423487666/artifacts/13130831706) has 16 gene × 3 strain × leave-one-out TSV; 194 EMT gene strain-agreement TSV; explicit source SHA-256 manifest and aggregate report. No FASTQ files or patient/gene-by-donor count matrix downloaded or exported.
- Script `scripts/ocular/ocular_gse301517_strain_loo_aggregate_audit.py`; on-demand workflow `.github/workflows/ocular-gse301517-strain-loo-aggregate.yml`.
- Read only the previously validated original [GSE301517 source-run #41186736618](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41186736618) archive #13037335888, and verified SHA256 of the frozen 4-core and 194-gene input. The original source CPM SHA256 was already verified in that archival pipeline.
- **Biological unit caveat:** three named TM strains `TM62`, `TM70`, `TM72`, with original technical samples per treatment arm `4`, `2`, `1` respectively. The three strains' actual donor independence is not confirmed. Mean effects give each **strain** equal weight; the raw replicate counts should never be treated as seven independent donors or a 14-donor comparison.
- **Estimator caveat:** 7-day values are `log2((mean CPM_DEX+1)/(mean CPM_Control+1))`, averaged over three strains. They are not model-derived gene log2FC, not directly comparable magnitudes to acute GSE313736 edgeR or published Mehrotra DESeq2, and have no independently computed 7-day FDR.

## Main 3-study candidate levels — real archived numerical values

| Gene | GSE313736 500nM 48h edgeR | Mehrotra 100nM 48h published | 7d TM62 | 7d TM70 | 7d TM72 | All leave-one-strain-out 7d signs match GSE313736? |
|---|---:|---:|---:|---:|---:|---|
| **MYOC** | +4.026 | +1.860 | +4.822 | +3.648 | +4.576 | **YES, 3/3** |
| **TGFBI** | −1.869 | −0.843 | −2.300 | −1.344 | −0.945 | **YES, 3/3** |
| **MMP2** | −0.589 | −0.357 | −0.835 | −1.271 | **−0.162** | **YES, 3/3** |
| ACTA2 | +0.469 | +0.374 | **+0.005** | −0.760 | −0.442 | **NO, 0/3** |

The three previously selected candidate genes **MYOC↑, TGFBI↓ and MMP2↓** maintain the direction in each of the three 7-day TM strain groups. Each exclusion of one group also retains the 48-hour direction; no one 7-day strain is needed to determine the *sign*. MMP2 is weak in TM72 (−0.162 on the CPM+1 scale), so this is **directional** not guaranteed effect-size, functional or statistical stability. ACTA2's TM62 +0.005 is effectively near zero and should not be interpreted as a convincing independent positive effect. Its 7-day aggregate is negative in two other strains, and both 48h sources were positive: biological context and study differences are confounded.

## 16-gene preselected panel sensitivity

- 7-day three-strain **equal-mean sign matches** the 500nM 48h original in **10/16**.
- All **three separate 7-day strain** log2 ratio signs match the original in **7/16**.
- All **three leave-one-strain-out equal mean** signs match original in **8/16**. One gene can pass LOO with two positive strain ratios and one mild negative, so the latter condition is less strict than three of three per-strain unanimity. The three definitions must not be collapsed.
- The canonical 3 (MYOC, TGFBI, MMP2) satisfy **all three** definitions.
- The sign-only comparison does **not** establish differential expression p-value in any third-study strain; normalized CPM is observational.

## Hallmark EMT 194-gene source-strain audit

| Criterion | Gene count | Fraction of 194 |
|---|---:|---:|
| Three study aggregate signs all match | **124** | **63.9%** |
| All three individual GSE301517 strains match GSE313736 acute sign | **103** | **53.1%** |
| **Both criteria: study mean + 7-day 3/3 unanimity** | **90** | **46.4%** |

The original normalized 7-day strain-alignment counts among the 194 matched genes were:
- 0 of 3 strains same as GSE313736: **12 genes**
- 1 of 3: **25 genes**
- 2 of 3: **54 genes**
- 3 of 3: **103 genes**

**Interpretation:** the 124/194 descriptive three-study statistic substantially overstates the fraction that also passes strict 7-day strain-level sign unanimity; only **90/194** passes both. Because the 194 Hallmark EMT genes were selected after the fact and the genes are correlated, these are **not 194 statistically independent replication trials**, and no binomial or pathway p-value is supported. The three source studies also differ in drug dose, exposure length, cultured strain composition and estimator.

## Thesis use and evidence gate

- **Main Figure 1** MYOC/TGFBI/MMP2 signed three-source cross-study candidate prioritization remains supported, with newly documented three-strain individual and leave-one-out sign checks.
- **Supplementary Table** provides 16 candidate genes × 3 strain estimates and all three LOO ratios; include accompanying table column for whether the same effect is seen in each source.
- **Figure 3 / supplement** state **124/194 aggregate signed** and now **90/194 strict strain-unanimity** as an additional exploratory attrition gate; do not relabel 90 genes as a new discovery.
- **GO**: original public source identity, same-gene signs and robust within-three-strain directional core on the published 7-day processed data.
- **HOLD**: 7-day strain-independent donor details and gene-wise p/FDR, full time/dose causal comparison, drug-vs-library technical confounding in GSE313736, protein and TM function.
- **STOP**: clinical steroid glaucoma risk predictor, mechanistic ECM/fibrosis causal claim or original source author's 80-gene exact reproduction.
- **Original author email not required or sent**. Exact 80-gene source overlap remains non-blocking HOLD.
