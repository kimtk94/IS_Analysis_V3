# Thesis Figures & Actual Validated Results — Ocular TM DEX (2026-10-31)

## Executed figure CI provenance

Final **[GitHub Actions #41208972011](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41208972011) SUCCESS**, [PNG/SVG/legends/table artifact #13045363678](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41208972011/artifacts/13045363678). Generated 3 PNG, 3 SVG, one Supplementary Table S1, figure legends and machine-readable source manifest with the SHA256 of all seven source files. Downloaded and verified the ZIP, checked all three images visually. The first figure run [#41208620110](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41208620110) **SUCCESS** computationally but captions visibly overlapped the bottom plot in rendered PNG; spacing was corrected and rerun. Do NOT use the older overlapped figures in manuscript.

Workflow \`.github/workflows/ocular-thesis-figures.yml\` has **workflow_dispatch only** after launch to avoid unnecessary action usage. Source \`scripts/ocular/ocular_thesis_figures_from_verified_artifacts.py\` uses Python/matplotlib and verifies frozen file SHA256 before generating images, allowing exact repeat of the figures even after source model versions change. The actual frozen input artifacts are:

| Figure input | GitHub Actions run | Source artifact ID | Role |
|---|---:|---:|---|
| GSE313736 48h + Mehrotra 48h + GSE301517 7-day gene-direction | [41186736618](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41186736618) | 13037335888 | 4-core + 194 EMT genes |
| Original BeamA confidence 0.75 vs 0.50 57-donor raw count paired edgeR | [41199532259](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41199532259) | 13041944372 | score sensitivity |
| Exact shared 55-donor BeamA vs scPred TM_Fibro raw paired edgeR | [41207583296](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41207583296) | 13045350639 | annotation sensitivity |

*The cross-annotation raw-source run #41206289766 produced valid two edgeR models but was flagged FAILED on the final post-hoc 16-candidate comparator because ERFE was filtered in scPred; the corrected aggregate-only reviewer #41207583296 is SUCCESS. Both facts are retained.*

## Figure 1 — independent 48h/7d DEX-associated transcript direction

| Gene | GSE313736 500nM 48h donor-paired edgeR log2FC | Mehrotra 100nM 48h published bulk TM log2FC | GSE301517 500nM 7d three-strain mean log2 CPM+1 ratio | Signed agreement all 3? |
|---|---:|---:|---:|---|
| **MYOC** | **+4.026** | **+1.860** | **+4.349** | **YES** |
| **TGFBI** | **−1.869** | **−0.843** | **−1.529** | **YES** |
| **MMP2** | **−0.589** | **−0.357** | **−0.756** | **YES** |
| **ACTA2** | **+0.469** | **+0.374** | **−0.399** | **NO** |

**Conservative Results paragraph draft:**
Three publicly available human trabecular meshwork transcriptomic experiments provided evidence of directionally concordant glucocorticoid-associated changes in MYOC, TGFBI and MMP2. MYOC increased after DEX exposure in all three experiments, whereas TGFBI and MMP2 decreased. ACTA2 showed positive effects in two 48-hour experiments but a negative log-ratio in the independent seven-day experiment, suggesting that its transcriptional response may depend on study context. Comparisons are descriptive because the sources differ in donor/strain composition, exposure duration, drug concentration, sequencing procedures, cell-state selection and analysis scale. In particular, the seven-day normalized CPM-based log-ratio is not an edgeR log2FC from a model matched to the 48-hour experiments, and no 7-day gene-level statistical p/FDR was independently estimated.

**What NOT to write:** "dose-dependent reversal", "temporal causal reversal", "validated novel genetic target", "DEX causes ECM fibrosis", "new prediction of steroid-induced glaucoma".

## Figure 2 — within-original-study annotation robustness

- **Panel A:** 57 identical paired donors: BeamA confidence ≥0.75 versus ≥0.50, direct raw UMI pseudobulk edgeR re-fitting. Across 16,050 shared genes log2FC Pearson **0.99911597**. Core MYOC +4.026/+3.908; TGFBI −1.869/−1.881; MMP2 −0.589/−0.595; ACTA2 +0.469/+0.419.
- **Panel B:** 55 **identical** paired donors: BeamA≥0.75 versus independently assigned scPred TM_Fibro cell annotation, two edgeR donor-blocked models. 13,975 common Ensembl genes, log2FC Pearson **0.94990493**; MYOC +4.006/+3.525; TGFBI −1.868/−1.941; MMP2 −0.591/−0.459; ACTA2 +0.481/+0.618. Of the prespecified 16 candidates, **15 were testable in both and 15/15 had the same direction and FDR<0.05**; ERFE was expression-filtered in scPred and must not receive a fabricated direction or p-value.
- Both comparisons are made **within the same GSE313736 source and share 10x libraries and donors; the label-specific sets can include overlapping cells**. Multiple tests do not count as independent cohort replication. Strong apparent significance across most genes requires careful discussion of potential condition/10x-library aliasing and post-treatment cell selection.

**Conservative Results paragraph draft:**
To characterize sensitivity to cell-classification choices, we recomputed donor-paired pseudobulks from the original raw UMI matrices using two BeamA classification-confidence thresholds and a second scPred TM fibroblast annotation. The 57-donor confidence-threshold analysis was nearly identical genome-wide, whereas the more demanding 55-donor cross-annotation analysis showed a correlation of 0.950 in log2 fold changes across 13,975 jointly testable Ensembl genes. The three leading candidates, MYOC, TGFBI and MMP2, retained direction and gene-level FDR significance under both sets of classification choices. These results support computational robustness of the RNA signatures to the specific tested annotation schemes but do not independently identify DEX biological causality because all comparisons reuse the same post-treatment expression experiment.

## Figure 3 — exploratory Hallmark EMT gene-sign counts

Gene set: 194 previously matched Hallmark EMT Ensembl genes; **141** same direction between the two 48-hour studies; **156** 500nM 48h vs 500nM 7d; **145** 100nM 48h vs 500nM 7d; **124/194 (63.9%)** share the same sign across all three. **Not independent gene Bernoulli outcomes or a calibrated pathway replication p-value.** No matrix deposition, gelatinase activity, tissue stiffness, outflow resistance or glaucoma clinical phenotype is observed.

**Conservative Results paragraph draft:**
Exploratory interrogation of 194 Hallmark EMT-annotated genes revealed broadly similar signed changes between transcriptomic experiments, with 124 genes retaining concordant directions in all three datasets. Because the gene set was examined following earlier enrichment findings and expression changes in gene networks are correlated, these fractions are descriptive rather than independent inferential tests. Moreover, the transcriptomic program cannot be equated with increased extracellular matrix deposition or functional glaucomatous outflow obstruction in the absence of independent biochemical and tissue-level evidence.

## Publication decision — exact currently supported claim

**Computational GO:** transparently reproducible candidate prioritization combining a 57-donor discovery cohort, 2 additional independent published TM experiments, study-context-dependent gene response, and sensitivity tests in the same human donors under two cell-annotation algorithms. The **three leading response-associated transcript candidates are MYOC, TGFBI and MMP2**. ACTA2 is context-specific (not continuously signed across the seven-day cohort).

**Scientific novelty limitation:** MYOC glucocorticoid responsiveness, ECM shifts and existing Singh 2026 ATAC/IOP prioritization have already been reported in source papers. Do not present these genes or the Hallmark EMT label as a newly discovered molecular mechanism. A strong thesis will instead focus on systematic independent *directional* replication, adverse selection/confounding audit, transparent transcript evidence tiers, and any additional orthogonal data not already reported.

**HOLD:** true batch-independent pharmacologic effect, inferred 48h-vs-7d temporal mechanism, original publication's 80 named high-effect DEG list, quantitative cross-study effect-size meta-analysis (CPM vs edgeR scale differences), independent ECM protein/outflow assay, and actual steroid-induced delta-IOP human clinical validation. Original GSE313736 Series abstract reports 63 up and 17 down high-effect genes; its official Series supplements do not expose their individual gene list, so exact identity-level 80-vs-83 overlap must not be invented.

## Immediate next research tasks

1. **Manuscript top-level structure**, not more exploratory bulk figure automation: Introduction (clinical steroid-response context and TM), datasets and unit-of-analysis methods, data-provenance and biases, 3-study direction evidence, paired annotation sensitivity, Discussion with limitations.
2. **Original study author 80-gene table access:** contact primary authors or locate legitimate publication supplements when available. Do not guess identity from summary counts.
3. **Novelty gate:** compare each claim to the original sc-eQTL/GSE313736 paper, Mehrotra 2026 bulk TM qPCR, and Singh 2026 ATAC/IOP sources. Add an exact claim-by-claim overlap table before submitting.
4. **Optional mechanistic extension** requires randomized same-donor vehicle/100/500nM at matching 48h and other timepoints, appropriate culture batch balance, qPCR MYOC/TGFBI/MMP2 and MMP2 zymography/secreted ECM assays. If unavailable, frame the work as computational transcriptomic reproducibility only.

The Figures 1–3 in the CI artifact are **ready for use as accurately labeled computational figures**, not complete manuscript proof or an article-ready mechanistic claim.
