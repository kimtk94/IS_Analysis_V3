# IS 646-assay Allele-Frequency, H3/H4 Prior and 2,225-gene Integration Audit

**Analysis date:** 2026-10-11 (Asia/Seoul)  
**Data policy:** existing BBJ ischemic stroke + GTEx v8 SNP-level input backup only. **No CKB, TPMI, new controlled-access, paid or externally gated GWAS database used.** Existing canonical data and research code checkout untouched; all new outputs in timestamped `results/is/audits/` subdirectories.  
**Scope:** descriptive scientific validation of the already archived 646 coloc ABF gene–tissue assays; NOT 646 independent lead loci, NOT 646 independent molecular mechanisms and NOT 646 validated causal genes.

## 1. Full input provenance

- Source files: **646/646 files retrieved from the user's already-existing Google Drive `MASTER_DEGREE/IS_COLAB/coloc_ready_v2/abf` backup** to an immutable derivative directory.
- `rclone check --checksum` found **646 files matching, zero differences** against the Drive originals.
- **646/646 archived coloc H0–H4 results reproduced** using direct SNP-level `coloc.abf` version 5.2.3 at `p1=p2=1e-4, p12=1e-5`. Maximum absolute posterior error **6.71685×10⁻¹⁵**. This is numerical reproducibility, not mechanism validation.
- Matched GTEx molecular-QTL LD was not supplied, and the coloc ABF **single-shared-causal-variant modeling assumption is not validated**. (ABF itself does not require an LD matrix; multi-signal conditional methods do.)

## 2. 1,224,224 input rows correspond to only 8,956 unique locus–SNP identities

Full 646 files contain **1,224,224 SNP–gene–tissue rows**, across just four existing BBJ IS regions. Many SNPs repeat across different genes and GTEx tissues; this multiplicity **must not** be read as 1.2 million independent variants.

| BBJ IS locus | Unique GRCh38 SNP identities | Unique palindromic SNPs | Palindromic sites near AF=0.5 in GWAS or GTEx |
|---|---:|---:|---:|
| L001 | 1,694 | 299 | 52 |
| L002 | 1,583 | 216 | 53 |
| L003 | 2,663 | 341 | 46 |
| L004 | 3,016 | 387 | 84 |
| **Total** | **8,956** | **1,243** | **235** |

Palindromic denotes the `A/T, T/A, C/G, G/C` allele pairs. Frequency near 0.5 was flagged if the original GWAS EAF or any sampled GTEx `AC/AN` fell in the approximate 0.42–0.58 window. These are **audit flags**, not automatic exclusions or proof of strand errors.

## 3. MAF versus effect-allele frequency: distinction is material

Mathematical consistency checks across all 1,224,224 local records:

- `gwas_maf = min(gwas_eaf, 1-gwas_eaf)`: **PASS** to 1e-9 tolerance.
- `eqtl_maf = min(ac/an, 1-ac/an)`: **PASS** to 1e-4 tolerance.
- All 1,224,224 input rows have GTEx `AC/AN <= 0.5` (0 rows above 0.5). However **423,138** rows have `gwas_eaf > 0.5`.
- Absolute minor-frequency difference `|gwas_maf-eqtl_maf|>0.1` is seen in **508,318/1,224,224** SNP–assay rows (~41.5%). This is *not* an independent-SNP rate, and 646 gene/tissue calculations reuse only 8,956 locus–variant identities.
- **646/646 assays** have >25% of their SNP rows showing this minor-frequency difference; **58/646** have >50%.
- More stringent `>0.2` MAF difference appears in **199,221** SNP–assay rows.

**What the source documentation actually says:** eQTL Catalogue documents `maf` as minor allele frequency, `ac` as **alternative allele count**, and `alt` as the effect allele; GTEx states eQTL normalized effect direction is ALT relative to reference. Citations:
- https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md
- https://www.ebi.ac.uk/eqtl/Data_access/
- https://www.gtexportal.org/home/faq

The observation that `AC/AN<=0.5` for every archived row is unusual and demands **local/source-processing provenance review**. It does **not** mathematically demonstrate that `ac` is *minor* allele count instead of ALT allele count, because a selected variant set could conceivably have only ALT-minor SNPs. The exact-reference coordinate/REF/ALT match alone does not independently attest original **BBJ GWAS effect-allele orientation**, and comparison of GWAS EAF to an unoriented MAF must **never be used to automatically flip effects**.

**State:** SNP reference-coordinate matching PASS; EAF/MAF numerical identities PASS; cross-study effect-allele orientation SOURCE_METADATA_CONFIRMATION_PENDING; automatic SNP removal / effect-sign flip NOT PERFORMED. Ancestry-specific genetic frequency differences between BBJ Japanese/EAS and GTEx donor populations may be substantial without representing data errors.

## 4. Direct 646 × 5 coloc prior sensitivity with H3/H4 comparison

Using the same complete per-SNP beta, SE and mathematical MAF as the archived notebook, re-ran `coloc.abf` for **3,230 models**. The p12 grid was defined prior to review and retained in full:

| p12 | PP.H4 ≥ 0.5 | PP.H4 ≥ 0.8 | PP.H4 > PP.H3 |
|---|---:|---:|---:|
| 1×10⁻⁶ | 0 | 0 | 1 |
| 3×10⁻⁶ | 2 | 0 | 7 |
| 1×10⁻⁵ **(archived baseline)** | 6 | 0 | 40 |
| 3×10⁻⁵ | 13 | 4 | 435 |
| 1×10⁻⁴ | 86 | 12 | 565 |

**Only one assay has H4 > H3 at all five p12 prior values** (CALHM2 cerebellar hemisphere); **none have H4 ≥ 0.8 across all five values**. The 12 assays with H4≥0.8 at the most favorable high prior represent **seven genes: FGF5, CALHM2, NEURL1 (alias NEURL), C4orf22, INA, SH3PXD2A and COL4A2**. These are all *prior-sensitive exploratory associations*, not causal gene validation.

Illustrative selected gene–tissue assays (cerebellar hemisphere):

| Gene | H4 at 1e-6 | H4 at 1e-5 | H4 at 1e-4 |
|---|---:|---:|---:|
| CALHM2 | 0.28278 | 0.79768 | 0.97526 |
| FGF5 | 0.26028 | 0.77870 | 0.97237 |
| NEURL1 | 0.17366 | 0.67759 | 0.95458 |
| C4orf22 | 0.16252 | 0.65994 | 0.95100 |

**Causal inference from prior inflation is forbidden:** increasing p12 without adding data can raise H4. There are 646 correlated gene–tissue analyses, and no multiplicity correction, matched molecular LD fine mapping or demonstration of single signal per locus. H3/H4 are hypotheses under the model, not experimentally validated molecular mechanisms.

**Technical repair:** initial manuscript-like JSON output rounded p12=3×10⁻⁵ to zero owing to `jsonlite` default `digits=4`. Re-ran with `digits=12`; the complete 3,230-row TSV is bytewise SHA-256 identical between v1 and corrected v2, and only corrected **v2** JSON should be used. Final corrected source:
`/srv/is-analysis/results/is/audits/is_646_direct_p12_5grid_20261011_v2/`

## 5. Scalar-QTL sample-N sensitivity (previously executed and carried forward)

Every one of the 646 GTEx assays contains variable SNP-specific `eqtl_n`, median 8 distinct sample-N values per gene–tissue file. Earlier original ABF model used the first SNP's `qtl_n_scalar` as the single QTL sample-N.

Using four **alternative scalar assumptions** (first, min, median rounded, max; 2,584 model fits):
- Maximum absolute H4 change: **0.0429922895**.
- Median assay-level maximum change: approximately **0.0014**.
- 0/646 changed H4 by >0.05, **1/646 crossed H4=0.5**, none crossed H4=0.8.
- INA, BA24 anterior cingulate cortex: H4 `0.53637` at original first N vs `0.49338` at minimum QTL N.

This is **scalar-N sensitivity**, not a per-SNP sample-N coloc model and not validation of imputation, ancestry effects or multi-signal fine mapping.

## 6. Full positional universe of 2,225 Ensembl IDs preserved

The upstream expanded table contains **2,225 positional Ensembl stable gene IDs** linked to **2,425 gene-region records** across 80 provisional GWAS regions (869 EAS/Japanese, 1,556 EUR).

- Existing GTEx/BBJ coloc work pertains to **43 gene IDs in four BBJ/EAS loci**.
- **2,182 additional candidate IDs** did not have a gene–tissue assay in this 646-pair archive; this is **NOT_TESTED_IN_THIS_BRANCH**, not tested-negative.
- Seven genes have H4≥0.8 in at least one assay **only at higher p12**; none are promoted to causal status.
- Historical `NEURL1` and GENCODE-v19 `NEURL` map to the **same ENSG00000107954**, avoiding a false 2,226th candidate.
- No EUR candidate was assigned a legacy BBJ GTEx coloc just because its symbol matched.
- Final gene overlay:
  `/srv/is-analysis/results/is/audits/is_2225_prior_af_integrated_20261011_v1/`
  - `IS_2225_GENE_ABF_PRIOR_AF_INTEGRATED.tsv`
  - `IS_646_GENE_TISSUE_PRIOR_AF_CONTEXT.tsv`
  - `IS_2225_PRIOR_AF_INTEGRATION_MANIFEST.json`

## 7. Figure assets

All figures generated directly from full TSV evidence, without statistical selection, using base R and preserved as editable SVG in GitHub:

1. [646 assay p12 H3/H4 counts](../../docs/figures/is/IS_646_DIRECT_P12_THRESHOLD_COUNTS.svg)
2. [Eight prespecified gene–tissue H4 prior curves](../../docs/figures/is/IS_PRIORITY8_DIRECT_P12_H4_SENSITIVITY.svg)
3. [2,225-gene QTL branch coverage](../../docs/figures/is/IS_2225_GENE_QTL_COVERAGE.svg)

Source `scripts/is/render_is_646_coloc_evidence_figures.R`; corresponding server derivative `results/is/audits/is_646_2225_publication_figures_20261011_v1/`.

## 8. Versioned code and science gates

Sources committed on `research/is-broad-discovery-20261009`:
- `scripts/is/audit_646_coloc_allele_frequency_provenance.py` (+5 offline tests).
- `scripts/is/run_646_coloc_direct_p12_sensitivity.R` (corrected `digits=12`).
- `scripts/is/build_is_2225_prior_af_integrated_context.py` (+4 offline tests).
- `scripts/is/render_is_646_coloc_evidence_figures.R`.
- Existing 646 replay, QTL-N scalar sensitivity and original 2,225 matrix scripts remain pinned unchanged.

All new input-specific scientific validations were performed on the personal server. The CI workflow was extended to cover new Python sources/tests. **GitHub Actions remote CI completion is separate and is not claimed.**

**Scientific states:**
- `TECHNICAL_REPRODUCTION_646/646 = COMPLETE`.
- `ALLELE_FREQ_NUMERIC_QC = COMPLETE`.
- `PALINDROMIC_ALLELE_ORIENTATION = SOURCE_METADATA_REVIEW`.
- `GTEx_COHORT_MATCHED_LD = NOT_VERIFIED`.
- `MULTI_SIGNAL_COLOCALIZATION_WITH_MATCHED_LD = PENDING`.
- `CAUSAL_GENE_ASSIGNMENT = NOT_ESTABLISHED`.
- Alcohol-related **ALDH2 BLOCKED**, **ADH1B EXPLORATORY**; no status promotion.
- `POSITIONAL_CANDIDATE_2225 = RETAIN_ALL`.

### Next tasks requiring no new DB

1. Review BBJ GWAS original effect-allele column and GTEx eQTL Catalogue `AC` preprocessing lineage in the **already saved local source/docs**; do not infer flips from MAF.
2. Stratify MAF variance and prior-sensitivity by locus, tissue, palindromic flags and allele frequency, and demonstrate the effect of *pre-specified* QC exclusions on n and H4 without retroactive H4-focused selection.
3. Extend healthy brain reference donor-level cell evidence carefully, retaining 6-donor and feature-not-present caveats. The sample is healthy tissue, not stroke patient vs control.
4. Update IS manuscript Figures/Results and dashboard with corrected full-prior narrative and robustly marked UNTESTED/BLOCKED states.
