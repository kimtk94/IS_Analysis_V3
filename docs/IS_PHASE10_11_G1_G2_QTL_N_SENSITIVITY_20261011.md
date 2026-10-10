# IS Phase10/11 G1 structural alleles and G2 SNP-varying QTL N sensitivity

**2026-10-11; 30 of the archived 646 original gene–tissue SNP input TSV files evaluated on server; no output or raw source overwritten.**

## Core finding — historical scalar N is the per-test maximum

- **646/646 original historical input-index rows** have qtl_n_first = qtl_n_max. Their per-SNP sample-size minimum ranges from **1.95% to 20.63%** below the corresponding maximum.
- The original 646/646 SNP-level H0–H4 Colab rerun was computationally completed previously; this new robustness analysis is **not** that same reproduction test.
- N stress scenarios preserve the same per-SNP beta, standard error, MAF, matched SNPs, and original coloc 5.2.3 priors (p1/p2=1e-4, p12=1e-5); only the single scalar QTL N in the quantitative-trait coloc dataset is varied: FIRST/MAX, MEDIAN, MIN.
- The first/MAX scenario reproduces the original historical H3/H4 for each of the 30 cached sources to absolute tolerance 1e-8. In this 30-source cache, MEDIAN = FIRST/MAX for all inputs.
- A scalar MIN is an **extreme bound-type stress scenario, not a correct or validated SNP-specific effective sample size model**. Do not replace historical posterior, reinterpret H4 as calibrated probability, or claim a gene becomes causal because a threshold is preserved.

## Numeric results on the 30 original source assays

| Gene / original GTEx tissue | QTL N: minimum / original maximum | Original PP.H4 | MIN-N PP.H4 | Change |
|---|---|---:|---:|---:|
| FGF5 — Brain Cerebellar Hemisphere | 170 / 175 | 0.7787 | 0.7770 | -0.0017 |
| C4orf22 — Brain Cerebellum | 204 / 209 | 0.5597 | 0.5557 | -0.0040 |
| INA — Brain Anterior cingulate cortex BA24 | 120 / 147 | 0.5364 | 0.4934 | -0.0430 |
| NEURL1 — Brain Cerebellar Hemisphere | 150 / 175 | 0.6776 | 0.6652 | -0.0124 |
| CALHM2 — Brain Cerebellar Hemisphere | 150 / 175 | 0.7977 | 0.7666 | -0.0310 |
| COL4A2 — Brain Amygdala | 112 / 129 | 0.3063 | 0.2857 | -0.0206 |

- **30 PASS / 0 FAIL / 616 NOT_CACHED** for scalar QTL N stress; the denominator remains all 646 original source tests.
- Original PP.H4 >= 0.50: **6/30**; minimum-N PP.H4 >= 0.50: **5/30**; one threshold crossing is **INA in anterior cingulate cortex BA24 (0.5364 -> 0.4934)**.
- Original PP.H4 >= 0.75: **2/30**; minimum-N PP.H4 >= 0.75: **2/30**. This is NOT a claim that both would remain above threshold under an actual genotype-level per-SNP-N statistical model.
- Maximum absolute H4 change = **0.0429923**, median absolute change = **0.00290378** (30 source assays). The sign and magnitude vary across tests.

## G1 internal structural allele audit

- **30 PASS_STRUCTURAL_ONLY / 0 FAIL / 616 NOT_CACHED**, representing **53,979 gene–SNP records**, NOT 53,979 unique genetic variants.
- All checked cached source QTL GRCh38 variant IDs agree with exact chromosome/position/REF/ALT, variant_id_hg38 and GTEx variant text; GTEx raw beta/SE/MAF agree with derived eQTL columns. GWAS MAF equals min(GWAS EAF, 1-GWAS EAF).
- **8,357 palindromic A/T or C/G records**, of which **730** have GWAS EAF in [0.4,0.6]. For population-different signals this is a risk/review flag, not a proven orientation error.
- **24,250 gene–SNP records** have |GWAS MAF - GTEx MAF| >0.1 (9,819 >0.2); this difference may reflect BBJ/GTEx ancestry, statistical properties or source differences and cannot independently prove an effect-allele flip.
- Existing harmonization string EXACT_REF_ALT_GRCH38 proves only internal tagged REF/ALT encoding consistency in available sources. A BBJ GWAS native effect-allele codebook and independent GRCh37->GRCh38 mapping are still required for true effect direction verification.

## Full-646 one-cell Colab (prepared, NOT executed for new sensitivity analyses)

[Open 646-input N stress + G1 internal allele source audit](https://colab.research.google.com/github/kimtk94/IS_Analysis_V3/blob/feat/is-p0-evidence-20261010/notebooks/is/IS_PHASE10_11_646_QTL_N_SENSITIVITY_ONECELL.ipynb)

- Reads the original 646 historical Drive SNP TSVs directly. Pinned original 646 index SHA256 and raw master SHA256, coloc 5.2.3, parsed exact embedded R source, no existing-result overwrite.
- Saves to **MASTER_DEGREE/IS_COLAB/results/IS_PHASE10_11_646_QTL_N_SCALAR_SENSITIVITY_V1** (separate from prior original full-646 reproducibility outputs).
- Outputs include all original assay statuses, N scalar sensitivity summary, source SHA256 per file, complete log, G1 structural allele status (another 646 input tests), G1 summary and full manifest with both code SHA256.
- Interpret `FULL_646_SCALAR_N_STRESS_COMPLETE` only if all original files are PASS. `G1_STRUCTURAL_646_PASS_EXTERNAL_EFFECT_ALLELE_UNVERIFIED` would still not prove independent GWAS effect direction.

## Publication gate decisions

- G0 original SNP numerical reproducibility: **COMPLETE 646/646**, unchanged.
- G1 exact internal source ID/beta/MAF invariants: **30/646 current server complete**; G1 external GWAS native effect allele, GRCh37/38 liftover and ancestry validation remain **PENDING**.
- G2 scalar QTL sample-size sensitivity: **30/646 current server complete**; validated per-SNP N modeling, full 646 robustness and GTEx molecular-cohort LD remain **PENDING**.
- G3 independent cerebrovascular stroke cell context/mechanistic replication: **PENDING**.
- Expanded 80 provisional GWAS components and 2,225 positional gene IDs are separate from original four BBJ loci. No gene causality promotion.

## Reproducibility

- R analysis: `scripts/is/is_646_qtl_scalar_n_sensitivity.R`
- Python allele identity audit: `scripts/is/audit_is_646_g1_source_alleles.py`
- Colab generator: `scripts/is/build_is_646_qtl_n_colab.py`
- Source-only numeric snapshot: `results/is/audits/is_646_n_scalar_stress_cache30_20261011_v1` (outside Git)
- Source-only allele audit snapshot: `results/is/audits/is_g1_source_structural_cache30_20261011_v1` (outside Git)
- Existing full-646 H0–H4 Colab and audited full source proof preserved without modification.
