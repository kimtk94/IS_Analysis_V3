**2026-10-11 CORRECTION: The full-646 V1 Colab erroneously expected the 30-source native-audited 8,309 SNP union. Actual full input contains 8,956. V1 failed before per-SNP-N common-sdY R execution. Use [V2 repair with 8,956 GWAS source audit](IS_PHASE10_11_NATIVE_BBJ_GWAS_8956_RECOVERY_20261011.md); previous V1 links below are retained for error reproducibility, not execution.**

# IS Phase10/11 — observed variant N assisted common sdY calibration

**2026-10-11 | Direct local 30/646 SNP-input pilot PASS; full original 646 one-cell Colab BUILT, NOT RUN. This remains statistical robustness analysis, not evidence of a causal gene.**

## Scientific basis

The pinned coloc v5.2.3 internal sdY.est estimates a common quantitative-trait standard deviation using the through-origin regression:

- x_i = 1 / SE_i squared
- y_i = 2 * N_i * MAF_i * (1 - MAF_i)
- sdY = sqrt( sum(x_i * y_i) / sum(x_i squared) )

The historical four-BBJ-locus / 646 GTEx assay pipeline supplied FIRST=MAX QTL N for every SNP. Coloc 5.2.3 sdY.est accepts a vector of per-SNP N, independently confirmed against installed R source and synthetic reference formula. The new method **changes only the estimated common sdY**, supplied explicitly to coloc.abf, while preserving observed beta, SE, SNPs, MAF, GWAS data and p1=p2=1e-4, p12=1e-5. It does NOT correct per-SNP standard errors, estimate genotype-based phenotype variance, fit individual-level models or provide matched GTEx genotype LD.

Official coloc documentation: https://chr1swallace.github.io/coloc/reference/sdY.est.html and https://chr1swallace.github.io/coloc/reference/check_dataset.html

## Actual 30 original cached gene–tissue experiments

- **30 PASS / 0 FAIL / 616 locally MISSING_INPUT** under unique-original-646 denominator. Original historical H3 and H4 reproduce within 1e-8 with first/max N estimated sdY explicitly passed to coloc.
- All 30 gene–tissue input SNPs have original eqtl_an = 2 * eqtl_n. These are reported allele-derived counts, not verified individual sample-level phenotype observation counts.
- **Maximum absolute H4 delta: 0.000046220853571421472**; **median absolute H4 delta: 0.0000057333564075778831**. Maximum absolute relative common sdY shift ~0.01158%.
- All six original PP.H4 >= 0.5 and both PP.H4 >= 0.75 tests remain above their respective cutoffs. No original H4>=0.5 test crosses the descriptive threshold.
- INA anterior cingulate BA24 PP.H4 changes from approximately 0.5363685 to 0.5363223. Under the earlier **all-SNP MIN-N scalar stress**, the same test fell to 0.4933763. The latter is an extreme sensitivity test, NOT a validated variant-specific effective-N correction.
- CALHM2 cerebellar hemisphere approximately 0.7976836 to 0.797648; NEURL1 same tissue approximately 0.677586 to 0.677574. These do NOT establish causal genes.
- Previous JSON serialization precision accidentally rounded a sub-0.0001 maximum difference to zero; current analysis emits digits=16, preserving actual tiny changes.

## New full-646 G1 native-BBJ SNP-catalog proof

Previous native original BBJ Sakaue/Kanai GWAS ZIP plus GRCh37-to-GRCh38 chain audit independently verified 8,309 unique GWAS variants from the 30 cached gene–tissue inputs across four loci: L001=1,694, L002=1,583, L003=2,016, L004=3,016.

All 8,309 native verified records were cross-checked against deterministic columns: locus, variant_id, match_key, gwas_beta, gwas_eaf, making the **expected catalog SHA256**:

d1e237d3c0ad8281adf3d9fe78f69091505eeeab6d152ae0bc54edd5e964853f

The new one-cell Colab reads **all 646 original source SNP TSVs** and reconstructs their exact unique GWAS catalog. The notebook FAILS CLOSED unless the full-646 variant set contains exactly those 8,309 source-verified variants with that same SHA256 and no conflicting beta/EAF/GRCh38 key per variant. **If this check passes, all 646 gene–tissue source assays can be mapped onto the already independently verified native GWAS variant universe.** Until then, native GWAS identity coverage should not be extended beyond the earlier 30 assays.

## Reproducible code and artifacts

- Main runner: scripts/is/is_646_qtl_variant_n_sdy_calibration.R, coloc 5.2.3 pinned, all 646 index rows reported with PASS/FAIL/MISSING.
- Catalog generator: scripts/is/build_expected_gwas_variant_catalog.py.
- One-cell Colab generator: scripts/is/build_is_646_variant_n_sdy_colab.py.
- Full Colab (prepared, NOT run): https://colab.research.google.com/github/kimtk94/IS_Analysis_V3/blob/feat/is-p0-evidence-20261010/notebooks/is/IS_PHASE10_11_646_VARIANT_N_SDY_ONECELL.ipynb
- Versioned independent results location on Google Drive: MyDrive/MASTER_DEGREE/IS_COLAB/results/IS_PHASE10_11_646_VARIANT_N_SDY_CALIBRATION_V1; refuses overwriting any previous results.
- Local 30 source result outside Git: results/is/audits/is_646_variant_n_sdy_cache30_20261011_v2.
- Unit tests test_is_646_variant_n_common_sdy.py include synthetic sdY formula, original posterior reproduction gates, single-cell notebook/source/hash contract and no overwrite.

## Publication release gates

- G0: historical four-BBJ-locus 646/646 original coloc H0–H4 reproducibility COMPLETE.
- G1a: internal SNP REF/ALT and GWAS/eQTL source consistency 646/646 COMPLETE. G1b: original GWAS native Allele2/BETA/EAF plus chain verified for the unique 8,309 SNPs in cached 30 assays; exact 646-assay catalog cross-check BUILT but UNRUN.
- G2a: 646/646 scalar MIN/MEDIAN/FIRST N stress COMPLETE. G2b: per-SNP observed n-assisted common sdY calibration 30/646 actual PASS; 646 source Colab READY but UNRUN. This is not per-SNP SE correction.
- G2c: true genotype/cohort-matched GTEx LD, multi-signal fine mapping, prior and gene/tissue multiplicity PENDING.
- G3: orthogonal stroke-specific cellular regulatory and mechanism evidence PENDING.
- Expanded 80 provisional GWAS regions / 2,225 positional genes are separate from the historical four-locus screen; no target causality promoted.
