# IS P3A: correct coloc case-control input contract (binary beta+SE does NOT require s)

**2026-11-23 KST**. Research: original ischemic stroke candidate locus G0015, 36 gene×QTL tissue groups, 15 ROSMAP genes. Source full chr10 original GIGASTROKE EAS/EUR effect-size units have already been verified for **725,707 source SNP rows** in [PR #325](https://github.com/kimtk94/IS_Analysis_V3/pull/325). All formal coloc analyses remain **0/36 HOLD** until separate source-quality and harmonization checks.

## Methodological correction from original coloc implementation (IMPORTANT)

We previously treated absence of GIGASTROKE file-specific case fraction s as an unconditional blocker to coloc.abf, **even though both original GIGASTROKE files supply beta (log odds ratio) and SE**. This is too restrictive. The actual original coloc R code has distinct pathways:

1. **Case/control direct effect data**: list(type="cc", snp=unique_allele_harmonized_IDs, beta=log_OR, varbeta=SE^2). When beta and varbeta are present, original check_dataset accepts this input WITHOUT s and WITHOUT N. The original process.dataset computes z=beta/sqrt(varbeta), and approx.bf.estimates selects the case/control log-OR effect prior. Supplying different hypothetical s values must not alter these approximate Bayes factors.
2. **Case/control p-value-only fallback**: list(type="cc", snp=IDs, pvalues, MAF, N, s). If beta or varbeta is missing, original check_dataset explicitly requires s (fraction cases), N and MAF. This fallback is unnecessary for current original GWAS files and should **not** be invoked merely because source original publisher case fractions are missing.
3. **Quantitative QTL**: list(type="quant", snp=IDs, beta, varbeta=SE^2, sdY=measured_SD), OR if sdY missing, MAF and a justified scalar N for estimating sdY. ROSMAP original QTL 15-gene source AN/2=560 and empirical R-equivalent sdY estimates were independently calculated under that scenario in [PR #326](https://github.com/kimtk94/IS_Analysis_V3/pull/326); real model N and normalized-expression SD still unapproved. That issue is NOT removed by the binary GWAS correction.

Source authority: original [coloc check_dataset implementation](https://github.com/chr1swallace/coloc/blob/27776b1bb14e19d1542798ac7477fa378fc2fc0d/R/check.R) at pinned SHA1 Git blob **7d244c2826ed36f363b986a26e4e0c612de9af1e**, original [process.dataset / approx.bf.estimates](https://github.com/chr1swallace/coloc/blob/27776b1bb14e19d1542798ac7477fa378fc2fc0d/R/claudia.R) Git blob **043783ac5ddec3477ff3b91e3e1faf271c51ce09**. These functions come from upstream immutable commit **27776b1bb14e19d1542798ac7477fa378fc2fc0d**. Source Git blobs are independently SHA1 validated at runtime. Original CRAN vignette [Coloc data structures](https://chr1swallace.github.io/coloc/articles/a02_data.html) recommends beta/varbeta over p-value approximations. Some narrative documentation lists s for cc without distinguishing input paths; the source code AND original R live tests settle this specific ambiguity.

## Native original R execution in GitHub Actions

The workflow fetches exactly two pinned source R files (under 40KB total) and directly executes original source functions with base R, **not an independent reimplementation**. Actual original R tests:
- check_dataset on type cc beta+varbeta no s -> accepted
- process.dataset original R on type cc beta+varbeta no s -> finite ABFs
- adding s=0.07 vs s=0.57 leaves original beta/SE-based ABFs identical
- p-value-only case/control with N/MAF but missing s -> explicitly rejected
- p-value-only case/control WITH s -> accepted
- quantitative beta/SE with known sdY accepted, while unknown sdY without N+MAF rejected.

No external patient data, GWAS SNP effect rows, or the full original GWAS download are required for this *methodological contract verification*. No formal coloc posterior H4 is generated or claimed.

## Remaining scientific requirements — corrected priority

1. Remove **case_fraction_missing_alone** from the hard stop list **only for an original GWAS beta+SE type cc model**; never impute a value or silently switch to p-value-only inputs. Still document case/control composition and ascertainment heterogeneity as evidence provenance and for comparisons to other models.
2. Verify quantitative eQTL source model N and trait scaling used in the 15 ROSMAP estimates (AN/2=560 still a provisional hypothesis) and the remaining 21 GTEx tissues. Validate model-derived sdY against original expression normalization or defensible study-specific method evidence.
3. Verify complete matched cis QTL/GWAS variants, liftover GRCh37↔GRCh38, allele flips, INFO/MAF, duplicates/multi-allelic, missingness, physical locus provenance, ancestry-matched LD where needed and prior sensitivity.
4. Coloc.abf single-causal-variant assumption must be explicit, with H3 vs H4 and prior sensitivity tests; coloc.susie requires appropriate shared-variant LD and credible fine-mapping, not simply SuSiE output from one trait. BBJ is included in original GIGASTROKE EAS discovery; do not claim BBJ independent replication.

**Scientific status: 36/36 still HOLD, formal coloc.abf/susie 0 executed.** This PR changes only a false *technical input requirement*, not approval of biological causality.

## Reproducibility

- Python source/hash and R-test orchestrator: scripts/is/audit_p3a_coloc_original_cc_beta_no_s_gate.py
- Base-R original-source runtime tests: tests/is/test_coloc_cc_beta_no_s_original_source.R
- Python negative contract tests: tests/test_is_p3a_coloc_original_cc_beta_no_s_gate.py
- GitHub Actions: .github/workflows/is-p3a-coloc-original-cc-beta-no-s.yml
- Aggregate evidence: artifacts/is/IS_P3A_COLOC_Original_CC_Beta_SE_S_NotRequired_HOLD.json

Do not label absence of s as a source failure or approval of exact meta-analysis per-SNP N. Do not supply invented case fraction, effective N or QTL sdY=1.
