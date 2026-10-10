# IS 30 reproduced ABF assays: joined SNP-input, scalar QTL-N, and p12 sensitivity

**Audit date:** 2026-10-11 KST. **Scientific gate:** technical reproducibility and assumption sensitivity only; no causal gene or independent colocalization validated.

The source audit used 30 numerically replayed ABF assays from 14 distinct Ensembl gene IDs (out of 646 archived total, 616 source-input missing). Each assay was joined by exact locus × dataset × gene ID to the already computed 646-assay scalar-N sensitivity and 646-assay archived posterior p12-reweighting grid; baseline original PP.H4 was required to agree within 1e-8 in all three inputs.

| Property | Direct joined result |
|---|---:|
| Posterior p12 range | 1e-6, 3e-6, 1e-5, 3e-5, 1e-4 |
| Assays whose reweighted H4 crosses threshold 0.5 | **18/30** |
| Assays whose reweighted H4 crosses threshold 0.8 | **12/30** |
| Assays whose scalar QTL-N range crosses H4 threshold 0.5 | **1/30** |
| Maximum H4 span due to p12 reweighting | **0.81675** |
| Maximum abs H4 shift across scalar-N models | **0.04299** |
| Assays with QTL sample N varying across SNPs | **30/30** |
| Assays where over 25% SNPs have cross-study MAF difference >0.1 | **30/30** |
| Source assays still unavailable | 616 |

Top existing baseline PP.H4 assays:
- CALHM2, cerebellar hemisphere: baseline **0.79768**; p12-grid **0.28278–0.97527**; scalar N **0.76664–0.79768**.
- FGF5, cerebellar hemisphere: baseline **0.77870**; p12-grid **0.26027–0.97239**; scalar N **0.77703–0.77870**.
- NEURL, cerebellar hemisphere: baseline **0.67759**; p12-grid **0.17366–0.95458**; scalar N **0.66523–0.67759**.
- C4orf22, cerebellar hemisphere: baseline **0.65994**; p12-grid **0.16250–0.95106**; scalar N **0.65555–0.65994**.

**Interpretation limits:** p12 scenario uses original archived posterior reweighting under assumed original p1=1e-4, p2=1e-4, p12=1e-5, not an independent full-SNP rerun under new priors. Scalar-N first/min/median/max scenarios are alternative *scalar* N assumptions; none substitutes variant-specific sample N. Discordant GWAS-vs-QTL MAFs may reflect East-Asian vs predominantly European ancestry rather than sample errors; donor-matched LD is unverified. ABF single-signal assumptions, tissue screening multiplicity and assay-level independence remain unverified. Thus H4>0.8 under a permissive prior must not be reported as validated causal shared signal.

**New reproducible code:** scripts/is/audit_is_30_abf_joint_sensitivity.py and tests/test_is_30_abf_joint_sensitivity.py.

**Outputs:** /srv/is-analysis/results/is/audits/is_30_abf_joint_sensitivity_20261011_v1/IS_30_ABF_JOINT_SOURCE_N_AND_PRIOR_SENSITIVITY.tsv and IS_30_ABF_JOINT_SENSITIVITY_SUMMARY.json, with hashes for three read-only source tables.

**Next:** source-defined single-SNP N should be incorporated into an independently validated coloc input computation only if algorithm and harmonization assumptions can support it; else prioritize tissue-aware matched ancestry and full tested variant coverage. Preserve entire 2,225-gene search universe.
