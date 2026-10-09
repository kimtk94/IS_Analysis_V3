# Ischemic stroke G0022 (chr12) multi-signal SuSiE-RSS PILOT
Date: 2026-10-09 KST
Research branch: `research/is-broad-discovery-20261009`

## Study purpose and firm limitations

The prior real-data PLINK2 EAS clumping screen found **3 AIS and 2 AS clump indices** in the large chromosome 12 G0022 discovery component. We tested whether local GWAS Z scores and signed EAS LD are technically compatible with a multi-signal SuSiE-RSS pilot. **We are NOT declaring 5 independent causal signals, 5 validated credible sets, or causal genes.**

GWAS datasets: GIGASTROKE EAS AIS `GCST90104545` and EAS all-stroke `GCST90104544`. Those studies may share participants and are not independent replications.

LD source: 1000 Genomes Phase 3 EAS, **504 subjects**, GRCh37, G0022 expanded regional stable-ID PGEN reference (101,913 source SNPs).

For each clump index, we extracted ±250kb from the 1KG genotype panel and joined SNPs to the prior 42-pair GWAS/reference-matched TSV. Enforced ALT-coded genotype dosage consistent with GWAS ALT-oriented effect estimate beta, dropped unmatched allele pairs and palindromic SNPs, checked missingness, reference MAF ≥1%, and absolute study/reference ALT EAF deviation ≤0.15. Created Z score vector and **signed** SNP correlation LD. Dosage export and LD were regenerated in full double precision after an initial 8-decimal LD export caused artificial negative eigenvalue warnings. Full-precision files eliminated that warning without materially changing diagnostics.

Hypothetical N values `20,000` and `256,274` were used for mathematical sensitivity; the second is the published total-sample count for EAS AIS and **not verified as the effective sample size of this binary trait**. Neither N should be cited as a validated case/control effective sample size. SuSiE v0.14.2, `L=5`, up to 150 IBSS iterations, requested 95% credible sets. Credible set **purity filtering was explicitly run with `Xcorr=R`, `min_abs_corr=0.5`**; the earlier version without Xcorr could not enforce purity and its unfiltered credible-set counts must not be used.

## Technical QC and pilot inference

One EAS 1KG-derived signed LD matrix per 5 local windows, with 504 reference subjects:

| Study | Window anchored at SNP position | Harmonized pilot SNPs | Both N scenarios converge? | LD inconsistency `s` (max across N) | Purity-filtered pilot CS count per N | Gate |
|---|---:|---:|---|---:|---:|---|
| EAS AIS | 112241766 | 514 | Yes | 0.00929 | 1 / 1 | **Exploratory**, effective N and full locus missing |
| EAS AIS | 113031474 | 561 | Yes | 0.01056 | 1 / 1 | **Exploratory**, effective N and full locus missing |
| EAS AIS | 110675363 | 254 | Yes | 0.00623 | 1 / 1 | **Exploratory**, effective N and full locus missing |
| EAS AS | 111629389 | 459 | **No** | **0.27794** | 3 / 3, **do not interpret** | BLOCKED very high LD mismatch and nonconvergence |
| EAS AS | 112930475 | 529 | Yes | **0.17796** | 3 / 3, **do not interpret** | BLOCKED LD mismatch |

Of 10 SuSiE-RSS runs (5 windows × 2 hypothetical sample sizes), **8 converged**. This is not statistical validation. The independent LD consistency diagnostic `susieR::estimate_s_rss` found **2 high-mismatch AS windows** that must be blocked. The **three AIS windows show lower estimated mismatch** and each produced one *purity-filtered local candidate credible set*, but credible set calibration still requires appropriate effective case/control N, larger GWAS windows, ancestry-specific ref sensitivity, and the complete region analysis.

`estimate_s_rss` documentation: https://stephenslab.github.io/susieR/reference/estimate_s_rss.html ; SuSiE diagnostics: https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html . High `s` indicates GWAS Z/reference LD inconsistency; do not interpret PIP/CS from such a fit as discovery.

## Position-only candidate interpretation

For exploratory AIS windows, the top-PIP variant (under the hypothetical `N=256274` sensitivity) and closest GENCODE v19 GRCh37 gene *body* are:

| AIS window | Top-PIP variant | Pilot top PIP | Nearest gene | Positional relationship |
|---|---|---:|---|---|
| 112241766 | `12:112241766:G:A` | 0.397 | **ALDH2** | Variant within gene body |
| 113031474 | `12:112930475:T:C` | 0.984 | **PTPN11** | Variant within gene body |
| 110675363 | `12:110675363:C:T` | 0.334 | **IFT81** | Variant ~18.8 kb beyond gene body; ATP2A2 ~43.2 kb away |

**A nearby gene is not a causal gene.** Distal regulation and ancestry-specific QTL colocalization may favor a different target. The much larger PIP at the PTPN11-overlapping variant can be an effect of local LD, reference mismatch or window truncation; it must not be promoted to definitive causality or protein-target priority based on PIP alone.

For AS, top-PIP positions had nearest CUX2 and RPL6; **their LD mismatch fails gate**, so they are not nominated as validated genes.

## Reproduce (must run from isolated IS branch)

```bash
cd /srv/is-analysis/worktrees/is-broad-discovery-20261009
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 scripts/is/prepare_g0022_susie_pilot.py
Rscript scripts/is/run_g0022_susie_pilot.R
Rscript scripts/is/diagnose_g0022_susie_ld.R
python3 scripts/is/audit_g0022_susie_pilot_gates.py
python3 scripts/is/annotate_g0022_pilot_signal_genes.py
python3 -m unittest discover -s tests -p 'test_is_*' -q
```

Outputs (not canonical):
`/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/susie_g0022_pilot/`
- Five `*/signed_ld.tsv.gz`, `z.tsv`, `variants.tsv`, `input_qc.json` files
- 10 exploratory SuSiE RDS and SNP-PIP tables
- `G0022_SUSIE_RSS_PILOT_SUMMARY.tsv`
- `G0022_GWAS_LD_CONSISTENCY.tsv`
- `G0022_SUSIE_PILOT_EVIDENCE_GATES.tsv` and summary JSON
- `G0022_PILOT_TOP_SNP_NEAREST_GENES.tsv`

## Publication-readiness gaps (next)

1. Verify EAS AIS and AS per-SNP effective case/control sample sizes. Do not infer from total study N. Identify overlap and harmonization across AS vs AIS.
2. Reconcile AS GWAS Z vs EAS LD: inspect SNP conditional Z outliers (`kriging_rss`), ancestry and strand flips. BLOCKED until fixed.
3. Re-run multi-signal SuSiE on whole merged G0022 region with validated LD panel, enough reference samples and larger locus/multiple causal effects; pilot truncation invalidates complete regional set assertions.
4. Obtain human brain/vascular eQTL and sQTL with locus/Ensembl IDs, fine-mapped QTL signals and true colocalization to decide among ALDH2/PTPN11/IFT81/ATP2A2 and other distal regulatory candidates.
5. Continue independent ancestry EAS/EUR replication without claiming that GIGASTROKE related cohorts are independent validation.
