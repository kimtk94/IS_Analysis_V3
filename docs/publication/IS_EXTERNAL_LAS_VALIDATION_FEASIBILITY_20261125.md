# External LAS case-control feasibility — planning, not observed replication

Use the fixed 24-marker family, exact rs2767060 C/T T allele and same phenotype only if the consortium supplies stroke-free controls, TOAST-defined LAS, model-compatible effect estimates and cohort non-overlap. Freeze the directional two-sided Bonferroni-24 endpoint before external readout. This evaluates the EAS LAS association; it does NOT independently reproduce the EAS–EUR difference, which requires an independent comparator with harmonized definitions and appropriate covariance. Case-only registry data require a separate protocol and cannot be substituted.

Illustrative asymptotic planning: independent unrelated cases/controls, additive log-odds, HWE genotype variance 2p(1-p), k controls per case. n_cases ≈ (z_(1-alpha/2)+z_power)^2 / [2p(1-p) beta^2 k/(1+k)]. This is a local-effect approximation, not an exact power calculator; genotype uncertainty, covariates, ascertainment and unequal missingness may raise requirements. Use hypothetical external control T frequency equal to the discovery reported EAF 0.0596 solely for planning, and discovery beta 0.3582 conditional on a comparable log-odds scale. 24-marker alpha=.05/24. All case counts rounded up to the next 100.

| Effect retained | Hypothetical beta | Approx OR if log-odds | Cases (80%, 4 controls/case) | Controls |
|---|---:|---:|---:|---:|
| 100% | 0.35820 | 1.431 | 1400 | 5600 |
| 75% | 0.26865 | 1.308 | 2400 | 9600 |
| 50% | 0.17910 | 1.196 | 5400 | 21600 |

See [36 planning scenarios](evidence/is_external_las_power_scenarios_20261125.json) spanning three frequencies, three effect sizes, two ratios and two power targets. Do not treat any scenario as an observed sample size or actual power. The GENESIS-K abstract reports 5,328 total participants, not 5,328 LAS cases, and stroke-free controls are unverified. The Taipei 716 clinical cases are not verified LAS counts or controls. TPMI deidentified subtype GWAS unavailable as previously documented. Request aggregate 24-marker effects/SE/P/allele frequencies, imputation and effective per-marker N, cases and controls by subtype, recruitment/overlap, ascertainment and model covariates via approved access. No request sent, no private data accessed, no new external result, no main merge.
