# CKM→Stroke GLGC EUR: independently verified technical GWAS progress

Evidence cut: 2026-11-16 (KST). Scope is **technical per-chromosome SNP-only QC**, not completion of scientific P1–P6.

## Observed trait–chromosome evidence

| Chr (GRCh37) | GLGC trait | Source rows | Retained SNPs | Multi-base / non-SNP | Symbolic / IUPAC non-SNP | Duplicate SNP keys | Large sample P–beta/SE mismatch | Technical status |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- | --- |
| 22 | LDL-C | 648,183 | 590,273 | 43,832 | 14,078 | 0 | 0/5,844 | PASS |
| 22 | HDL-C | 636,267 | 592,304 | 43,875 | 88 | 0 | 0/5,864 | PASS |
| 22 | logTG | 651,154 | 593,058 | 44,019 | 14,077 | 0 | 0/5,871 | PASS |
| 21 | LDL-C | 640,388 | 581,991 | 44,910 | 13,487 | 0 | 0/5,762 | PASS |
| 21 | HDL-C | 629,325 | 584,370 | 44,955 | 0 | 0 | 0/5,785 | PASS |
| 21 | logTG | 642,765 | 584,192 | 45,087 | 13,486 | 0 | 0/5,784 | PASS |
| 20 | LDL-C | 1,063,703 | 975,190 | 68,077 | 20,436 | 0 | 0/9,655 | PASS |
| 20 | HDL-C | 1,043,874 | 975,740 | 68,134 | 0 | 0 | 0/9,660 | PASS |
| 20 | logTG | 1,068,738 | 979,949 | 68,356 | 20,433 | 0 | 0/9,702 | PASS |

Raw source rows and retained normalized SNPs reconcile after *explicit* exclusions in each of these 9 technical results. The dataset is GLGC Graham2021 European ancestry-specific LDL, HDL, and log triglycerides, with the ALT allele considered the effect allele and reference build GRCh37 (`POS_b37`).

**Numerator of independently verified per-trait real-data QC: 9/66.** Three chromosomes (22/21/20) have 3/3 real-data trait technical PASS **individually**. Only chr22 has already completed its official independent chromosome-level aggregate. chr21 and chr20 are **pending independent metadata aggregation**, so the formally consolidated 3/22 result is **not yet a completed CI PASS**.

Subtotal over these nine chromosome/trait evidence records: **6,457,067 retained SNP rows** (these are *per phenotype*, not unique variants merged across phenotypes). None of this validates the entire 22-autosome data or original full-source SHA.

## Evidence sources and software-version provenance

- **chr22 full-chromosome run:** [#43431051168](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43431051168) complete, all 5 jobs PASS; original source QC commit `91de63c4916a4f078ca5490b35c530b1b1709b31`. Its `ckm-glgc-chr-aggregate` confirms chromosome technical PASS with `scientific_status=HOLD`.
- **chr21 and chr20:** [#43848671279](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43848671279) complete **per-trait for these two chromosomes**, source QC commit `fea2ce570e5abcad7f1dc80e5287e8639baed514`. chr19 remains queued. Nine-job run and its three aggregate jobs are incomplete.
- **Coverage snapshot (metadata only):** [#43857732975](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43857732975) submitted, not yet certified PASS. It downloads existing source evidence from the above runs and independently aggregates chr21 and chr20, asserting `PARTIAL`, 3/22 chromosome technical results, 9/66 trait-chromosome results, 6,457,067 retained SNP entries, and `scientific_status=HOLD`.
- **Coverage logic:** `scripts/ckm_stroke/audit_glgc_autosomes.py` with tests `tests/test_ckm_stroke_audit_glgc_autosomes.py`. [Contract CI #43854658801](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43854658801) is queued at the evidence cut.

## Required before continuing research

1. Complete **chr19** LDL-C, HDL-C, logTG independently, plus the full first-batch three chromosome aggregates. Ensure **9/9** per-trait and **3/3** aggregate PASS.
2. Run the next gated batch **chr18,17,16** only after all predecessor evidence is validated; do not pile up multiple large-source CI runs during GitHub Actions queue congestion.
3. Continue until all 22 autosomes have explicit technical source metadata, with per-phenotype totals, no duplicates, and harmonization-ready effect-allele/frequency fields. Existing scratch normalized SNP files were removed after each run; summary-only metadata is insufficient to perform downstream analyses.
4. Before P1 science PASS or P2 LDSC / MR, establish independent original source integrity and full GWAS data availability, allele harmonization, LD reference ancestry/build, cross-cohort sample-overlap checks, and appropriate quantitative phenotype/statistical-scale handling.

**Interpretation:** This supports that specific GLGC chromosome input segments and their SNP-only parsers function correctly. It does **not** support CKM→stroke causal inference, cross-trait colocalization, or prioritizing candidates.

Tracking: [Issue #285](https://github.com/kimtk94/IS_Analysis_V3/issues/285) · [Draft PR #284](https://github.com/kimtk94/IS_Analysis_V3/pull/284).

## chr19 extreme-significance source P-value underflow (new 2026-11-16 evidence)

Original GLGC chr19 jobs from [#43848671279](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43848671279) show source `pvalue == 0` for **LDL-C 332**, **HDL-C 8**, **logTG 29** rows (total **369**); these rows were *excluded* from the reported valid SNP-only counts. Those are not scientific negative associations. They require independent recovery/verification against original `pvalue_neg_log10` for subsequent instrument discovery and must not silently disappear in downstream MR input. Excluding these records still permits a deliberately scoped **technical subset QC PASS**, but all scientific and source-SNP-export readiness remains HOLD. New per-chromosome summaries expose `excluded_p_underflow_zero` explicitly; older pre-fix aggregates without this field must be treated as **unreported**, not assumed zero.

The existing sampled β/SE-to-P diagnostic performs comparisons only where both P values avoid severe floating-point underflow; the number of sampled rows is not necessarily the number of numerically comparable rows.
