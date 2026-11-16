# GLGC 2021 EUR: recovering the meaning of source `P=0` SNP records

**Evidence cut 2026-11-16. This is a source-statistics numerical validation, not a causal association or MR readiness result.**

## Problem

GLGC chr19 LDL-C original source SNP-only technical QC detected **332 canonical SNPs** whose original `pvalue` column was numerically `0`. The initial subset QC excluded these rows. This must not be interpreted as `P=1` or as evidence of no association; extreme double-precision tail probabilities can underflow to zero.

The subsequent technical full-autosome audit established **22/22 autosomes, 66/66 lipid trait×autosome QC PASS** and **128,925,825 retained SNP association rows**, but this technical audit alone did not recover original zero-P records.

## Actual official-source evidence

1. **Original GLGC indexed source audit:** [Actions #43872372345](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43872372345) — **SUCCESS / 5 unit tests PASS**. It read the official chr19 LDL-C indexed GWAS original 1,081,460 rows using verified BGZF/TBI and recovered **332/332** canonical SNP rows with `pvalue=0`, valid `N`, ALT-allele AF, rsID, beta, SE and original `pvalue_neg_log10` field. No missing original source logP, duplicate canonical SNP keys, parse errors or row cap. Reported `−log10(P)` range **324.064207720935 to 8410.26925143796**.
2. **Independent log-domain numerical diagnostic on the already downloaded 332-row artifact:** [Actions #43872824243](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43872824243) — **SUCCESS / 6 unit tests PASS**. This run did NOT redownload the full chromosome. It directly compared source `−log10(P)` to beta/SE-derived two-sided Wald Gaussian tail with SciPy `log_ndtr` evaluated in log space.
3. **Concordance:** 332/332 available comparisons; absolute difference in `−log10(P)` median **3.410605131648481 × 10⁻¹³**, p95 **4.547473508864641 × 10⁻¹²**, maximum **5.4569682106375694 × 10⁻¹²**. Number of absolute discrepancies above 1, 5, or 50 log10 units: **all zero**.

Log-domain statistic:

\[
z = \frac{\beta}{SE},\qquad
-\log_{10}(p_{\mathrm{Wald,2s}}) =
-\frac{\log(2) + \mathrm{logNdtr}(-|z|)}{\log(10)}.
\]

This avoids computing `p=2*norm.sf(|z|)` directly in double precision when it would round to zero.

## Interpretation

**Technical numerical result: PASS** for these 332 source rows. Their original minus-log10 P values are internally consistent with the corresponding beta and SE to machine precision. These source `P=0` entries are therefore consistent with **extremely small probabilities suffering ordinary floating-point underflow** rather than absent or biologically negative associations.

This is NOT an independent replication of the GWAS association. Original logP may have been derived from the same underlying beta/SE, so internal agreement does **not** establish independently validated biology, exact model assumptions, original whole-genome file checksum, valid MR instruments, ancestry harmonization or stroke causation.

## Required engineering and scientific gates

- Preserve source `NEG_LOG10P` as a primary association significance field. Do **not** replace large values with `0`, `1`, `NaN`, or an arbitrary conventional floating point `P` without a documented cap and corresponding flag.
- Extend the same source-logP and log-domain β/SE validation to the remaining original `P=0` candidates across all 22 autosomes and three GLGC traits. The all-autosome technical audit's exclusion accounting identifies **2,632** other/zero-P candidate rows, while earlier aggregates have incomplete dedicated underflow fields. Source-log-level independent verification remains outstanding for most candidates.
- Define an explicit complete normalized source-SNP format containing `CHR POS EA NEA BETA SE RSID N EAF NEG_LOG10P` and per-row flags such as `SOURCE_P_UNDERFLOW`. This is a **source representation step**, not harmonization.
- Verify original upstream GWAS source identity/checksums and cross-commit code parity. Confirm transformed trait scales, EUR LD reference ancestry/build, allele matching and ambiguous palindromes before LDSC/MAGMA or MR.
- Keep `scientific_status=HOLD`, and do not claim a completed CKM–stroke causal analysis.

Source code: `scripts/ckm_stroke/glgc_zero_p_probe.py`, `scripts/ckm_stroke/glgc_zero_p_logdomain_wald.py`. Tested input/diagnostic artifacts are attached to the above two GitHub Actions runs.

[All-22 original technical audit](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43871252114) · [Research Issue #285](https://github.com/kimtk94/IS_Analysis_V3/issues/285) · [Draft PR #284](https://github.com/kimtk94/IS_Analysis_V3/pull/284).
