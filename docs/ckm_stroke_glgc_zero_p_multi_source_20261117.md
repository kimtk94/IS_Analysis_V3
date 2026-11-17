# GLGC EUR source zero-P original -log10P numerical validation (2026-11-17)

**Status: technical source-statistics PASS for 1,204 additional variants; scientific HOLD.**

## Original verified source baseline

- All-22 technical QC: [Actions #43871252114](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43871252114), 22/22 chromosomes, 66/66 trait×chromosome technical PASS, with 128,925,825 retained LDL-C / HDL-C / logTG SNP association rows across phenotypes.
- Earlier zero-P original chr19 LDL-C verified 332/332 from [source #43872372345](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43872372345) and [Wald #43872824243](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43872824243), maximum numerical difference 5.457 × 10^-12 log10 units.
- Overall full-autosome audit reported 2,632 excluded other/zero-P **candidate** records. Twelve early trait-chromosome aggregates did not have a dedicated underflow field, so unresolved "other" records must not automatically be called source-zero-P until checked.

## Independent four-case source/Wald pilot

[New source/Wald Actions #43990911870](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43990911870) **SUCCESS**, all 4/4 jobs PASS, each with **12/12 offline tests PASS**. Each job downloaded and verified its previous original source-QC artifact (original source SHA, chromosome, trait, exact original-row and zero-P counts) before re-reading official public GLGC BGZF/TBI chromosome rows. Two jobs ran in parallel at most.

| Chromosome × lipid trait | Original association rows | Original P=0 canonical SNPs | Original -log10P valid | Max |delta log10P| | Duplicate P=0 SNP keys |
| --- | ---: | ---: | ---: | ---: | ---: |
| chr8 HDL-C | 2,603,353 | 286 | 286/286 | 5.0022 × 10^-12 | 0 |
| chr8 logTG | 2,653,852 | 341 | 341/341 | 5.2296 × 10^-12 | 0 |
| chr2 LDL-C | 4,039,086 | 378 | 378/378 | 6.8212 × 10^-13 | 0 |
| chr2 logTG | 4,055,066 | 199 | 199/199 | 3.4106 × 10^-12 | 0 |
| **Total** | **13,351,357** | **1,204** | **1,204/1,204** | **5.2296 × 10^-12** | **0** |

All four cases have **zero absolute source–Wald differences exceeding 1 log10 unit**. The original zero-P source evidence and per-variant diagnostics are retained in four GitHub Actions artifacts for 30 days:

- ckm-glgc-pzero-chr8-HDL_C-evidence
- ckm-glgc-pzero-chr8-TG-evidence
- ckm-glgc-pzero-chr2-LDL_C-evidence
- ckm-glgc-pzero-chr2-TG-evidence

All files are sourced from publicly available GWAS summary statistics, not individual participant data.

## Interpretation and remaining conditions

- **Verified 1,536 source zero-P rows in five chromosome×trait cases** (332 previous chr19 LDL-C + 1,204 new). Against 2,632 all-autosome "other/zero-P" candidates, this is 58.4% verified. **Remaining 1,096 candidate rows require source-level checking**; classification as zero-P is unconfirmed until original source read.
- Internal source -log10(P) and beta/SE-derived Wald log-tail concordance is **NOT independent GWAS replication**; they may be calculated from the same association estimates. No source whole-file SHA, no stable full normalized SNP input, no harmonization/LD matching and no MR causal inference have been established. Overall scientific_status=HOLD.
- Next GO criterion: reconcile remaining potential zero-P records by exact trait/chromosome and original source row, preserve a float-safe NEG_LOG10P field with SOURCE_P_UNDERFLOW flag, and retain full original N/EAF/rsID/REF/ALT/BETA/SE fields before cross-study build, effect orientation, LD score, sample-overlap and stroke-subtype sensitivity analyses.
- Updated glgc_zero_p_logdomain_wald.py now enforces chromosome identity instead of emitting a misleading hard-coded chr19 diagnostic for all cases.

**Links:** [Issue #285](https://github.com/kimtk94/IS_Analysis_V3/issues/285) · [Draft PR #284](https://github.com/kimtk94/IS_Analysis_V3/pull/284).
