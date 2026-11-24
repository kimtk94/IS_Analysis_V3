# IS P3A — BBJ 2021 vs EAS5 source evidence freeze (2026-11-02)

## What was actually run

The **historical**, completed 2026-10-18 GitHub Actions comparison
[run 39124421779](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39124421779),
job 117666449564, independently streamed the Sakaue 2021 public BBJ ischemic
stroke GWAS ZIP and compared SNP-position/allele keys, beta and p-values to the
existing EAS5 archived overlap set using
[scripts/is/compare_bbj_2021_zip_to_eas5.py](../../scripts/is/compare_bbj_2021_zip_to_eas5.py).

The source ZIP ingest succeeded in
[run 39123883155](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39123883155).
Recorded size: 1,611,973,884 bytes. SHA-256:
a8032286cf2206a02c746b138a263930dc7d371c57b94352820005f8a32c8d2b.

The original job's exact numeric JSON has been transcribed into
[the historical evidence record](../../research/is/evidence/bbj2021_eas5_source_comparison_log_20261018.json).
**This freeze is not a fresh comparison of the primary source files.**
The original comparison V1 output did not include the EAS5 archived input
SHA-256; do not invent it retroactively.

## Actual matched subset

| Source-statistics metric | Historical value |
| --- | ---: |
| Official GWAS autosomal source rows streamed | 13,435,541 |
| EAS5 archived candidate variant keys | 25,441 |
| Allele-pair coordinate matches | 25,441 |
| Matched beta and p (either sign) | 25,441 |
| p within V1 log10 tolerance | 25,441 |
| Effect-A2 hypothesis compatible | 25,441 |
| Effect-A1 hypothesis also compatible | 20 |
| A2-only compatible under V1 tolerances | 25,421 |

**Important direction nuance:** The 20 A1-compatible cases are not 20 clear
A1-effect counterexamples. The original beta match criterion used
relative tolerance 0.015 and absolute tolerance 0.00005. For very small betas
both signs can satisfy the numerical criterion. Because every variant satisfies
the A2 hypothesis, those 20 are dual-compatible, leaving 25,421 A2-only
compatible. This is *numerical compatibility*, not official allele-definition proof.

This 100% match is conditional on the **25,441 variant overlap subset**
and the archived variant-key construction; it does not claim that all
13.4 million source rows were present in EAS5 or that all variants are
properly normalized against an independent reference FASTA.

## Cohort metadata available, but not adopted

The official [BBJ PheWeb 2021 download catalog](https://pheweb.jp/downloads)
lists **ischemic stroke: 22,664 cases, 152,022 controls, 174,686 total**.
The [2021 original study](https://www.nature.com/articles/s41588-021-00931-x)
identifies BBJ PheWeb and hum0197 as published summary statistic repositories.
Separately, [BIGA's BBJ harmonization notes](https://www.bigagwas.org/documentation)
say that the original format uses Allele2 as the effect allele and GRCh37
positions. BIGA is a secondary reprocessor, not an original source manifest:
the official release-specific file-format definition, true reference orientation
and case/control analysis model must still be cross-checked before adopting N
or disease definitions into EAS5.

**No permission to use these case/control counts in coloc yet.** BBJ is also
included in the GIGASTROKE East Asian discovery evidence; BBJ–GIGASTROKE is
not independent replication (see
[PR #216](https://github.com/kimtk94/IS_Analysis_V3/pull/216)).

## CI-only reproducibility gate

- [Python fail-closed historical evidence validator](../../scripts/is/audit_bbj_2021_source_freeze.py).
- [Six negative/positive regressions](../../tests/test_is_bbj_2021_source_freeze.py).
- [GitHub-hosted workflow](../../.github/workflows/is-bbj-source-evidence-freeze.yml).
- Aggregate JSON artifact with historical_numeric_concordance_screen=PASS
  and scientific_status=HOLD.

This is strictly **historical evidence integrity**, not fresh primary-data
replication or new MR/coloc. No SSH, source data transfer or credentials used.

## Gate status and next genuine scientific action

- Historical BBJ 2021/EAS5 numeric compatibility: **PASS**.
- Independently verified official effect-allele and genome-build semantics: **HOLD**.
- Freeze SHA-256 of original EAS5 archive and reproduce comparison on
  GitHub-hosted data, if an authorized non-server access path is available: **HOLD**.
- Exact EAS-AIS/GIGASTROKE source SNP beta/p/SE/N/case-fraction identity: **HOLD**.
- Per-locus cis-QTL harmonization, EAS LD reference and 25/36 candidate
  coloc.abf: **HOLD**.
- Causal gene claims, PP.H4 and independent replication: **NOT APPROVED**.

Do not claim this PR establishes new source identification by itself.
