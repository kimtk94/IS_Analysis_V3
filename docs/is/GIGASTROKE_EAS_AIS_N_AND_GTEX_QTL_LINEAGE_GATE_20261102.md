# IS P3A — GIGASTROKE EAS-AIS case/control N and GTEx cis-eQTL lineage (2026-11-02)

## Scientific scope
Prior actual original-data GWAS matching remains valid: official original
GCST90104545 file integrity, genome-wide 6,725,216 rows, numeric beta/SE/p/EAF
and **6,273 / 6,273 matched unique SNPs** across the four-gene EBI
reprocessed full-cis QTL universe (5,326 nonpalindromic, 947 palindromic
numeric matches). [Actions #41439023099](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41439023099).

**This report is a separate gating analysis**, not retrospective promotion of
any previous MR or colocalization.

## 1. EAS-AIS reference sample size

GIGASTROKE East Asian **Any Ischemic Stroke** is official GWAS Catalog
GCST90104545, not GCST90104544 (Any Stroke). The published EAS subtype table
provides 19,032 AIS cases and 237,242 controls, total 256,274. Calculations
*conditional on treating the literature numbers as one uniform case-control
set* give:
- Literature-level case fraction, s = 19,032 / 256,274 ≈ 0.0742642640.
- Illustrative effective N (4 * Ncase * Ncontrol / Ntotal) ≈ 70,474.4101.

**Neither s nor Neff is approved for a given original GWAS variant.**
The official original file has beta, standard_error, p_value, effect and
other allele, effect_allele_frequency, but no per-variant sample size in
its verified 11-column header. Source companion YAML may or may not provide
case_count/control_count/sample_size. A multi-cohort meta-analysis can have
variant-specific missingness, heterogeneity and effective sample sizes.
The GWAS Catalog [YAML field definitions](https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics/)
say case_count and control_count can exist for case-control studies but do
not guarantee these fields exist for every historical pre-GWAS-SSF study.

The official YAML's accession, original data filename, MD5 and genome build
were already verified against the full source file. New independent CI now
extracts **actual** case_count/control_count/sample_size fields and will
report an explicit HOLD if they are absent, inconsistent, or study-level only.
If counts happen to agree with published values, the CI still refuses to
approve variant-level case fraction or uniform effective N.

## 2. GTEx v8 native vs reprocessed eQTL

These are scientifically **distinct analysis products**:
- [Original GTEx Portal V8](https://gtexportal.org/home/methods): GTEx v8
  FastQTL cis-eQTL with nominal gene-variant tests, 1 Mb cis windows,
  adaptive permutation and original sample QC/normalization.
- [EMBL-EBI eQTL Catalogue release 6](https://www.ebi.ac.uk/eqtl/Release_notes/):
  the label **GTEx** denotes its uniformly reprocessed GTEx source data
  (QTS000015 and QTD000161 here), including its own quantification, imputation,
  covariates and statistical workflow; \`.all.tsv.gz\` files contain the
  reprocessed complete nominal cis-eQTL tests.
- [eQTL Catalogue Studies](https://www.ebi.ac.uk/eqtl/Studies/): the separately
  listed **imported GTEx_v8** corresponds to the original external analytical
  product. The eQTL Catalogue explicitly warns that imported and reprocessed
  workflows can yield systematic differences in eQTL results.
- [2026 r8 pre-release notes](https://www.ebi.ac.uk/eqtl/Release_notes/)
  mention a GTEx v10 update for **QTS000015**, with dataset IDs retained.
  Consequently **QTD000161 by itself is NOT a sufficient immutable release
  identifier**. Must pin release version, compressed file content hash and
  donor/sample/quantification metadata.

Source-level access to the current QTD000161 URL, if verified, does not
automatically prove its bytes are the 2023 r6 input archived in the user's
prior 2026-10-11 four-gene analysis. Nor does r6 reprocessed QTL equivalence
prove it equals the original GTEx v8 FastQTL inputs used by historical
colocalization.

## 3. CI evidence and decision table

New CI pipeline:
- \`scripts/is/audit_is_eas_ais_n_eqtl_lineage.py\`
- \`tests/test_is_eas_ais_n_eqtl_lineage.py\`
- \`.github/workflows/is-eas-ais-n-eqtl-lineage.yml\`
- [GitHub Actions #41440738101](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41440738101)

This reads only bounded public metadata:
1. Exact accession's official original GWAS companion YAML, recording
   sample_size, case_count, control_count **only when those fields exist**,
   with actual metadata SHA256 and original source MD5 verification.
2. Official eQTL Catalogue study index + release notes (including original
   imported GTEx_v8 versus QTS000015 reprocessed), with source-page hashes.
3. Bounded original QTD000161 gzip prefix check, **not** a full QTL re-analysis.

Until the actual CI log and artifact are reviewed, results for new metadata
gates must not be claimed PASS. Whatever the metadata screen result,
**scientific_status = HOLD** for exact GWAS case fraction, original GTEx v8
FastQTL equivalence, ancestry-matched LD, locus colocalization and causal gene
promotion.

## 4. Work remaining to approve proper coloc

1. Pin the exact four-gene molecular QTL statistical source and version
   (if original FastQTL: retrieve original GTEx V8 per-tissue cis tested file;
   otherwise present EBI r6 reprocessed QTL as an explicitly distinct
   hypothesis/test).
2. Verify complete cis tested-SNP universe and effect allele orientation
   for *that specific* QTL model.
3. Resolve GWAS study-level versus SNP-specific case N and controls / case
   fraction; if one fixed N must be assumed, explicitly justify and run N/s
   sensitivity rather than silently adopting a published cohort total.
4. Confirm properly matched ancestry-aware LD and credible-set assumptions.
5. Only then compute preregistered MR/coloc sensitivity; never treat
   overlapping BBJ/GIGASTROKE EAS GWAS as independent replication.

**Scientific conclusion: GWAS source crosswalk PASS (limited four-gene
6,273-SNP universe); original molecular-QTL model equivalence and
case/control N for coloc HOLD.**
