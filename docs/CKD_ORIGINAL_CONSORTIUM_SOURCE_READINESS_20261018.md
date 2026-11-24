# CKD original consortium source readiness — 2026-10-18

**Executed:** [read-only GitHub Actions 39112953062](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39112953062), SUCCESS. This is **file inventory only**.

Known server entry points:
- `/srv/is-analysis/data/ckd/analysis_ready/UKBPPP_ST16_cis_independent.tsv.gz` (pQTL independent signals, **not necessarily all locus-level raw summary variants**).
- `/srv/is-analysis/results/ckd/stage2/stage2_pqtl_download_manifest.tsv` and `stage2_pqtl_schema.tsv` (next targets for resolving original cis-pQTL sources).
- `/srv/is-analysis/data/ckd/analysis_ready/EUR/eGFRcrea.tsv.gz` and `/srv/is-analysis/results/ckd/stage1/harmonized/EUR/eGFRcrea.tsv.gz` (transformed and/or analysis-ready summary; original download metadata must be traced).
- `/srv/is-analysis/data/ckd/rawdata/support/EUR/ckdgen_gorski2017_egfrcys_eur/CKDGen_1000Genomes_DiscoveryMeta_eGFRcys_overall.csv.gz` (support **eGFRcys**, NOT the primary **eGFRcrea** outcome).
- `/srv/is-analysis/data/ckd/analysis_ready/EAS/eGFRcrea.tsv.gz` (analysis-ready EAS, not native source independently verified).

Four MAF outliers already located:
- SDCCAG8 rs35646972, rs6663222
- INHBC rs11834807
- UMOD rs11860566

**Current scientific status:** Original consortium allele/MAF provenance is **PENDING**, and pQTL-MAF vs external EUR reference discordance has not been explained. File inventory does not demonstrate raw input coverage for a particular SNP. It is invalid to assert that UMOD's MAF discrepancy caused SuSiE nonconvergence.

**Next audit (no full download needed):** Read stage2 pQTL manifest and schema to locate actual original per-variant QTL files; read eGFRcrea GWAS manifest and inspect source columns; query only the four target variants in those files and compare build, effect/other allele, EAF, INFO, and outcome betas. Do not substitute eGFRcys for eGFRcrea.
