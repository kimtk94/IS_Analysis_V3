# IS P3A GWAS source matching — provenance first (2026-10-18)

## Current state (do not upgrade without reproduced matches)

Harness acceptance [#193](https://github.com/kimtk94/IS_Analysis_V3/issues/193) completed PASS, but `P3A_APPROVED_GWAS_SOURCES.tsv` is **MISSING_APPROVED_REGISTRY**. Consequently, `accepted_sources=[]`, neither BBJ nor EAS-AIS sample size is adopted, and 25 source-stat-positive P3A candidates remain blocked. This is correct.

## External sources — candidates, not approved matches

1. **BBJ ischemic stroke**, `bbj-a-129`, OpenGWAS: https://opengwas.io/datasets/bbj-a-129 . Metadata reports East Asian, GRCh37, 17,671 cases / 192,383 controls, 8,885,705 SNPs. These sample sizes do **not** become study parameters until archived ALT beta/p match to downloaded accession-identical official GWAS.
2. **GIGASTROKE East Asian AIS**, candidate `GCST90104545` in GWAS Catalog. The PGS Catalog enumerates accession `GCST90104545` among stroke GWAS and indicates 256,274 EAS individuals. This aggregate is *not* a verified exact case/control split for archived EAS-AIS. Primary Mishra et al. (2022) publication: https://www.nature.com/articles/s41586-022-05165-3 . Data statement: https://pmc.ncbi.nlm.nih.gov/articles/PMC9524349/ .
3. Official GWAS Catalog source file layout and release metadata: https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics/ ; accession-specific FTP download. OpenGWAS VCF downloads can require authenticated token and time-limited URL; never store credentials in public GitHub Issues.

## Independent original-file matching acceptance protocol

- **Primary file**: download source-specific original release using authorised access, record accession, phenotype label, ancestry, assembly/build, release date, source URL, local absolute path under approved data root, size, SHA-256, and selected columns.
- **Variant coordinates**: independently normalize original row to `chr:pos:REF:ALT` GRCh37 using an explicitly approved reference; exclude multiallelic, invalid, ambiguous strand or incompatible build records. Don't infer GRCh37 from a filename.
- **Signed effects**: verify the effect-allele orientation. `beta_ALT` requires a documented flip when the original effect allele equals REF; record flip and REF/ALT checks.
- **Numerical reproduction**: choose at least 20 *distinct, nonpalindromic* matching variants distributed across at least two loci; compare archived ALT beta and p with raw transformed GWAS values. Preserve matching/discordant/ambiguous counts and absolute/log-p error diagnostics. Real thresholds must be justified by the original release's rounding precision.
- **Statistical metadata**: only after independent source identity is accepted, adopt Ncases/Ncontrols and case fraction for the *same phenotype and sample set*. Case/control fractions from unrelated PGS evaluation samples or general East-Asian stroke analyses are invalid substitutes.
- **P3A `coloc.abf`**: requires trait-consistent QTL/GWAS summary stats including beta/SE or appropriate p/MAF method, sample sizes/case fraction, independent allele checks, predeclared full cis-region coverage, prior sensitivity and single-causal-variant caveats. No PP.H4 produced by this protocol.

### Implementation status

The repository has `scripts/is/harness_archived_evidence.py` with a fail-closed optional approval manifest and synthetic tests. The actual independent source-to-archive comparator and source-file access are **not verified as installed**: the attempted automatic integration was blocked by tool safety checks and was not claimed as deployed. Do not mark either study as matched or run coloc from source names alone.

Next permitted step: ingest exact original raw GWAS files through an authorised route, then perform an isolated read-only SNP comparison in a reviewed workflow before attaching its verification artifact to Harness.
