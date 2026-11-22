# IS P3A: 36 genuine QTL assays, publisher-original GRCh37 GWAS and reciprocal SNP check

**2026-11-22 | MASTER DEGREE / ischemic stroke only.** 20 distinct genes, 36 actual gene–tissue QTL assays, a single chr10 locus (G0015). This is source technical QC, not a formal Coloc result or 36 independent GWAS loci.

## Authoritative source and CI evidence

- [Actual EAS+EUR original publication and three original QTL files GitHub Actions #45031834139](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45031834139): **SUCCESS**. New fail-closed 20/20 tests; repository-wide 320 executed, 319 pass, one skipped because Rscript unavailable.
- Molecular QTL: eQTL Catalogue Release 7, QTD000131 GTEx aorta, QTD000171 GTEx cortex, QTD000434 ROSMAP DLPFC, **GRCh38 ALT-effect**. Original indexed association source query over shared 20-gene Ensembl105 reference union at GRCh38 chr10:101,837,472–105,353,833, 3,516,362 bp.
- Reference: Official UCSC hg19→hg38 plus hg38→hg19 chains (published MD5 35887f73fe5e2231656504d1f6430900 and ff3031d93792f4cbb86af44055efd903); two actual UCSC chr10 FASTA files with previously checked original MD5 and SHA-256.
- Original GIGASTROKE published GRCh37 files were downloaded to ephemeral GitHub runner and **entire publisher MD5, SHA-256, file size, original row count** verified: EAS AIS GCST90104545 260,689,208 compressed bytes and 6,725,216 rows; EUR AIS GCST90104540 283,482,604 compressed bytes and 7,482,032 rows. No original SNP effect results uploaded.

## Actual source SNP alignment across all 36 assays

For each canonical nonpalindromic single SNP in the original QTL source: require unique GRCh38→GRCh37→GRCh38 return position, plus strand and chromosome identity, exact GRCh38 and GRCh37 FASTA REF; then compare GIGASTROKE GWAS original effect_allele and other_allele with QTL ALT and REF. **Flip GWAS signed beta when the original GWAS effect allele is QTL REF** and orient original EAF to ALT. Do not equate QTL minor allele frequency with ALT allele frequency; source ac/an can provide separate observed QTL ALT frequency, but not effective cohort N.

| Source QTL | Assay groups | Valid nonpal QTL SNP×gene rows | EAS signed-allele matched | EAS absent GWAS pos | EUR signed-allele matched | EUR absent GWAS pos | EUR discordant-allele rows |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| QTD000131 / GTEx aorta | 10 | 44,693 | 25,593 | 19,100 | 35,895 | 8,798 | 0 |
| QTD000171 / GTEx cortex | 11 | 45,570 | 27,158 | 18,412 | 38,521 | 7,049 | 0 |
| QTD000434 / ROSMAP DLPFC | 15 | 58,081 | 36,846 | 21,235 | 53,866 | 4,207 | 8 |
| **TOTAL** | **36** | **148,344** | **89,597 (60.4%)** | **58,747** | **128,282 (86.5%)** | **20,054** | **8** |

- All 148,344 valid nonpal QTL SNP×gene×tissue records mapped reciprocally and passed both reference FASTAs: **0 reference or chain blockers** across the 3 QTL sources. These are NOT 148,344 distinct genomic SNPs (same SNP may occur under multiple genes/tissues). The source dataset-level unique GRCh38 pos+REF counts are aorta 7,312, cortex 7,044, ROSMAP 6,736 (overlapping between tissue datasets).
- Original target GWAS genomic positions found: EAS **4,359**, EUR **6,444**. No original target GWAS duplicate position or invalid beta/SE/EAF/allele record.
- EAS matching effect allele identical to QTL ALT: **48,711** gene-assay records; reference allele original effect requiring beta sign flip **40,886**. EUR matching ALT **72,615**, reference effect requiring sign flip **55,667**.
- All 36 groups contained at least one matched nonpal QTL×GWAS variant in **both** EAS and EUR. EAS matched 60.4% of candidate QTL nonpal records; EUR 86.5%. Coverage is descriptive and not an ancestry-matched LD or replication statistic.
- EUR mismatches: **8 ROSMAP DLPFC SNP×gene assay records** across eight gene rows; could reflect repeated views of one genomic SNP rather than eight distinct loci. Hold rather than force allele direction.
- Difference between observed QTL ALT frequency (ac/an) and alt-oriented GWAS outcome EAF greater than 0.2 occurred in **14,054 EAS matched gene-assay records and 0 EUR**. Different QTL and stroke GWAS ancestral composition makes these descriptive only, not independent genotype error proof.
- Prior 3.52Mb expanded original source also includes 44,454 valid source INDEL/multibase association records: not normalized/converted in this SNP-only analysis. Palindromic A/T and C/G, including original historical lead, similarly remain a separate HOLD class.

## Literature corroboration of phenotype sample counts (NOT per-variant N approval)

Later studies referencing Mishra et al. GIGASTROKE 2022 report **EAS AIS 19,032 cases / 237,242 controls (N=256,274)** and **EUR AIS 62,100 cases / 1,234,808 controls (N=1,296,908)**. External corroboration: https://pmc.ncbi.nlm.nih.gov/articles/PMC13015613/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC10081428/ . These numbers must still be checked against exact original publisher study release phenotype/cohort metadata and overlap before use as approved model N or case fraction. GIGASTROKE EAS contains BBJ study contributions, so do not claim BBJ-independent replication.

## Remaining scientific gates: original issue #307 stays OPEN

1. Exactly matched original source-study phenotype positions, actual gene-specific nominal tested cis variants, and QTL study-specific exposure effect donor N cannot be inferred from GTEx/ROSMAP metadata cohort N 387 / 205 / 560. The 3.52Mb reference union is a positional reference survey, not proof every original tested cis variant is present.
2. GWAS exact cohort/phenotype Ncase/Ncontrol, per-variant representation and BBJ overlap require primary release approval; no arbitrary effective N or imputed standard errors.
3. Preserve eight discordant EUR ROSMAP records, genuine source INDELs and palindromic strand alleles for later independent normalization and sensitivity. Never silently approve them.
4. **Formal coloc.abf eligible/executed 0/36; formal coloc.susie eligible/executed 0/36**; no PP.H4 calculated, no causal gene claims and no independent replication label.

Implementation: scripts/is/audit_p3a_36_real_gwas_qtl_reciprocal_snp.py; tests/test_is_p3a_36_real_gwas_qtl_reciprocal_snp.py; .github/workflows/is-p3a-36-reciprocal-snp-gwas.yml. GitHub artifacts contain only original publisher metadata, source hashes, group-level QC counts and explicit HOLD statuses.
