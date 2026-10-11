# Priority-four EBI–BBJ EUR503 full-cis LD reference coverage (2026-10-11)

Reference: 1000 Genomes phase3 GRCh37 regional VCF, EUR503 individuals; UCSC indexed distribution; re-extracted chr4:80,087,872–82,356,697 and chr10:104,113,502–106,353,104. The two regional VCFs contain 62,011 and 59,329 variant records, respectively.

| Locus | ABF unique matched variants | Earlier EUR reference exact SNPs | New EUR reference exact SNPs | Missing |
|---|---:|---:|---:|---:|
| BBJ_IS_L001 (FGF5/C4orf22) | 2,836 | 1,773 (62.52%) | **2,834 (99.93%)** | 2 |
| BBJ_IS_L002 (CALHM2/NEURL) | 2,830 | 1,773 (62.65%) | **2,828 (99.93%)** | 2 |

Missing exact chr:pos:REF:ALT keys: chr4:80478866:A:C, chr4:82079041:C:T, chr10:104510493:T:C, chr10:105222770:G:A. Presence in the exact indexed panel was checked, not imputed or extrapolated.

**Important:** This is a completeness audit for the EUR 1000G *reference* genotypes. It is not GTEx donor-cohort QTL LD and does not validate fine-mapping, a causal SNP, or coloc.susie. It does not verify genotype missingness, monomorphism, PSD conditioning, phased allele direction, or study cohort representativeness. Those checks must precede exploratory LD-based conditioning. All 2,225 parent IS positional genes remain in scope. Original 646 archived ABF posteriors and recent 120 pilot scenarios remain unchanged.

Code: `scripts/is/audit_is_priority4_ebi_bbj_eur_ld_coverage.py`.
Audit files: `/srv/is-analysis/results/is/audits/is_priority4_ebi_bbj_eur_ld_coverage_20261011_v2/` (per-variant exact match table and JSON/TSV summary).
Regional EUR VCFs: `/srv/is-analysis/results/is/audits/is_priority4_ebi_bbj_eur_fullcis_20261011_v1/`.
