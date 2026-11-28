# ASTAGR phenotype qualification and external context decision — 2026-11-28

Sources: Moon et al., JID 2025 doi:10.1016/j.jid.2025.02.130, public article https://www.jidonline.org/article/S0022-202X%2825%2900290-8/pdf (Methods: Participants, Skin wrinkle measurement, Phenotype data processing, GWAS). Catalog metadata actual readback run 46132165733 job 140939499803, harmonised metadata SHA256 127e28d6fc93be954e8aaa129882da4662e81094a62f861b993043b2225314fe. Publisher text accessed; screenshot access returned HTTP 403.

## Qualified endpoint

ASTAGR GCST90502802: 2,424 women after QC, 20–59 years, four countries (Korea, Vietnam, Thailand, Indonesia); pooled results. Noncontact three-dimensional Primos CR measures arithmetic mean surface roughness Ra (micrometers) in a 12-mm diameter area around either eye. The exact distinction between inferior eye and lateral canthus is not specified in the Methods. GWAS uses rank-based inverse-normal transformation of age-adjusted residuals, PLINK 2 linear regression with 20 genetic PCs and batch covariates. Metadata also lists age. MAF >=0.05, call rate >=95%, HWE P >=1e-6, related individuals removed. Sequencing coverage 30x in three countries and 10x in Indonesia. The article states ancestry-specific power was limited. No exact same-site endpoint or overlap verification with Lee et al. 2022 GCST90094903/904.

Decision fixed before ASTAGR SNP effects are read: external **periocular roughness context**, not site-specific or Korean-only independent replication. Beta scale not matched, do not pool raw coefficients, do not infer a significant cross-site effect difference from differential P-values. Report all six fixed SNPs, exact coordinates/alleles, source EAF, nominal P, descriptive Bonferroni-six, missing and direction; no proxy substitution. If a variant is absent, MAF filtering is possible but not proven unless source evidence supports it. No broad 26-phenotype discovery in this analysis. External null does not rule out a site-specific effect.

## Additional candidate sources

- Lee KC et al. 2022, PMID 35313536, https://pmc.ncbi.nlm.nih.gov/articles/PMC8933694/ : 749 Korean women aged 30–50, canthus fine wrinkle measurements, Kolmar Korea recruitment 2019–2020. Possible lateral-eye context; exact marker coverage, outcome scale, sample overlap and summary statistics access unverified. No public full GWAS confirmed.
- Lee DS et al. 2022, doi:10.3390/app122211422: 261 Korean women, P&K Skin Research Center July 2021–January 2022, lateral eye and glabellar roughness/depth, age-adjusted regression. Possible site-specific context; small N, no full fixed marker response or cohort non-overlap confirmed.
- Prior Korean facial morphology GWAS PMID 37734486 / doi:10.1016/j.ygeno.2023.110710: facial morphology, not automatically wrinkles. Deprioritize without matching phenotype.
- EDAR prior functional work doi:10.3390/biom14030279 (2024), from LG H&H; EDAR wrinkle biology already studied. Novelty depends on distinguishing GCC2 regulation from existing EDAR association, not discovering EDAR as a wrinkle gene.

Hold: original Korean phenotype scale, joint outcome covariance, Korean study-specific dosage LD and conditional regression, complete regional molecular QTL and genuinely independent same-site replication. No contact sent. No main merge.
