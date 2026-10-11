# IS priority-four EUR503 LD blocks and shrinkage sensitivity

**2026-10-11 KST — reference-only numerical QC, NOT causal fine-mapping.**

We partitioned the existing 503-individual 1000 Genomes EUR full-cis biallelic SNP dosages (GRCh37 REF/ALT verified) into fixed non-overlapping **250kb physical windows**. Fixed intervals are technical tiling choices, **not independent association signals or LD recombination blocks**. Every PLINK `--export A` genotype was converted to genomic ALT dosage, missing values were rejected and monomorphic SNPs excluded from matrix construction.

| Metric | BBJ_IS_L001, chr4 | BBJ_IS_L002, chr10 |
|---|---:|---:|
| Exact EUR503 SNPs in exported source | 2,834 | 2,828 |
| Monomorphic SNPs | 9 | 13 |
| Polymorphic SNPs | 2,825 | 2,815 |
| 250kb tiles, at least two polymorphic SNPs | 10 | 10 |
| Largest tile polymorphic SNPs | 541 | 427 |
| Empirically rank-deficient tiles | 10 | 10 |
| Max estimated condition number under alpha=0.05 | 2,710 | 3,093 |

A diagnostic **correlation shrinkage** family `R_alpha=(1-alpha)*R+alpha*I` was considered for `alpha=0,0.01,0.05,0.10,0.20`. The spectrum was calculated efficiently from the **503×503 standardized genotype Gram matrix**, whose nonzero eigenvalues match those of the SNP correlation matrix; zero eigenvalues were inferred from rank deficiency. This provides eigenvalue/condition-number diagnostics **without falsely treating a large, singular full-locus matrix as invertible**.

**Scientific caveat:** the method computes eigenvalue bounds, not a fine-mapping likelihood, posterior, SuSiE convergence test or conditional association regression. The numerical conditioning from `alpha=0.05` is a **sensitivity scenario, not an endorsed choice**. Shrinkage does not make EUR 1000 Genomes genotypes equivalent to the eQTL Catalogue GTEx donor cohort LD. Ancestry differences, QTL effective expression N and variant-level alleles still matter. In particular, the original ischemic-stroke BBJ GWAS is primarily EAS, while this panel is EUR.

Data and code:
- `scripts/is/audit_is_priority4_eur503_ld_blocks_shrinkage.py`
- Server: `/srv/is-analysis/results/is/audits/is_priority4_eur503_ld_blocks_shrinkage_20261011_v1/`, containing block/alpha grid, locus summary, source SHA256 manifest
- Raw dosages: `/srv/is-analysis/results/is/audits/is_priority4_eur503_fullcis_signed_ld_qc_20261011_v1/`

No existing GTEx V8 or EBI-BBJ ABF PP.H4 results, GTEx donor LD, or the original **2,225-gene** discovery universe were changed.
