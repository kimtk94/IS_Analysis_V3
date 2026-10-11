# IS priority four: expanded EUR503 full-cis signed LD QC

Date: 2026-10-11 KST. The source is public 1000 Genomes Phase 3 **EUR503 external reference genotypes**, NOT original GTEx individual-level QTL donor LD.

## Results

| QC metric | chr4 FGF5/C4orf22 | chr10 CALHM2/NEURL |
|---|---:|---:|
| Exact source SNP keys in VCF + exported PLINK raw | 2,834 | 2,828 |
| Samples | 503 | 503 |
| SNPs with any missing ALT genotype | 0 | 0 |
| Monomorphic reference SNPs | 9 | 13 |
| Variants informative for signed r against the GWAS lead | 2,825 | 2,815 |
| Reference BBJ lead | GRCh37 4:81182554:T:C | 10:105459834:T:C |
| Lead ALT allele frequency in EUR503 | 0.26739563 | 0.39363817 |
| 120-SNP reference signed correlation matrix lowest eigenvalue | -5.57e-15 | -1.63e-15 |
| 120-SNP reference signed correlation matrix maximum eigenvalue | 14.3555 | 16.4767 |

The minimum eigenvalue is negative only at numerical roundoff and satisfies the tolerance of -1e-8 for **these chosen 120-variant slices**. This does not validate inversion, full 2,800-SNP LD rank, sample adequacy, condition number or the original QTL cohort's LD. As the genotype panel has 503 individuals, full 2,800-variant LD is **necessarily rank deficient** and cannot be assumed invertible; a shrinkage/regularization sensitivity design is required before any full-locus conditional analysis.

## Allele direction

PLINK 2 `--export A` may report REF-counted or ALT-counted dosage. Every SNP column's counted allele was interpreted and converted to biallelic genomic **GRCh37 ALT** dosage. Unknown allele suffixes raise an error; source genotypes outside 0–2 are invalid. Missing, invariant and low-N SNPs are separately marked instead of assigning them an LD correlation of zero.

## Coverage and inference limits

Previous audit found exact allele coverage of 2,834/2,836 (chr4) and 2,828/2,830 (chr10). This run uses exactly those matched SNPs, not the four missing keys. No genotype was imputed. The inherited 646 coloc posteriors, new 120 ABF sensitivity runs and the 2,225-gene discovery universe remain unchanged.

The reference is useful for **exploratory** lead-tagging and possibly sensitivity-only EUR LD conditioning. It cannot be used to claim `coloc.susie` with in-study GTEx LD. Required research gates remain actual EBI QTL model/effective expression N, original cohort ancestry mix/LD, conditioning assumptions, convergence, credible set QC, priors and SNP-level effect orientation.

## Reproducibility

Input EUR503 regional VCF: `/srv/is-analysis/results/is/audits/is_priority4_ebi_bbj_eur_fullcis_20261011_v1/` (chr4 and chr10).
Extracted 503×(2,834/2,828) hardcall ALT dosage: `/srv/is-analysis/results/is/audits/is_priority4_eur503_fullcis_signed_ld_qc_20261011_v1/`.
Full per-SNP signed LD to BBJ lead, MAF, missingness, monomorphism and 120-variant matrix PSD summaries: `/srv/is-analysis/results/is/audits/is_priority4_eur503_fullcis_signed_ld_qc_20261011_v2/`.
Code: `scripts/is/audit_is_priority4_eur503_fullcis_signed_ld_qc.py`.

**No formal SuSiE-coloc was run.**
