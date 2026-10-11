# IS priority-four — BBJ vs GIGASTROKE EAS AIS exact-allele effect direction

**2026-10-11 KST.** Descriptive audit of published ancestry-aligned summary statistics; **not** independent causal replication, meta-analysis, formal heterogeneity test, or gene-level inference.

| Gene | Shared exact SNPs, non-palindromic | Same sign of genomic ALT beta | Direction fraction | QTD nominal p<1e-4 | BBJ p<5e-8 | GIGASTROKE EAS AIS p<5e-8 | Both stroke datasets p<5e-8, identical SNP |
|---|---:|---:|---:|---:|---:|---:|---:|
| FGF5 | 2,634 | 2,216 | **84.13%** | 5 | 1 | 1 | 1 |
| C4orf22 | 2,633 | 2,215 | **84.12%** | 0 | 1 | 1 | 1 |
| CALHM2 | 2,566 | 2,292 | **89.32%** | 3 | 11 | 53 | **0** |
| NEURL | 2,535 | 2,281 | **89.98%** | 0 | 11 | 53 | **0** |

Sources: previous frozen `IS_PRIORITY4_EBI_BBJ_GIGASTROKE_SNP_EVIDENCE.tsv`, where each comparison already underwent GRCh38→GRCh37 lift, FASTA REF alignment and BBJ/GIGASTROKE ALT beta harmonization. The full molecular QTD cis universe is 20,689 gene–SNP rows; the above counts require both stroke datasets and non-palindromic alleles.

**Scientific interpretation:** A SNP-direction fraction is strongly influenced by LD among SNPs and by two genes sharing the same stroke region. There are **two genomic loci**, not four independent association studies. A same-direction SNP is not an independent replication. The present work does not determine whether BBJ and GIGASTROKE summary samples overlap. Consequently do not attach a binomial p-value to SNP sign concordance and do not calculate a naive fixed-effect pooled p-value.

At chr10, the BBJ-significant and GIGASTROKE-significant lists do **not** contain any identical SNP satisfying both thresholds within this non-palindromic common SNP set, although locus-level direction agreement is high. This does not exclude a genuine regional trait association; it prevents the misleading claim that identical genome-wide-significant SNPs have replicated here.

At chr4, one common significant SNP is counted for FGF5 and C4orf22 because these two genes share a locus; it is not two separate replicated GWAS observations.

**Further formal comparison gates:** establish BBJ/GIGASTROKE study sample overlap and trait definition correspondence, use EAS LD-aware cross-GWAS regional/conditional analyses with full variant coverage and multiple-testing control, and separately authenticate QTL donor LD/model and N. GCTA-COJO remains blocked by inadequate reference genotype N (public EAS504), as documented in `IS_PRIORITY4_BBJ_EAS_COJO_READINESS_20261011.md`.

Files: `scripts/is/audit_is_priority4_bbj_giga_direction.py`; server `/srv/is-analysis/results/is/audits/is_priority4_bbj_giga_direction_20261011_v1/` includes 4-gene TSV summary and JSON interpretation manifest. Neither the 646 inherited coloc outputs, 120 new EBI–BBJ ABF scenarios, nor all 2,225 IS candidate genes were modified.
