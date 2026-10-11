# IS priority-four — BBJ vs GIGASTROKE EAS AIS regional lead structure

**2026-10-11 KST.** Exploratory two-locus summary GWAS comparison; **not in-study conditional GWAS, fine-mapping or causal inference**. Crucial selection boundary: the SNP universe here comprises only **exact allele-matched, non-palindromic SNPs included in EBI QTD000161 full-cis molecular tests and both BBJ and GIGASTROKE AIS GWAS**. A "lead" is therefore an *intersection-restricted* most significant SNP, not the validated overall GWAS locus lead.

| Locus | Exact shared variants, union over two gene labels | BBJ most significant shared variant | p BBJ | GIGASTROKE EAS AIS most significant shared variant | p GIGASTROKE | EAS504 lead–lead signed r | EAS504 lead–lead r² |
|---|---:|---|---:|---|---:|---:|---:|
| chr4 FGF5/C4orf22 | 2,697 | `4:81182554:T:C` | 1.212e-8 | **same** | 2.888e-8 | 1 | 1 |
| chr10 CALHM2/NEURL | 2,626 | `10:105459834:T:C` | 8.762e-10 | `10:105599770:T:G` | 5.584e-11 | **+0.08586** | **0.00737** |

Chr10 BBJ lead's cross-association in GIGASTROKE: **p=7.792e-7**. GIGASTROKE's lead cross-association in BBJ: **p=9.246e-8**. Neither cross-association reaches the two-sided genome-wide significance threshold of 5e-8 as an identical SNP. Conversely, both marginal p-values are smaller than 1e-6, making it inaccurate to claim no supporting regional association.

The two chr10 lead SNPs are only weakly tagged in the public 1000G ancestry-aligned reference: **EAS504 r=+0.08586 (r²=.00737)**, **JPT104 r=+0.15410 (r²=.02375)**, **non-JPT EAS400 r=+0.06789 (r²=.00461)**. These are observed **reference-genotype correlations**, not correlation estimates from BBJ or GIGASTROKE cases and controls. Weak pairwise LD is compatible with multiple marginal signals but does **not prove two independent causal effects**, particularly without EAS study matched conditional statistics and fine-mapping.

FGF5 and C4orf22 share a single chr4 genomic GWAS locus and its lead SNP. CALHM2 and NEURL share a single chr10 locus even though the two GWAS report different strongest variants within the restricted SNP universe. Avoid fourfold pseudoreplication.

Additional descriptive direction agreement between GWAS effect-aligned SNPs was 84.0% in chr4 and 89.0% in chr10, but nearby SNPs are LD-correlated: these percentages are not replication probabilities or p-values.

## Inputs, scope, reproducibility

Frozen source: `/srv/is-analysis/results/is/audits/is_priority4_fullcis_bbj_giga_snp_exact_20261011_v1/IS_PRIORITY4_FULLCIS_EBI_BBJ_GIGASTROKE_SNP_EVIDENCE.tsv` (20,689 QTL gene–variant rows).

EAS reference genotypes: `/srv/is-analysis/results/is/audits/is_priority4_eur503_eas504_conditional_20261011_v1/`; 504 subjects, GRCh37 genomic ALT orientation.

New script: `scripts/is/audit_is_priority4_locus_lead_structure.py`; research-server output `/srv/is-analysis/results/is/audits/is_priority4_regional_leads_20261011_v1/`, with lead SNP summary, cross-GWAS lookup and provenance manifest.

This audit **does not** establish full-gene cis or full-GWAS regional conditional signals, sample independence, BBJ–GIGASTROKE overlap, or QTL donor cohort LD. Existing 646 inherited coloc posteriors, newer 120 EBI–BBJ ABF sensitivity results, and the 2,225 positional IS gene discovery universe were left unchanged.

**Next analytical gate:** extract **studywide raw GIGASTROKE EAS AIS GWAS summary variants across both loci** and identify the truly strongest *unrestricted* signal, comparing that lead to BBJ's unrestricted locus lead. Only then test whether ancestry-appropriate cohort LD and formal conditional analysis support one or multiple signals. This will prevent restriction-to-QTL bias.
