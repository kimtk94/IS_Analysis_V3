# Korean LD source qualification and aggregate request — 2026-11-28

Status: HOLD_STUDY_MATCHED_LD. Public metadata review only, no data application or email sent, no individual genotypes accessed. Fixed post-readout locus: six SNPs in the prior Phase 3 audit; not a new discovery panel.

## Source decisions

| Source | Public evidence | Decision |
|---|---|---|
| Lee et al. 2022 original Korean wrinkle GWAS | DOI 10.1016/j.jid.2022.03.027; PMID 35421404; GCST90094903/904, reported 11,079 each | Priority 1: ask original study team for study-specific signed genotype correlation or conditional regression, phenotype scale, model/covariates and overlap. No public study LD confirmed. |
| Korea4K | DOI 10.1093/gigascience/giae014; 4,157 WGS in article; EGA EGAD00001015348 describes 3,776 WGS datasets and 3,913 CRAM files; DAC EGAC00001003442 | Priority 2: independent Korean reference, ask for aggregate LD under approval. Article sample count not usable LD N. Public frequency download not a phased LD panel. No individual-level access on Actions. |
| KOVA2 | https://pmc.ncbi.nlm.nih.gov/articles/PMC9628380/ ; 5,305 healthy Koreans, mixed WES/WGS | Public frequency context only until joint genotypes or approved signed LD demonstrated; mixed capture and locus coverage must be checked. |
| KRGDB/Korean Imputation Service | https://pubmed.ncbi.nlm.nih.gov/32133509/ ; https://coda.nih.go.kr/usab/kis/intro.do | 1,722 WGS KRGDB resource and separate imputation service. Neither service is proof of downloadable signed LD. Ask operator about aggregate locus request. |
| Korea1K | https://pmc.ncbi.nlm.nih.gov/articles/PMC7385432/ ; 1,094 WGS | Approved access to individual genotype data, not public LD. Overlap with Korea4K must be assessed before calling independent. |
| KoGES | https://openkoges.com/pheweb/about | Public GWAS summary statistics, not genotype correlation; no matched wrinkle endpoint certified. |
| ASTAGR | DOI 10.1016/j.jid.2025.02.130; 2,762 women across four Asian countries | Independent outcome candidate only; verify site-specific phenotype, Korean subset, overlap, genotype coverage and access. No replication data received. |
| Korean facial-trait study | DOI 10.1016/j.ygeno.2023.110710; 85 facial traits, KoreanChip | Phenotype inventory lead only; do not assume wrinkle equivalence or independent cohort. |

Korea4K access policy: https://ega-archive.org/datasets/EGAD00001015348 (approved purpose and redistribution restrictions; raw data must remain in approved environment). Project page https://koreangenome.org/wiki/korea4k-genomes. KOVA2 site https://www.kobic.re.kr/kova/. Original GWAS https://pubmed.ncbi.nlm.nih.gov/35421404/. All links public metadata checked 2026-11-28.

## Frozen request

Use rs4584990 C/T (GRCh38 chr2:108872562), rs144918519 G/A, rs260643 A/G, rs115957039 T/C, rs6542772 C/T and rs3827760 A/G as in `skin_allele_resolved_phase3_ld_readback_20261128.json` (blob 70746dd0956e65416f33828184ea772ce0ac809b). Do not treat rsID alone as identity; request GRCh37/38 coordinate, REF/ALT, coded allele, allele frequency, per-marker usable N and INFO. Confirm exact alternate allele on the requested build. Freeze full nominal region GRCh38 chr2:108372562–109372562 from prior GCC2 audit, not only significant variants.

Preferred original-GWAS deliverable: signed **dosage correlation r**, allele-ordered matrix with diagonal 1, per-pair N, missingness policy, imputation QC, ancestry definition, phenotype scale/transform and covariate specification. Prefer covariate-compatible score covariance or original joint regression (GCC2 candidate conditional on EDAR; reciprocal conditional estimates), including beta, SE, P, N and model specification for both wrinkle phenotypes. Report all six markers, including missing/failed. Phased haplotype r² alone does not supply signed dosage covariance; independent Korean reference is sensitivity only, not study-matched.

Before any calculation: check exact variant order, alleles, matrix symmetry, positive semidefiniteness, diagonal, pair N, allele frequency parity, model/scale compatibility, collinearity, and overlap. Do not silently repair a non-PSD matrix or use arbitrary ridge regularization. Any incompatible inputs => HOLD. Conditional results are exploratory post-readout, not confirmatory. Regional coloc requires >=100 shared and >=80% both denominators, matched QTL regional nominal data and ancestry/conditional review as previously frozen.

## Draft contact text (not sent)

Dear research team,

We are assessing the regulatory context of the GCC2–EDAR region in the published Korean facial-wrinkle GWAS (GCST90094903/904). We would like to ask whether an authorized aggregate-only analysis would be feasible. For six fixed variants and the previously defined regional window, could you provide allele-aligned signed dosage correlations with per-marker/pair sample counts and quality information, or fit the corresponding joint association models in your approved environment? We would also appreciate the exact outcome scale, covariate/model definitions and sample overlap between the two sites. We do not request individual participant records or raw genotypes. Any results would be used only within the agreed data-sharing and publication conditions.

## Stop rule

No full cohort genotype download, data application, request transmission, or further conditional/coloc claims without appropriate access. Keep CHB/JPT/CEU LD descriptive. If study-matched LD and compatible effects cannot be obtained, retain locus ambiguity in manuscript, rather than claiming independence of GCC2 from EDAR.
