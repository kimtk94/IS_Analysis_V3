# IS: EAS cis-eQTL variant-level feasibility, exact N, MAF and eCAVIAR signed-Z alternative

Updated 2026-11-16. **Internal methods/evidence and GitHub CI handoff, not a claim of ischemic stroke biological associations.**

## Verified source facts from completed real runs

| Check | Direct source evidence | Status |
|---|---|---|
| EUR+EAS GIGASTROKE IS × 1000G Phase3 22chr coordinate QC | [Original 22chr cross-run #43033560457](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43033560457) | PASS source coordinate/allele set; **not LD correlation matrix** |
| Japanese ImmuNexUT nominal eQTL actual records | [Original ZIP64 eQTL CI #43847910968](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43847910968) | PASS original 13-col CD16p Mono first 79 records; 75 P>0.05 |
| First complete ImmuNexUT gene cis block | [Original data #43848328705](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43848328705) | PASS WASH7P first 835/835 source-declared pairs, 775 P>0.05; **only one source gene** |
| Official cohort 416 WGS SDRF MD5 | [Source cohort #43433142595](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43433142595) | PASS, 416 source records, 79 healthy, 337 immune disease, **NOT cell-specific eQTL N** |
| Independent original eQTL+cohort reconciliation | [Cross-artifact #43849319075](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43849319075) | **SUCCESS tests 7/7, SOURCE_ROW_QC_PASS_COLOC_HOLD** |
| ImmuNexUT actual eQTL first gene ref/alt | WASH7P `chr1:54712-54716` `REF=TTTTC ALT=T` | actual source **indel** + multi-rsID; standard SNP join is insufficient |
| eQTL 13-field statistic availability | `Gene_id Gene_name CHR TSS_position Number_of_variants_cis Variant_ID REF ALT Variant_CHR Variant_position_start Variant_position_end nominal_P_value slope(ALT)` | **NO SE, MAF/EAF, per-celltype N** |
| GWAS and eQTL source builds | GIGASTROKE original GRCh37, ImmuNexUT original GRCh38 | liftover and reference-allele normalization **NOT_RUN** |

## New source audits submitted as GitHub Actions

1. **Exact celltype-N provenance, original expression source**: [E-GEAD-397 original publisher filelist](https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-397/E-GEAD-397.filelist.txt) has `E-GEAD-397.processed.zip` 286,026,411 bytes (publisher MD5 `8255d3aa3c0ccb6f703bd98d8593ffe7`), nested `count.tar` 1,310,883,840 uncompressed bytes and `clinical_diagnosis_age_sex_v2.txt` 8,256 bytes. [GitHub run #43850338467](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43850338467) registered with only bounded ZIP ranges and independent **small clinical file MD5**. **Do not conflate expression source sample count or 9,852 total celltype samples with genotype+expression+QC-complete N for CD16+ monocytes.**
2. **EAS GIGASTROKE outcome frequency distribution**: [GitHub run #43851027085](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43851027085), publisher-original full 260MB source EAS IS GRCh37 MD5, 6,725,216 full autosomal rows, EAF coverage/missingness, minor-frequency histograms, palindromic MAF≥0.4 sensitivity. **Outcome-GWAS EAF is not the Japanese immune-cell eQTL participant MAF** and cannot substitute by default in coloc.
3. **Numerical signed-Z input preflight for eCAVIAR**: [GitHub run #43851922254](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43851922254), `scripts/is/audit_ImmuNexUT_eCAVIAR_signed_z.py` and synthetic tests. Read real WASH7P first 835 SNP+indel association rows from the bounded original E-GEAD420 ZIP64 source; compute *derived* two-sided asymptotic signed Z as `sign(slope(ALT))*norm.isf(nominal_P_value/2)`. Preserve original GRCh38 REF/ALT unnormalized and report assumed normal/test calibration; no original SE or real GWAS-joined z/LD/CLPP. **WASH7P is a source-format sentinel, not an IS prioritized gene**.

All three runs were QUEUED at last direct check. Workflow existence and synthetic unit tests alone do not establish remote source PASS.

## Two distinct formal colocalization options

### A. coloc.abf / coloc.susie (primary after statistical eligibility)

[Official coloc input documentation](https://chr1swallace.github.io/coloc/articles/a02_data.html) prefers `beta+varbeta`; for a quantitative trait with unknown `sdY`, a defensible `N+MAF` may also be required. The alternative `pvalues+MAF+N` approximation is permitted in coloc.abf but **per-celltype N and matching MAF are not presently known**. The 13-col ImmuNexUT eQTL source does not publish `slope_se`. No assumptions about N=416, MAF derived from GWAS patients, or untested SE conversion are valid as evidence.

### B. eCAVIAR based on marginal signed Z and population-matched LD

[Original eCAVIAR paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5142122/) and [published statistical workflow](https://escholarship.org/uc/item/0rp7n0tg) describe per-locus **GWAS marginal Z scores, eQTL marginal Z scores, and corresponding LD correlation matrices**. For ImmuNexUT the original `slope(ALT)` can supply signed direction and the nominal two-sided P value a *derived* asymptotic Z magnitude. This is **an alternative input route** when beta SE, MAF, N are unavailable, not a validated eCAVIAR result or license to ignore uncertainty/calibration.

Before **any** CLPP: (1) prespecify GIGASTROKE IS locus, (2) lift GRCh37 GWAS variants to GRCh38 with chain, independently check matching reference bases and left-normalized indels, (3) harmonize outcome beta/SE and EAS QTL ALT-directed signed Z to a single reference allele, (4) intersect exactly the same variants and report variant loss, (5) compute a real local EAS PLINK genotype LD matrix from verified matching build/ancestry reference, (6) inspect matrix order/allele direction, (7) test 1 versus 2 causal variants, and (8) report uncertainty and CLPP cutoff provenance. Proxy population LD is a limitation; study genotype not publicly accessible.

Source technical completeness for **one** WASH7P gene (835/835) is insufficient to show tested-cis coverage of any preselected ischemic stroke locus. EAS GWAS case-control ascertainment, power, QTL cohort 337 autoimmune patients, tissue differences and cell-type effective N remain confounders. The primary novel claim must concern robust *ischemic stroke-specific locus-to-gene prioritization* and not generic EAS ancestry matching.

## GO/HOLD matrix

| Gate | Status | Required next proof |
|---|---|---|
| G1 genuine EUR/EAS GWAS & 1000G source structural matching | **PASS** | Actual local LD values and independent build validation |
| G2 genuine EAS all-nominal gene-level public eQTL source | **PILOT PASS** | IS locus genes, multiple cell types and full tested-variant coverage |
| G2 exact eQTL effective N, MAF, beta SE | **HOLD** | Independent provenance, never derive from WGS cohort totals alone |
| G2 eCAVIAR signed-Z input transform | **RUN SUBMITTED; NOT YET RESULT** | Real 835 variant Z numerical test and original P/beta consistency |
| G3 GRCh38↔GRCh37 allele/haplotype harmonization | **NOT_RUN** | Mapping test + unique reference bases + matched variant counts |
| G4 eCAVIAR/coloc.abf/coloc.susie | **BLOCKED** | Same set of true variants + matched LD, checked summary stats |
| G5 novelty, MAGMA, LDSC-SEG | **NOT_RUN** | Secondary once locus and ancestry model is calibrated |

## Added outcome-GWAS first-pass locus window screen (descriptive, not independent loci)

The EAS original-GWAS frequency QC script now also screens all original GIGASTROKE `GCST90104545` EAS IS variants with `p<5e-8` into 1 Mb **fixed genomic windows** (GRCh37), retaining the top beta/SE/P/EA/OA per window and the number of significant variants in that bin. It does **not** perform LD clumping, SuSiE, multi-ancestry prioritization or independence testing. Top signals are **NOT** defensibly called distinct causal loci, and genes cannot be assigned solely by distance. All source totals and published GWAS MD5 must pass.

- Latest descriptive screening CI: [#43852757654](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43852757654), **queued** as last verified.
- Only after this actual source CI is complete should preselected stroke GWAS intervals be annotated via independent GRCh37 gene coordinates, harmonized to matching ImmuNexUT GRCh38 loci, and assessed for per-locus tested-variant coverage.
- The derived GWAS outcome EAF/MAF is **not Japanese eQTL N/MAF**, and fixed 1 Mb bins are **not LD-defined**.
