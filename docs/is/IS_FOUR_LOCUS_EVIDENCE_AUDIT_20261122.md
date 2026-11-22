# IS four legacy loci: original-file readback and scientific-claim boundary (2026-11-22)

**Research scope:** MASTER DEGREE BBJ ischemic stroke genetics, **not** Sales prospecting or the CKM→Stroke P0–P6 pipeline.  
**Verdict:** original *analysis input table integrity* locally checked; `coloc` historical, formal independent replication / MR causality **NOT APPROVED**.

## Real private-source file inspection

On 2026-11-22 the archived four user-controlled Google Drive input `stats.tsv` files were **actually read individually** (SHA256 computed on full original file bytes). The filenames are repeated inside distinct archived `P##__BBJ_IS_...` folders. For privacy, Google Drive file IDs, signed download URLs, and individual variant rows are not copied into this public repository.

| Gene / archived folder prefix | Source rows = distinct SNPs | Source-byte SHA256 short | QTL N written in *every* input row | A/T or C/G palindromes | |MAF_GWAS−MAF_eQTL| > 0.2 |
|---|---:|---|---:|---:|---:|
| FGF5 / P01, BBJ_IS_L001 | 1,639 | `65135b9e984c` | 175 | 293 | 438 |
| ALDH2 / P11, BBJ_IS_L003 | 1,796 | `e4c5dae71eea` | 584 | 233 | 454 |
| SH3PXD2A / P14, BBJ_IS_L002 | 1,557 | `bf02be38e1bf` | 584 | 212 | 106 |
| COL4A2 / P09, BBJ_IS_L004 | 2,977 | `c074094bfe12` | 170 | 380 | 329 |
| **Four input tables** | **7,969** | Full hashes in the JSON evidence | Fixed per file | **1,118** | **1,327** |

All four have exactly `snp, gwas_beta, gwas_se, gwas_maf, eqtl_beta, eqtl_se, eqtl_maf, qtl_n` columns; all numeric GWAS/eQTL beta and SE values are finite, SE > 0, reported MAF in (0,0.5], no duplicated full variant strings, and one chromosome per file. **A constant `qtl_n` column is not source-proof of donor-level effective association N or of the sample-size model used in earlier coloc.** The 1,327 high cross-ancestry MAF gaps are flags for further allele / cohort reconciliation, **not proven allele flips**.

Each archived source folder also has `gtx_eur_ld.f64.bin` and `bbj_eas_ld.f64.bin`. The Google Drive *metadata byte sizes* for each are exactly `8×(number of stats.tsv rows)^2`, consistent with a dense float64 square matrix. This is **shape-size evidence only**: raw LD values, byte hash, SNP ordering, matrix PSD and original eQTL-donor/BBJ-in-study LD match were **not revalidated**. The `GTEx_V8` folder label and effect/allele build cannot replace original production provenance.

Machine-readable source hashes, observed aggregate counts and fail-closed scientific flags: [frozen audit JSON](evidence/IS_FOUR_LOCUS_PRIVATE_INPUT_QC_20261122.json). The optional `--source-root` argument of `scripts/is/audit_four_locus_private_inputs.py` can re-read the original four files from an authorized locally staged directory without publishing their variant rows. **GitHub Actions does not have private Google Drive source access:** the accompanying CI checks the frozen aggregate contract and unit-test behavior, **not** newly downloaded raw files.

## Why legacy posteriors do not prove gene causality

| Gene | Previously reported inherited GTEx-labeled ABF PP.H4 | Additional related evidence | Formal manuscript claim now |
|---|---:|---|---|
| **FGF5** | **0.7786969** | Same inherited-source multi-signal SuSiE PP.H4 **0.7522052**; source matched to 1,639 original table rows in this snapshot, but posterior **not re-computed here** | Historical exploratory evidence; no independent replication or causal claim |
| **ALDH2** | 0.0790399 | GWAS rs671 signal and EAS four-SNP full G0022 fine-map; strong local LD and coding/alcohol/BP pleiotropy alternatives | Variant-centric candidate; **not** an established gene-level MR effect, target direction or QTL-coloc |
| **SH3PXD2A** | 0.3508156 | Hypothesized cell-specific immune/vascular regulation | Candidate tissue/cell hypothesis only |
| **COL4A2** | 0.3307780 | Hypothesized cerebrovascular mural/stromal mechanism | Candidate tissue/cell hypothesis only |

Historical values above are **read from the existing IS scientific chronicle and master**, not computed using the newly read TSV files. Inherited 646 assays cover 43 genes × 16 tissues but are **not 646 independent experiments**; the primary historical label `GTEx_V8` lacks verified original fastQTL provenance and some scalar sample N and prior assumptions remain. Independent MR instruments / pleiotropy controls, causal mediation, sample independence, and cohort-matched QTL LD are unverified.

**Separate independent QTL *source*, not an independent stroke replication:** 2026-10-11 reprocessed EBI QTD000161 × BBJ original full-cis *exploratory* `coloc.abf` FGF5 H4 **0.9132** at assumed QTL `N=175`, `p12=1e-5`; its 15 tested N/prior scenarios span **0.370–0.991**. Source model, SNP universe and N assumptions differ from inherited H4 0.7787. Controlled same-source nested-vs-full SNP comparison reported **0.9260 vs 0.9132** at the same stated baseline, but this is not biologically independent. Neither model proves a single shared causal variant.

The **new 2026-10-18 P3A** gate on a *different archived EAS5/BBJ/GIGASTROKE summary source* had **36 candidates**, **25 source-stat-positive** and **0 independently approved or executed new P3A coloc analyses**. This is not a claim that **all historical 646 or 120 exploratory ABF runs were never executed**.

## Next scientific acceptance conditions (not yet satisfied)

1. **Original data model:** confirm exact archived `GTEx_V8` dataset accession, genetic build, donor ancestry, source QTL beta/effect allele orientation, source GWAS cohort N/case ratio, and genotype-level LD provenance. Compare each original per-variant GWAS/QTL effect and SE to authoritative sources using source hash and signed allele audit.
2. **Locus/statistical reproducibility:** verify palindromic/MAF conflicts and full cis region SNP coverage; read binary matrices' source hashes, SNP key order, symmetry/PSD and actual donor LD matching. Rerun legacy ABF / SuSiE with pinned R/coloc/susieR versions, priors, credible-set convergence and H3-vs-H4 sensitivity.
3. **MR and replication:** select truly independent cis instruments with allele direction/LD and IV strength, explicitly report any single-IV Wald-only limitation, handle ALDH2 rs671 alcohol/behavior/BP pleiotropy and possible QTL assay-binding artifacts; use independently sourced stroke outcome cohorts (GIGASTROKE EAS discovery contains BBJ and therefore is **not** standalone BBJ replication).
4. **Publication language:** use `historical/exploratory` and source-verified descriptive SNP claims only; never convert technical PASS, reference LD shape, high H4, or ordinary effect concordance to a proven causal gene or therapeutic target.

**Research decision:** FGF5 first for coloc provenance work because it has historical and separately modeled exploratory ABF results; ALDH2 for detailed coding/behavior/ancestry mechanism (without naïve MR), SH3PXD2A and COL4A2 remain cell-context alternatives. This is a priority recommendation, not new causal ranking.
