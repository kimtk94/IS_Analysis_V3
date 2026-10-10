# IS G0022 — Japanese vs whole-East-Asian LD, full-QTL availability and actionable research gates

**Research date:** 2026-10-10 KST
**Status:** GENOTYPE_LD_SENSITIVITY_VERIFIED / VALIDATED_STROKE_QTL_COLOCALIZATION_NOT_ESTABLISHED / MANUSCRIPT_CAUSAL_TARGET_BLOCKED
**Purpose:** Directly test whether EAS504 LD proxies are stable for a Japanese molecular QTL target and determine whether downloadable Japanese QTL datasets support unbiased full-cis colocalization.

## 1. Inputs & reproducibility constraints

Read-only:
- /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/reference_export.traw (5,060 genotype markers, 504 EAS participants)
- /srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/integrated_call_samples_v3.20130502.ALL.panel (populations: JPT104, CHB103, CHS105, CDX93, KHV99)
- Same G0022/GCST90104545 EAS AIS 4,649 variant source and 4,649² signed dosage-based LD matrix (stored matrix 172,905,608 bytes)
- rs671 chr12:112241766 G>A. PLINK .traw COUNTED allele matches REF; ALTERNATIVE dosage=2−COUNTED. Verified allele letters and file row order.

Original LD lead row had 19 markers with EAS r²≥0.50 to rs671 (including anchor). Only these 19 selected; no model re-fitting, no use of GWAS phenotype for subgroup computation.

## 2. Actual 1000G rs671-A reference allele frequencies

| Population in 1000G | N | ALT A frequency |
|---|---:|---:|
| EAS total | 504 | 0.173611 |
| Japanese in Tokyo (JPT) | 104 | 0.240385 |
| Han Chinese Beijing (CHB) | 103 | 0.160194 |
| Southern Han Chinese (CHS) | 105 | 0.271429 |
| Chinese Dai (CDX) | 93 | 0.043011 |
| Kinh Vietnamese (KHV) | 99 | 0.136364 |
| EAS ex-JPT | 400 | 0.156250 |

These frequencies come from the actual source dosage matrix for reference individuals, not from GIGASTROKE study samples. EAS/JPT contain the same JPT genotypes (nested, not independent).

## 3. Real LD r² to rs671: whole EAS vs Japanese subset

| 95% EAS AIS CS marker | Whole EAS n504 | JPT n104 | CHB n103 | CHS n105 | KHV n99 |
|---|---:|---:|---:|---:|---:|
| rs11066015 / ACAD10 (112168009) | 0.987528 | 1.000000 | 1.000000 | 0.959393 | 1.000000 |
| rs671 / ALDH2 (112241766) | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| NAA25 112468206 | 0.974237 | 1.000000 | 1.000000 | 0.912227 | 1.000000 |
| HECTD4 112736118 | 0.869614 | 0.920413 | 0.935947 | 0.736123 | 0.882109 |

CDX 93 r²=1.0 for these three proxy markers in this smaller subgroup; does not imply perfect LD for all populations or true genetic identity.

For 18 non-rs671 markers with valid JPT and EAS r² estimates:
- Median |JPT r² − EAS r²| = **0.0628233**
- Maximum |JPT r² − EAS r²| = **0.170661**
- Recomputed EAS n504 dosage LD vs original full 4,649² signed-LD matrix max absolute r² difference = **1.04×10⁻¹⁴** for selected 19 markers. This provides an exact row-order and dosage-contract check.

## 4. Interpretation

1. The 1000G EAS sample pools genetically differentiated EAS populations; rs671 ALT frequency varies substantially even within East Asia. Using pooled EAS LD as a **proxy for Japanese JCTF causal signal fine-map** risks misfit; JPT104 reduces population mismatch but increases finite-reference LD sampling error.
2. Japanese JPT104 r²≈1.00 among ALDH2/ACAD10/NAA25 markers (for this exact reference subset): Japanese LD alone CANNOT distinguish these variants. A larger Japanese source might yield precise but near-collinear effects, so do not assume causal resolution from PIP solely.
3. Subgroup LD estimates are NOT independent cohort replication or an alternative LD analysis of original GIGASTROKE study participants. No rerun SuSiE, credible set, QTL coloc, MR, drug directionality or mediator effects today.
4. The previous EUR AIS source audit showed 3,303/4,649 EAS SNPs matched overall but 0/4 four CS SNPs, 0/19 markers with pooled EAS rs671 r²>=.5 included. **This specific EUR source is not an immediate escape path for G0022 fine mapping**, even if other EUR genome associations might be useful elsewhere.

## 5. East-Asian molecular-QTL and LD source access gates (as of 2026-10-10)

| Source | Access / sample | Strength | Fatal limitation / next action | Gate |
|---|---|---|---|---|
| JCTF NBDC hum0343.v3 / NHA000193 | Unrestricted 1.7GB ZIP; 1,019 Japanese eQTL / 1,384 pQTL | Direct rs671 Japanese ALDH2 QTL association and PIP, measured Olink | Reports rows only if P<0.05 OR PIP>0.001. Genome-wide or gene cis *unfiltered* p and beta universe not released in this ZIP. Request complete tested cis-QTL association stats and SNP ascertainment list; Olink epitope/affinity artifact controls | SUPPORTED_MARGINAL_QTL, NOT_FULL_COLOC_ELIGIBLE |
| JCTF NBDC hum0343.v2 / NHA000172 | Unrestricted 749MB ZIP, Japanese eQTL/sQTL, 465 donors | Older Japanese RNA/sQTL cohort comparison | Need evaluate v2 source filtering and exact rs671 presence before a coloc claim; it is a smaller or overlapping JCTF predecessor, NOT independent external replication without cohort provenance | ALTERNATIVE_SOURCE_UNAUDITED |
| MAEEA Asian blood eQTL (2026) / Zenodo 21296030 | 2,024 EAS (1,005 Chinese + Japanese); 2.8m reported cis associations; source repository reported restricted in existing server metadata | Cross-Japanese/Chinese non-EUR molecular validation | Seek depositor approval for complete tested SNP–gene cis statistics, not just significant variants. No access to restricted file without authorization | PERMISSION_BLOCKED |
| NARD2 genotype reference (2023) | 14,393 East Asian WGS reference; public imputation website | Larger ancestry-aligned imputation source | Web imputation availability ≠ raw individual-level LD matrix availability. Need service compatibility with *LD extraction* or approved de-identified haplotype reference access. No claim that NARD2 source genotypes are locally available | RAW_LD_ACCESS_NOT_VERIFIED |
| 1000G EAS 504 / JPT 104 | LOCAL genotypes inspected, per-marker allele audited | Reproducible preliminary G0022 EAS/Japanese subgroup LD comparison | Both small relative to 4,649 markers, EAS ancestry mixture, LD rank≤503 and JPT≤103; no in-sample disease GWAS LD | EXPLORATORY_ONLY |

**Sources:**
- JCTF 2024 metadata https://humandbs-production.ddbj.nig.ac.jp/en/dataset/NHA000193
- JCTF 2022 metadata https://humandbs-production.ddbj.nig.ac.jp/en/dataset/NHA000172
- Japanese Wang 2024 https://pmc.ncbi.nlm.nih.gov/articles/PMC11525184/
- MAEEA Wang 2026 https://academic.oup.com/nsr/advance-article/doi/10.1093/nsr/nwag566/8779978
- MAEEA source record https://zenodo.org/records/21296030 (repository marked restricted in prior server audit)
- NARD2 14,393 https://pmc.ncbi.nlm.nih.gov/articles/PMC10411914/
- NARD2 SNU Korean institutional access explanation https://snurnd.snu.ac.kr/?q=node%2F1435

## 6. A complete East-Asian cis-QTL request contract (minimum needed)

Required output without p/PIP thresholding for **every tested SNP–molecular phenotype association** in at least the G0022 region:
- Variant ID with GRCh37 or GRCh38 build, chromosome, base-pair position, REF/ALT, effect allele, allele-frequency and per-variant INFO/missingness/exclusion indicator
- Gene Ensembl/version and tissue; molecular phenotype ID, platform (RNA, Olink probe ID, peptide LC-MS if available)
- Per-variant beta, SE, nominal p-value, exact N for that QTL, sample covariates/normalization and ascertainment strategy
- Availability of matched signed LD/reference genotype or finely mapped independent molecular credible sets, per-signal locus boundaries, local ancestry/relatedness
- Source provenance for genotype/phenotype sample overlap vs BBJ and GIGASTROKE; sex/medical condition (COVID infection), protein epitopes and imputed vs measured variants.
- Data use permission and conditions; no unauthorized restricted-data access.

**If unfiltered source unavailable:** fine-mapping/coloc gate remains BLOCKED. JCTF filtered summary is permissible as descriptive variant-level supporting evidence, not a new positive causal-gene colocalization.

## 7. Next predeclared analysis decision

**Priority 1:** Determine whether the 2022 JCTF v2 releases *all* cis-SNP–gene associations or filtered significant rows only. Do not auto-download 749MB just to produce posterior unless source columns/filter pass first.

**Priority 2:** Contact MAEEA and JCTF custodians for unfiltered G0022 cis source and ask NARD2 for a *scientifically permitted* signed LD/correlation export, not just imputation.

**Priority 3:** Fine-map EAS AIS with LD finite-reference sensitivity, compare all five 1KG EAS population subgroups, and obtain independent *AIS* cohort validation without BBJ reuse.

**Manuscript gate:** At present, fine-map PIP and Japan ALDH2 molecular marginal associations are independently documented, but shared causal variant, protein-mediated disease effect and therapeutic inhibition/activation remain **NOT_ESTABLISHED**.

## 8. Machine-readable reproduction

Code: scripts/is/audit_g0022_jpt_eas_ld_stability.py

Output: /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/jpt_eas_ld_sensitivity_v1
- G0022_RS671_SUBGROUP_LD_SENSITIVITY.tsv
- G0022_RS671_SUBGROUP_LD_SENSITIVITY_SUMMARY.json

Run:
    python3 scripts/is/audit_g0022_jpt_eas_ld_stability.py --traw /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/reference_export.traw --panel /srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/integrated_call_samples_v3.20130502.ALL.panel --variants /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/variants.tsv --eas-signed-ld /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/ld.rowmajor.f64 --outdir /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/jpt_eas_ld_sensitivity_v1

No primary data, canonical GWAS, other research projects, or LIVE services modified.
