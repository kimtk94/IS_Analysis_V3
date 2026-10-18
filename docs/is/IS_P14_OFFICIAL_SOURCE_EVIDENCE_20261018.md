# IS P1.4 — Official-source provenance and unresolved evidence

Date: 2026-10-18
Status: **DOCUMENTARY EVIDENCE FOUND; INDEPENDENT REPLICATION NOT ESTABLISHED**

## Verified sources

### GIGASTROKE 2022 (Nature)
- Paper: https://www.nature.com/articles/s41586-022-05165-3
- Consortium discovery meta-analysis: 29 population-based cohorts or biobanks plus 25 clinic-based case-control studies.
- Overall up-to counts: 110,182 stroke cases and 1,503,898 controls; AIS cases = 86,668.
- **Do not treat the above as the ancestry-specific EAS AIS sample size.** Obtain study/phenotype/ancestry-specific counts in supplementary tables.
- Distinguish consortium discovery, internal cross-ancestry follow-up, and independent external replication.
- BBJ also appears in East Asian polygenic score evaluation. This does not itself establish BBJ's inclusion in the particular GIGASTROKE EAS AIS discovery file, nor prove cohort non-overlap. Require cohort-level lineages and sample IDs before independent replication claims.
- Relevant source documents: GIGASTROKE Supplementary Table 1 (cohorts), Table 2 (genotyping), and per-ancestry/per-phenotype summary statistics catalog GCST90104534–GCST90104563.

### eQTL Catalogue GTEx v8
- Official studies: https://www.ebi.ac.uk/eqtl/Studies/
- QTS000015 represents reprocessed GTEx v8 across 49 tissues, total 15,178 samples and 838 donors (source studies table).
- eQTL dataset QTD000161 is associated with this GTEx source according to public catalogue mirror: https://filer2.niagads.org/browse_faceted_with_search.php?assayType%5B%5D=eQTL&count=20&genomeBuild%5B%5D=hg38&start=580
- **838 is not the tissue-specific effective QTD000161 eQTL N.** Verify precise tissue label, number of donors with paired genotype/expression data, analysis version, tested SNP universe, effect AF semantics and genotype LD source before recomputing coloc or SuSiE.
- Distinguish EAS/JPT public reference-panel LD from genotype LD calculated on GTEx QTL study donors; do not equate them.

### BBJ Japanese GWAS
- Original population-scale BBJ disease GWAS: https://pmc.ncbi.nlm.nih.gov/articles/PMC7968075/
- A shared biobank label is not a substitute for cohort-level overlap attestation. Confirm BBJ GWAS source exact version, sample inclusions and original supplementary materials.

## P0/P1 observed machine-verified facts

- P1 full universe 2,225 gene IDs / 2,425 gene-region pairs; provenance run https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39124908298
- P1.4 metadata audit https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39125943738
- P0 four-gene comparators CALHM2, FGF5, NEURL, C4orf22 lack verified cohort independence; QTL study N / matched LD are not established in P0 scorecard.
- Existing source-variant direction checks are **descriptive**. Same allele-aligned sign among correlated SNPs is not independent genetic replication.
- Formal multisignal colocalization stays BLOCKED_QTL_STUDY_LD_AND_N.

## Required independent validation gates
1. Extract GIGASTROKE EAS AIS cohort roster, BBJ participation and leave-one-cohort-out/source sample overlap evidence; do not infer from PGS evaluation subset.
2. Extract BBJ 2021/target release exact phenotype sample counts, controls and lineage.
3. Map QTD000161 exact tissue, donor N and QTL model metadata from authoritative catalogue accessions.
4. Separate GWAS cohort-matched LD, 1000G EAS/JPT reference LD, and QTL donor-derived LD.
5. Only after independent data provenance is affirmed: allele harmonization (strand/palindromic AF), effect direction, overlap and appropriate region-level colocalization.

No causal claims. No changes to raw data or LIVE config.


## 2026-10-18 addendum — direct evidence from GIGASTROKE original full text

Original primary paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC9524349/

- GIGASTROKE explicitly defines **internal across-ancestry discovery-sample follow-up** versus **external replication**. The latter totals 89,084 stroke cases (85,546 AIS) and 1,013,843 controls; 10.1% of stroke cases in the external follow-up are East Asian. This is NOT the case/control N of EAS AIS discovery GWAS.
- BBJ is explicitly named as an East Asian **iPGS predictive evaluation** cohort: 1,470 prevalent AIS cases and 40,459 controls; the iPGS model derivation set is 577 cases / 9,232 controls. These are predictive-model sample counts, NOT confirmed GIGASTROKE EAS AIS discovery N.
- Independent BBJ2 first/second-cohort context: https://pmc.ncbi.nlm.nih.gov/articles/PMC11451717/ describes BBJ1 PRS derivation and independent BBJ2 evaluation at 1,470 cases / 40,459 controls. Independence of BBJ2 from BBJ1 in this PRS task does **not** establish independence of every BBJ discovery GWAS release from GIGASTROKE.
- The paper's associated supplementary data include a downloadable 56-table XLSX, and its Methods point to study-specific recruitment for exact discovery and replication cohorts. It must be crosswalked to the **actual GIGASTROKE EAS AIS summary-statistics accession** and the exact BBJ IS release before the cohort-independence gate can be marked PASS.
- GIGASTROKE independent external follow-up is reported for stroke associations; it should not be reinterpreted as independent replication of our specific QTL/gene-level mediation claims.

## Frozen conservative gate interpretation

- P1.4 cohort overlap: **PENDING_DATASET_SPECIFIC_CROSSWALK**, not automatically FAIL because BBJ appears in a predictive risk score assessment.
- EAS AIS discovery study N: **PENDING_SUPPLEMENTARY_TABLE_EXTRACTION**.
- BBJ 2021 IS release overlap with GIGASTROKE: **PENDING_COHORT_AND_RELEASE_CROSSWALK**.
- QTD000161 study-specific N / donor LD: **PENDING_PRIMARY_QTL_METADATA**.
- Multisignal coloc: **BLOCKED_QTL_STUDY_LD_AND_N**.


## 2026-10-18 direct supplemental workbook verification — BBJ in GIGASTROKE EAS discovery

**Official source file**: https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05165-3/MediaObjects/41586_2022_5165_MOESM4_ESM.xlsx

- SHA256: `7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4`
- Reproducible extractor: `scripts/is/audit_gigastroke_supplement_cohorts.py`
- Successful Actions run: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39126549153
- `Suppl.Table1`, spreadsheet row **62**, ancestry EAS, study `BIOBANK_JAPAN*`:
  - Incident & prevalent combined: **AS = 17,651; AIS = 17,651; CES = 758; SVS = 4,915; LAS = 1,331; controls = 159,610**.
  - The row is under the original discovery study sample section, distinct from the separate follow-up section.
- `Suppl.Table2`, spreadsheet row **78**, study `BIOBANK_JAPAN_1`:
  - genotyping Illumina HumanOmniExpressExome / HumanOmniExpress + HumanExome; imputation reference BBJ1k+1KGp3v5; REGENIE v2.0.2, approximate Firth logistic regression.
  - Rows 79 and 80 list `BIOBANK_JAPAN_1r (for PGS)` and `BIOBANK_JAPAN_2 (for PGS)` **separately**. PGS datasets are not interchangeable with GWAS discovery data.

**Updated interpretation**:
1. **CONFIRMED: a BBJ cohort participated in GIGASTROKE EAS discovery**. Therefore a naïve analysis of BBJ discovery GWAS vs GIGASTROKE EAS discovery GWAS must NOT be labelled independent replication.
2. **UNRESOLVED: identity and participant overlap** between the exact BBJ 2021 IS release under study and the GIGASTROKE BBJ1 17,651/159,610 cohort. A same-named biobank does not by itself prove that all individuals overlap, and absence of matched sample IDs prevents quantifying overlap.
3. `BBJ_vs_GIGASTROKE_EAS` becomes `NONINDEPENDENT_COMPARISON_RISK_BBJ_PRESENT_IN_META` for independent-replication gating (conservative), rather than generic `PENDING`. **Do not use this comparison as independent replication evidence.**
4. Preserve descriptive allele-aligned effect-direction analyses; note their status as internally overlapping-source comparisons, not external validation.
5. Replication route: use GIGASTROKE EAS **leave-BBJ-out** meta-analysis with documented exclusion, or independent externally sampled EAS stroke GWAS, after exact cohort provenance is verified.

This supersedes earlier uncertainty in this document about whether BBJ contributes to the GIGASTROKE discovery meta-analysis. It does **not** establish full case-control overlap with the specific BBJ 2021 release or the QTL coloc assumptions.
