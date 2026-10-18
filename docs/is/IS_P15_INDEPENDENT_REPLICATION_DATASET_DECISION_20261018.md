# IS P1.5 — Independent East Asian GWAS replication dataset decision

Updated: 2026-10-18
Status: **No BBJ-excluded EAS-AIS public summary-statistics release verified yet; do not claim independent replication.**

## Established facts with primary sources

1. GIGASTROKE primary data availability lists GWAS Catalog accessions **GCST90104534–GCST90104563**. This documents released ancestry/phenotype meta-analyses, **not** a confirmed public BBJ-excluded EAS-AIS meta-analysis. https://www.nature.com/articles/s41586-022-05165-3
2. Publisher Supplementary Table 1 EAS discovery contains **BIOBANK_JAPAN** AIS **17,651 cases**, controls **159,610**. Publisher workbook SHA256 **7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4**; verification workflow https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39126549153
3. GIGASTROKE explicitly states that for its **East Asian iPGS training and testing** the separate BBJ training set (577 AIS/9,232 control) and evaluation set (1,470 AIS/40,459 control) **were not included in GIGASTROKE discovery meta-analysis**. This PGS split does NOT negate BBJ discovery participation; no valid replication inference from these PGS cohorts without access to matched GWAS results. https://pmc.ncbi.nlm.nih.gov/articles/PMC9524349/
4. Another public BBJ ischemic stroke study **bbj-a-129**, reported as Ishigaki et al. 2019, has 17,671 cases / 192,383 controls / N=210,054, GRCh37. These numbers **differ** from GIGASTROKE BBJ discovery 17,651 / 159,610, thus cohort overlap cannot be quantified from sample-size counts and labels alone. No independent study proof. https://gwas.mrcieu.ac.uk/datasets/bbj-a-129/
5. KoGES public Korean GWAS resources exist, including the original 76 trait summary statistics (2022, https://zenodo.org/record/7042518) and another public 2025 Korean GWAS source (https://zenodo.org/records/15132424). **Ischemic-stroke-specific phenotype, case counts and GWAS availability not verified**. The original paper involved some **KoGES/BBJ meta-analyses**, so each usable trait file requires explicit source/phenotype provenance. Do not classify these as an independent stroke GWAS without inspection. https://www.sciencedirect.com/science/article/pii/S2666979X22001380 and https://doi.org/10.1038/s41467-025-59950-5

## Dataset eligibility gate

| Dataset | Role | Current status | Next action |
|---|---|---|---|
| BBJ original IS GWAS | Discovery anchor | Confirm release accession, effect AF, source cohorts | Anchor, **not** independent validation |
| GIGASTROKE EAS AIS published meta-analysis | Descriptive concordance / shared-source sensitivity | **NONINDEPENDENCE_RISK_CONFIRM_BBJ_IN_DISCOVERY** | Never call independent replication against BBJ |
| GIGASTROKE EAS AIS BBJ-excluded leave-one-biobank-out | Preferred replication candidate | **NOT_VERIFIED_PUBLICLY_AVAILABLE** | Search accessions / request summary stats from authors with documentation of excluded BBJ participants |
| GIGASTROKE external EAS follow-up | Possible validation candidate | Cohort identities and downloadable variant-level summary stats **PENDING** | Distinguish follow-up results from discovery or PGS model evaluation |
| KoGES / other Korean stroke GWAS | Potential alternative | **STROKE_TRAIT_AND_ACCESS_NOT_VERIFIED** | Audit phenotype list and constituent studies before adoption |
| BBJ 2019 vs BBJ 2021 | Same-source sensitivity only | **POTENTIAL_COHORT_OVERLAP** | Not external replication |

## Required analysis once independent summary statistics obtained

- Identify cohort/sample exclusion with supporting methods or signed author confirmation.
- Confirm phenotype definition (AIS, AS, subtypes), ancestry EAS, case/control N and overlap.
- Enforce GRCh37/GRCh38 liftOver lineage, SNP alleles/AF/strand, beta/SE/effect units; exclude ambiguous palindromes unless AF disambiguates.
- Compare *independent loci* rather than treating LD-correlated SNPs or multiple mapped genes as independent observations.
- Descriptive effect concordance may proceed despite source overlap but never label it replication.
- Freeze release ID, timestamp, SHA256, exact cohort inclusion roster, licenses, N, genome build, and pipeline commit.
- No gene-level causal claim or formal QTL colocalization without matched tested-region coverage, QTL donor N, and appropriate LD.

## Decision

P1.5 independently replicated genes: **0 established**.
Continue original 2,225-gene universe and existing exploratory candidate workflow, but mark external GWAS replication **BLOCKED_NEEDS_UNSHARED_COHORT_DATA**. Existing BBJ vs GIGASTROKE signals remain descriptive and potentially valuable for within-discovery evidence triage.
