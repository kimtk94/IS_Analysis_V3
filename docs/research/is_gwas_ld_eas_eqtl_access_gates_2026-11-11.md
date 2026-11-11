# IS GWAS–LD–eQTL source and statistical decision gates (2026-11-11)

**Purpose:** internal source-audited master's study design, not ready-to-submit course text. The Research Proposal class template prohibits generative-AI-authored submissions.

## Current main research question

Among human ischemic stroke (IS) risk loci, how does the inferred locus-to-gene regulatory assignment change between EUR- and EAS-ancestry-matched cis-eQTL, with ancestry-matched LD, harmonization, source ascertainment, and equivalent tested cis variants?

This is more specific than merely noting East Asians are underrepresented. The Sept 2026 MAEEA study (DOI [10.1093/nsr/nwag566](https://doi.org/10.1093/nsr/nwag566)) already demonstrates ancestry-matched eQTL advantages for broad traits, so stroke locus-level **incremental novelty must be shown**.

## Full original GIGASTROKE outcomes (already independent checksum QC PASS)

| Ancestry | Accession | Cases / controls (published) | Original genome build | Full-file QC |
|---|---|---|---|---|
| EUR IS | GCST90104540 | 62,100 / 1,234,808 | GRCh37 | original whole gzip CRC, per-SNP structure and publisher MD5 MATCH |
| EAS IS | GCST90104545 | 19,032 / 237,242 | GRCh37 | original whole gzip CRC, per-SNP structure and publisher MD5 MATCH |

Source [complete original structural audit](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43003883431) and [independent publisher MD5 check](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43005166208). Standardized columns include CHR, POS, EA, OA, beta, SE, P, EAF; neither `rsID` nor per-SNP sample size N appears in observed columns. Clinical/TOAST versus administrative ICD and per-subcohort case contribution remain partially unresolved.

## Reference LD: exact new implementation

Source [Gazal S-LDSC reference Zenodo record 7768714](https://zenodo.org/records/7768714):

- EUR `1000G_Phase3_plinkfiles.tgz`; provider MD5 `a7773ab485827b533cb300c76356d76b`, ~288MB.
- EAS `1000G_Phase3_EAS_plinkfiles.tgz`; provider MD5 `abf52fba6416622ea0757ca9aca51a87`, ~323MB.
- These source archives are from the S-LDSC EUR/EAS 1000G Phase3 reference set; **source build GRCh37 is documented/assumed, not independently proven by a full variant reference-coordinate liftover**.
- The CI [original GWAS vs 1000G chr1+2 pilot](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43021975612) transfers original complete compressed GWAS and LD reference files, verifies both publisher MD5 values, parses PLINK BIM coordinate allele pairs and scans entire GWAS while matching chromosomes 1 and 2 only. This pilot is **not** whole-genome matching, nor does it compute LD matrices.
- SNP classes are separated as (a) same chrom/pos nonpalindromic unordered allele pair, (b) palindromic A/T or C/G, (c) complement strand but unresolved direction, (d) multiallelic reference locus, (e) same coordinate discordant alleles, (f) absent reference coordinate. Never recode a class as allele-harmonized effect orientation.
- GWAS source full-file no-variant-level-N remains statistical-model limitation; case/control-derived effective N is an illustrative study-level approximation, **not per-variant N**.

## EAS cis-eQTL source prioritization

| Resource | Disease/tissue/ancestry | Actual access | Statistical caveat | Decision |
|---|---|---|---|---|
| [MAEEA/CBEA whole blood](https://doi.org/10.5281/zenodo.21296030) (2026 DOI 10.1093/nsr/nwag566; 2,024 EAS) | EAS whole blood meta-eQTL | Zenodo record metadata fetched but **0 attached files in previously audited API result**; source description claims all cis-SNP–gene summary stats | Need effect, SE, EAF, alleles, tested all cis variants, N, ancestry/build and LD; no completed raw download | **P0, SOURCE_NOT_YET_ACCESSIBLE (HOLD)** |
| [NHA000193 JCTF](https://humandbs-production.ddbj.nig.ac.jp/en/dataset/NHA000193), 2024 | Japanese COVID-19 blood eQTL; 1,019 eQTL donors (1,405 overall study cases) | Public unrestricted 1.7GB QTL ZIP, readme | **Selected pairs only** P<0.05 or PIP>0.001, expression conditional on infection, so **do not run standard full-cis coloc naively** | Sensitivity/QTL annotations ONLY |
| [NHA000051 Japanese immune-cell eQTL](https://humandbs-production.ddbj.nig.ac.jp/en/dataset/NHA000051), 2017 | 105 healthy Japanese donors, five immune cell types and whole blood | Public unrestricted 1.8GB ZIP, GRCh37 genotype basis | Verify whether per-gene nominal *all tested* cis variants are actually included; n=105 gives weak power, potential imbalance 21M/84F. Cannot label full-coloc ready without download and row-level QC | **P1 feasibility / small-cohort sensitivity, NOT yet coloc** |
| [Chinese Immune Multi-Omics Atlas (CIMA)](https://ngdc.cncb.ac.cn/omix/release/OMIX009786), 2026 | 428 Chinese healthy adult PBMCs; 70 cell contexts | **Controlled-access** 39.49GB cis-eQTL summary archive, requires data access committee | Even with access, patient population/cell context differs and QTL sample power varies by cell type | **ACCESS_HOLD** |
| [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/Data_access/) | EUR tissue/cell context | Full region query and downloadable `all.tsv.gz` available for selected datasets (do not assume *every* study retains all nominal associations) | Need actual dataset version, genome build, allele conventions, variant mapping; choose EUR comparator matched tissue | **P0 source metadata pilot** |
| [GTEx v8](https://www.gtexportal.org/home/v8ReleasePage) | chiefly EUR normal multi-tissue eQTL | Median TPM downloaded/QC; *all pairs* GTEx GCS are requester-pays for full association files | Normal TPM ≠ eQTL; full cis associations may incur GCS project access/egress cost | Expression context GO; eQTL-source cost HOLD |

Other results: [GTEx actual 54 tissue median TPM](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43010162288) and [human stroke brain 8/8 CD45+ immune single-cell data](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43010007530) are **biological expression context only**, not variant–expression effects in genotyped patients. No new eQTL colocalization or independent EAS biological association inferred from them.

## Statistical analysis gates

**G0 — Phenotype ascertainment (PARTIAL_HOLD):** preserve any ischemic stroke vs cerebral infarction/TOAST/lacunar cases and distinguish clinically adjudicated vs ICD electronic medical record cohort definitions. Do not assume identical EAS vs EUR phenotype composition.

**G1 — GWAS provenance and SNP/LD (IN_PROGRESS):** official GWAS MD5/CRC PASS; verify 1000G EUR/EAS source MD5, actual BIM allele overlap, GRCh37 coordinate matching and strand/ambiguous variants by ancestry. If chr1-2 passes, expand to 22 autosomes; record coverage per chromosome and compare selected loci before coloc.

**G2 — Full cis-eQTL access (HOLD):** record actual physical access, full cis variant testing (including null associations), variant effect orientation, beta/SE or adequate P+MAF+N, genome build and gene ID; candidate-only thresholded eQTL not sufficient for standard coloc. A locus lacking a complete all-variant summary remains UNTESTABLE.

**G3 — Fine-map and coloc (BLOCKED):** avoid single causal variant assumption where multi-signal. Use ancestry-specific LD and conditional fine-mapping/SuSiE-coloc when feasible, predefine priors and locus/variant intersection. Do not interpret PP4 or spatial expression alone as causal proof.

**G4 — MAGMA gene-based / gene-property and LDSC-SEG (BLOCKED):** only after correct variants, ancestry reference LD, N strategy and LDSC genome-wide h2 intercept/ratio QC. LDSC-SEG estimates tissue/cell gene-specific heritability enrichment, not gene causal identity. A source reused across tests is not an independent replication.

**G5 — Novelty/final decision (HOLD):** compare prior stroke single-cell/brain eQTL coloc (2024), EAS-specific proteogenomic MR (2025), MAEEA ancestry-matching results (2026); assess whether any robust, new locus-to-gene changes justify a master's thesis. With no full EAS cis-eQTL, pivot to EUR validated stroke-regulatory study with EAS GWAS directional replication, explicitly **NOT cross-ancestry eQTL portability**.

The class template prohibits use of generative AI in the submitted Research Proposal: this document is an internal evidence ledger and execution specification only.
