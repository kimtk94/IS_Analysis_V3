# Human ischemic-stroke expression datasets as disease-context alternatives to ROSMAP (2026-11-11)

**Research method note / not a course submission.** The distinction is essential: human IS scRNA-seq contains *disease-associated expression states*, whereas ROSMAP provides genotype-linked brain cell-type eQTL. Human IS scRNA-seq does **NOT** replace ROSMAP QTL statistical function unless paired donors' genotype and full cis-eQTL summary data are actually available.

## Highest priority clinical source: GSE331114

- **Official GEO:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE331114 ; deposited 2026-05-20; Homo sapiens; 8 brain scRNA-seq samples, **AIS n=4 and non-stroke controls n=4**.
- 10x Genomics processed per-sample matrix files in GEO, raw FASTQ in SRA. Individual GSM metadata report **CD45+ immune-cell enrichment/isolation** from human brain parenchyma. This is **not** an unbiased whole brain atlas of neurons/astrocytes.
- NCBI study summary mentions skull marrow AIS 7/control 7 and brain AIS 4/control 4; GEO series **GSE331114 has 8 brain samples**. Never claim 22 samples are present in this GEO series.
- Genome alignment reported **GRCh38**, while main GIGASTROKE EUR/EAS GWAS files are **GRCh37**. Expression gene symbols can be mapped with stable Ensembl IDs, but variant-level matching (not provided) requires provenance.
- Primary biological context: infiltrating/activated T cells, myeloid inflammatory response in post-ischemic brain. **Do not** infer baseline pre-stroke genetic susceptibility or cell-specific eQTL from case-control differential expression.
- Risk: cohort size n=4+4 and possible brain tissue acquisition/clinical and control differences; use pseudobulk by donor, not cells as independent n. No genotype paired to tissue reported in source metadata.

## Second-tier sources

| Source | Disease/material | Design | Re-analysis role | Main limitation |
|---|---|---|---|---|
| [GSE267023](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267023) | Human IS PBMC | 16 sample records; acute and chronic phases per series description | time/stage-specific immune expression | no proven independent control sample count; donors/time points must be deduplicated |
| [GSE285659](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285659) | Human large-vessel occlusion intravascular PBMC | ischemic hemisphere intravascular vs systemic arterial control **within stroke patients** | arterial compartment contrast | not independent stroke cases vs healthy controls, multiple GEX/ADT libraries per donor |
| [Amini et al., 2020](https://doi.org/10.1002/acn3.51154) | Blood DNA genotype + mRNA measured in 137 acute IS patients / 138 vascular-risk matched controls | SNP × stroke diagnosis interaction eQTL | **most directly stroke-conditioned QTL**, pending accessible summary SNP–gene coverage | post-stroke response, not disease susceptibility; report top cis/trans interactions but full QTL matrix file not shown |
| [Human peri-infarct cortex 2026](https://doi.org/10.1016/j.expneurol.2025.115550) | Human peri-infarct cortex 2 IS / 3 trauma surgical controls, 23,638 nuclei | snRNA-seq | potential glial/oligodendrocyte context if data acquired | raw/data/code available **upon request**; very small sample; comparator trauma tissue |
| [GSE174574](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174574) | **Mouse** MCAO brain n=3 vs sham 3 | scRNA-seq | translational mouse sensitivity only | not human or human QTL |
| [GSE337349](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE337349) | **Mouse**, 21-source integrated atlas 103 samples | reprocessed cross-study scRNA | complementary injury time course | mouse and different models, no human human-stroke genetics |
| [STARNET / Athero-Express](https://www.nature.com/articles/s44161-025-00628-y) | human vascular/carotid plaque genetics and expression | candidate vascular eQTL for LAA mechanism | more proximal vascular etiology than AD brain | tissue is atherosclerotic plaque not infarct brain; controlled/dataset rights to review |

## Inference layers and replacement decisions

- **Tier 1, ischemic stroke response:** GSE331114 (AIS human brain CD45+ immune) → disease-specific *immune cell expression* validation. GSE267023 / GSE285659 → peripheral vascular/immune sensitivity.
- **Tier 2, genotype–RNA regulation:** GTEx v8 and/or ROSMAP / PsychAD brain cell-type eQTL retain proper variant-level locus-to-gene statistical role. Amini 2020 stroke-conditioned interaction eQTL may add disease-state evidence **only if full variant-level summary statistics become available and correctly harmonized**.
- **Tier 3, ancestry:** EAS eQTL (MAEEA or genuinely complete alternative) independently required. Stroke patients recruited in China do **not** establish genetically confirmed EAS ancestry unless donor genotypes or reported ancestry are verified.
- **Tier 4, optional annotation:** MAGMA gene-based + gene-property and LDSC-SEG may identify gene/tissue-wide background; same GWAS inputs mean no independent replication.

**Design correction:** Replace unqualified "ROSMAP substitute" with **"stroke-specific immune transcriptional validation replacing ROSMAP expression-context only"**. Keep ROSMAP or equivalent brain eQTL for risk-locus colocalization. Selection is **conditional GO**, requiring real sample/raw count access QC, counts/labels, cell state mapping, pseudobulk at donor level and comparison against already published interpretations.

## Gates

1. **S0:** Official GEO SOFT metadata confirms GSE331114 eight true brain samples, 4/4 labels, processed count links, and whether donors are independent. For the others resolve number of unique donors vs library sample records.
2. **S1:** Download sample-level **features + barcodes metadata** first, then a **small single matrix** with checksum and shape. Do not announce whole-data completion from FTP listing.
3. **S2:** Define donor pseudobulk and multiple-testing strategy with min-cell criteria, classify CD45+ immune composition bias and possible severe-stroke surgery selection bias.
4. **S3:** Test whether GWAS/coloc-prioritized genes show detectable expression and group-related shifts; distinguish **descriptive expression response** from genotype-regulation and susceptibility.
5. **S4:** For full stroke-conditioned eQTL fallback, verify 2020 study raw genotype/SNP–gene summary access; if absent, use its reported results only as literature context.

### Novelty warning

Published 2024 and 2026 studies already perform stroke GWAS × brain cell eQTL with colocalization and single-cell expression (e.g., DOI 10.1016/j.compbiolchem.2024.108086; PMID 42003665). Neither addition of a 2026 stroke expression dataset nor cell-type annotation alone ensures new scientific contribution. Compare exact gene, locus and tissue claims against published work. Use study evidence to critique tissue appropriateness, not to fabricate causal evidence.

### Automated verification

`scripts/is/audit_stroke_versus_rosmap_sources.py` checks public **official GEO family SOFT** sample groups and processed-file metadata only.
`tests/is/test_stroke_versus_rosmap_sources.py` tests sample/format logic.
`.github/workflows/is_stroke_alternatives_audit.yml` records a status artifact. No patient-level genotypes or restricted material downloaded.
