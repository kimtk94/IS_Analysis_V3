# IS Phase10/11 — G1/G2 full 646-case audit (2026-10-11)

**Verified: full 646/646 source-level MIN/MEDIAN/FIRST QTL scalar-N sensitivity; full 646/646 G1 input structural allele consistency. Neither proves a causal gene.**

## Source provenance

Actual [Google Drive full output folder](https://drive.google.com/drive/folders/1kJJl2PV0aJsfnrRQx3OORrXv49p7pHd4), [646-row N table](https://drive.google.com/file/d/16BOYFZ0pNKRu74pZmyqaXDfRFPSAd7dJ/view), [full N summary](https://drive.google.com/file/d/1mxG0p6l3rplWu5W5yQinYbOPGevD7puI/view), [full 646 source SHA256s](https://drive.google.com/file/d/1ZwnZ4XsCeo2ldMys9gVRDmYaEDBPZDwi/view), [execution manifest](https://drive.google.com/file/d/11lbLEb4nADjMLpNfTLYkmYtOhpMg1-xr/view), [646 G1 statuses](https://drive.google.com/file/d/19l11e4nPWzAbHRqZ_MSubWr8li3QEMh4/view) and [G1 summary](https://drive.google.com/file/d/11djVyV0cFG7yOpAFeXigbukiLXQSIb4u/view).

Independently cross-checked all 646 unique locus/tissue/gene keys, input SNP counts, source SHA256 entries, baseline H3/H4 and scalar N against original 646-G0 replay, original index and raw master: **all 646 agree**. New full N status TSV SHA256 `e118206ad83b70357aefaff700456d6d4292fc6455d6900f77f6ad2ad8c85ddf` matches Colab manifest. R script SHA256 `6f759ae76d25b00c4f308b0046ff03d27d668c2019f611fe446664b6819f796f`; G1 source script SHA256 `a7bca1a61b48bbe63770d04530f1158aa2e449c1ff2323ab86bb4240928cd650`. Pinned coloc 5.2.3, p1=p2=1e-4, p12=1e-5; 646 PASS, 0 FAIL, 0 MISSING both for numeric scalar N and internal G1 QC.

## G2 full original 646 scalar-N sensitivity

All 646 input assays satisfy **FIRST = MEDIAN = MAX QTL N**, with MIN N below maximum by 1.95–20.63%, but only median **0.737%** of SNP rows per assay have N lower than maximum. Thus forcing all variants to MIN N is an extreme diagnostic, **not a valid SNP-specific effective-N model**.

| H4 criterion | FIRST/MAX baseline | MIN scalar QTL N |
|---|---:|---:|
| H4 >= 0.50 | 6 | 5 |
| H4 >= 0.75 | 2 | 2 |
| H4 >= 0.80 | 0 | 0 |

**INA** GTEx anterior cingulate BA24 is the only H4>=0.5 test crossing below the descriptive threshold: **0.5363685 → 0.4933763**. CALHM2 brain cerebellar hemisphere **0.7976836 → 0.7666429**; FGF5 **0.7787 → 0.7770**. Across 646 assays **535 H4 values increase**, **111 decrease**, median absolute change **0.00141486764** and maximum absolute change **0.04299228949**.

## G1 structural QC

646/646 full original gene–tissue SNP input files passed internal GRCh38 REF/ALT/variant identity and GTEx beta/SE/MAF preservation plus GWAS EAF↔MAF arithmetic tests. **1,224,224 repeated gene–SNP rows**, including 164,317 palindromic records (16,791 GWAS EAF 0.4–0.6). Inter-population absolute GWAS–GTEx MAF differences >0.1 for 508,318 rows and >0.2 for 199,221. These are **not unique SNP counts or proven effect-allele flips**.

A previous independent BBJ native Allele2/BETA/AF and GRCh37→38 chain check covered **8,309 distinct GWAS SNPs from 30 cached gene–tissue tests**, not all variants of the 646-source collection.

## Manuscript-release boundaries

- G0 numeric replay **COMPLETE 646/646**, historical **four original BBJ GWAS loci**, 43 genes and 16 GTEx tissues; **not** the broader 80 provisional GWAS windows / 2,225 positional genes.
- G1 internal source QC **COMPLETE 646/646**, but full unique-SNP native GWAS allele/chain/FASTA and independent GTEx source transform checks **PARTIAL**.
- G2 scalar MIN/MEDIAN/FIRST robustness **COMPLETE 646/646**, but genuine per-variant effective-N model, GTEx-matched genotype LD, multi-signal colocalization and gene/tissue multiplicity **PENDING**.
- G3 stroke-specific cellular regulatory mechanism and independent molecular validation **PENDING**.
- **No causal gene or protein mediation claim is established.**

Historical [30-source pilot](IS_PHASE10_11_G1_G2_QTL_N_SENSITIVITY_20261011.md) is retained. Independent plot and full audit JSON/Markdown are available as conversation artifacts.
