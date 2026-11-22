# IS P3A — primary expression-matrix public-source inventory and GWAS beta-unit evidence

**Date:** 2026-11-23 KST. **Project:** MASTER DEGREE ischemic stroke, 36 authentic gene×tissue QTL assays, 20 original genes, one chr10 G0015 locus. **Formal coloc 0/36 eligible or executed; scientific HOLD remains.**

## Verified original source audit

[GitHub Actions original-source CI #45040743755](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45040743755) **SUCCESS**: 11/11 focused tests plus full repository smoke testing **390 tests, 389 PASS and 1 unrelated Rscript SKIP**. The workflow loads prior immutable original summary evidence (QTD000131/QTD000171/QTD000434, 148,344 nonpal canonical SNP×assay records, all 36 authentic IDs, GIGASTROKE EAS and EUR AIS original 2022 discovery counts and BBJ EAS discovery inclusion) and independently fetches three public source documents.

| Official public source | Original accessed URL | Actual bytes | SHA-256 source integrity |
|---|---|---:|---|
| Original eQTL Catalogue qtlmap studyFile / sample inclusion specification | https://raw.githubusercontent.com/eQTL-Catalogue/qtlmap/master/docs/usage.md | 12,249 | 7f3143a28ed67df7a2c56dc5403b269cff5cb04390f8cf3fce538c9eb663bbe0 |
| 2021 eQTL Catalogue processed count public Zenodo repository metadata | https://zenodo.org/api/records/4678936 | 11,721 | 818a5ebd7f448dd685f3f23d90d42f719a520a8fc27db233cb624bee21619682 |
| Original 2022 GIGASTROKE paper full text (Mishra et al.) | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9524349/fullTextXML | 505,858 | 7f1105825c4111f73f8b8c36f55897e5c1a19b39b0b240eabc9da3822ac89f79 |

### 1. Original published processed counts actually exist

The original eQTL Catalogue publication Kerimov et al. links [Zenodo 10.5281/zenodo.4678936](https://zenodo.org/records/4678936), and its **public original API response has 32 files**. Four explicitly refer to the two source cohorts:

| Study | Original public filename | Exact deposited bytes | Publisher MD5 |
|---|---|---:|---|
| GTEx v8 | GTExV8_merged_counts_pseudo.tsv.gz | **735,218,485** | ec951a9499c88e4f9071f54f12a8990d |
| GTEx v8 | GTExV8.tsv | **1,769,117** | 0e29ac0664ad10403edf7a833d0a4503 |
| ROSMAP | ROSMAP_merged_counts_pseudo.tsv.gz | **25,484,860** | d215dc801f687e15d5b938356d80d183 |
| ROSMAP | ROSMAP.tsv | **53,064** | 038f2fc144129d401144f3924109e1df |

This is **actual original repository metadata read from Zenodo**; it is not inferred from an article paragraph. The count matrices were NOT downloaded, and the sample identifiers were NOT retrieved or uploaded. The current public original record represents the 2021 compendium, while the candidate QTL sources are **eQTL Catalogue Release7 QTS000015/QTD000131 GTEx aorta, QTD000171 GTEx cortex, and QTS000025/QTD000434 ROSMAP**.

None of the 32 repository **filenames** is explicitly labelled as the QTLtools/qtlmap analysis-stage fully normalized / inverse-normal-transformed matrix; **filename absence alone does not establish that an equivalent analysis matrix is unavailable**. In particular, **raw or processed counts are not necessarily the same matrix** as genotype-compatible normalized expression traits actually tested in Release7. The Zenodo 2021 counts file cannot be used directly as a direct estimate of the Release7 specific normalized trait sdY without complete preprocessing, tissue sample selection and qtlmap source-version/sample metadata crosswalk.

### 2. Exact original QTL regression sample inclusion rule

The official public [eQTL Catalogue qtlmap studyFile specification](https://github.com/eQTL-Catalogue/qtlmap/blob/master/docs/usage.md) explicitly requires:
- count_matrix = **normalized phenotype input matrix** (phenotype_id and sample columns),
- pheno_meta = exact molecular phenotype position/ID metadata,
- sample_meta = **sample_id, genotype_id and qtl_group** metadata,
- genotype VCF matching the genotype ID,
- qtlmap **uses the intersection of count_matrix and sample_meta samples**. Genes may be filtered according to expression/covariate QC.

Thus the earlier observed SNP-specific genotype allele count AN/2 (GTEx aorta 366–387, GTEx cortex 186–205, ROSMAP 560) does not independently verify expression-sample inclusion, covariate degrees of freedom, or the target gene's actual phenotype sdY. **All 36 genes still have measured normalized phenotype sdY = NOT VERIFIED and source regression effective N = NOT VERIFIED.** The earlier SE×sqrt(2f(1−f)N) residual scale is explicitly a covariate-ignorant diagnostic, not measured trait sdY.

### 3. GIGASTROKE 2022 primary paper confirms mixed cohort model origins

Official publisher paper DOI **10.1038/s41586-022-05165-3** was re-verified directly by EuropePMC original fullTextXML. Its methods describe **individual-cohort additive genetic association primarily by logistic regression, with Cox models in some longitudinal population-based cohorts**, followed by a **fixed-effect inverse-variance-weighted (IVW) meta-analysis** of stroke GWAS estimates.

This **does not establish that beta in the final 2022 GWAS Catalog original files** for EAS AIS GCST90104545 and EUR AIS GCST90104540 is uniformly the direct **natural-log odds ratio** needed for a case-control effect-size colocalisation interpretation, especially for mixed logistic and hazard-ratio study inputs. The release's original YAML supplies trait/build/data file publisher MD5 but **no direct study-specific beta unit field**. The paper source is a method-level observation, **not exact file per-SNP scale proof**.

The [official coloc data structure guide](https://chr1swallace.github.io/coloc/articles/a02_data.html) confirms that original GWAS signed beta and varbeta with type=cc need not invoke an arbitrary SNP-level GWAS N_eff/s in the effect-size input route, **if the effect-size units and other assumptions are appropriate**. P-only mode needs N/s, and quantitative QTL beta/varbeta requires a defensible direct sdY or a scientifically appropriate N/MAF route.

**Original GIGASTROKE discovery cohort totals** were previously independently verified from 2022 published Supplementary Table 1: EAS AIS 19,032 cases/237,242 controls; EUR AIS 62,100 cases/1,234,808 controls. BBJ contributed 17,651 EAS AIS cases and 159,610 controls and cannot be used as independent EAS replication.

### 4. Decision matrix, all 36

For each of the 36 original QTL gene×tissue combinations (20 genes, one chr10 locus), the verified public source matrix records:
- exact QTL study ID and original tissue;
- existing original nonpal QTL SNP×gene association record count and per-SNP AN/2 heterogeneity;
- genotype AN/2 metadata verified from prior source work;
- **Release7 actual qtlmap-normalized phenotype matrix ID and genotype-linked sample intersection NOT verified**;
- **direct, source-specific gene expression sdY NOT measured**;
- **exact GWAS Catalog EAS/EUR beta log-OR unit NOT verified**;
- prior BBJ overlap restriction remains explicit;
- **coloc.abf eligibility/execution 0/36, coloc.susie eligibility/execution 0/36; no H3/H4, no causal gene claim, no invented N_eff or sdY**.

## Next evidence priorities

1. **Inspect the two small public GTExV8.tsv and ROSMAP.tsv metadata files**, identify available study/tissue/sample columns and donor aggregation, and determine what can be reconciled to the 2021 count files, then separately establish whether an equivalent **Release7 normalized qcnorm input and qtlmap sample_meta** are available and source-identical. Do not conflate v8 original eQTL files with Catalogue Release7 transforms or infer normalized sdY from raw counts.
2. **Trace the exact GIGASTROKE AIS publisher release effect-scale provenance** (original per-cohort/whole-meta beta definitions; OR versus HR, or explicit source descriptor) and test all required effect-size unit transformations before admitting any CC beta/varbeta Coloc analysis.
3. Once exact quantitative phenotype sdY/sample model or well-founded alternate estimator AND appropriate GWAS effect units are accepted, reassess all 36 original groups jointly. R-based coloc.abf and SuSiE require separate methodological review and prior sensitivity; do not fabricate posterior probabilities.

Source script: scripts/is/audit_p3a_36_source_expression_and_gwas_beta_units.py. Tests: tests/test_is_p3a_36_source_expression_and_gwas_beta_units.py. Workflow: .github/workflows/is-p3a-36-source-expression-gwas-beta-units.yml. Scientific issue [#307](https://github.com/kimtk94/IS_Analysis_V3/issues/307) remains OPEN.
