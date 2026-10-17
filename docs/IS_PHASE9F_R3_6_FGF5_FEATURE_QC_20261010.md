# IS Phase9F-E R3_6 — FGF5 feature-space verdict (2026-10-10)

## Execution provenance

- Notebook: IS_Phase9F_E_Human_Vascular_R3_6_FGF5_FEATURE_QC_PARSECHECK_ONECELL.ipynb
- [Colab source](https://colab.research.google.com/drive/1fn3B3eNRfYB69RROIPtel7aOW7TW3x5S)
- [R3_6 Google Drive results](https://drive.google.com/drive/folders/1pq7IzN4rCa4oREHuAmIsfXMvU8Synl52)
- GSE256493: adult/control temporal lobe endothelial/perivascular Seurat reference.
- SHA256 of original RDS: 337c3d82c070c7ea71ad21114af847e5f4c76ab199e70c2ebcdce6b03db104dd.
- R parser PASS before loading; Seurat loaded; RNA assay 17,349 genes x 80,515 cells; eleven author cell-type categories; six donor labels TL5, TL6, TL8, TL7, TL9, TL10.
- Execution PHASE9F_E_HUMAN=COMPLETE; original processed outputs preserved.
- Drive FGF5 diagnostic: https://drive.google.com/file/d/1EwbtzvfwYRsDs78ZPRVeCkaXxEeKiGkk/view
- Drive run manifest: https://drive.google.com/file/d/1Xgvj5dvG08X9PESDpAUsJ12pmq9EkV8_/view

## Data-supported FGF5 result

| Criterion | R3_6 |
|---|---|
| RNA assay features | 17,349 |
| Exact FGF5 symbol RNA feature matches | 0 |
| ENSG00000138675 RNA feature matches | 0 |
| Alternate feature-metadata match count | 0 |
| Available gene/symbol/feature metadata fields | none |
| Raw result | UNRESOLVED_IN_DISTRIBUTED_SEURAT_RNA_ASSAY |
| Scientific status of cell localization | NOT_ASSESSABLE_IN_THIS_REFERENCE |

This verifies that the matched FGF5 feature is absent in the **released,
processed Seurat RNA assay**. Zero feature matches should **not** be
reported as zero-expressing cells, absence of FGF5 biology, absence of the
GWAS signal, or failed colocalization. The raw FASTQ/gene-filtering stage
has not been audited and is not inferable from this Seurat object. The
metadata-field count was zero, so no additional Ensembl-to-symbol metadata
alias resolution could be performed.

## Independent expression context (external; not GSE256493)

- Human Protein Atlas [FGF5 cerebellum](https://www.proteinatlas.org/ENSG00000138675-FGF5/brain/cerebellum):
  GTEx Brain-Cerebellar Hemisphere mean **2.8 nTPM**, **277 samples**;
  HPA Brain cerebellum reported **24.4 nTPM**. Distinct source assays:
  do not combine these values.
- Human Protein Atlas [FGF5 blood vessel](https://www.proteinatlas.org/ENSG00000138675-FGF5/tissue/blood%2Bvessel):
  consensus and GTEx blood-vessel RNA summary **0.0 nTPM**.
- This resolves the 'FGF5 absent from all brain' concern, not the
  gene-specific stroke causal mechanism or specific brain endothelial expression.
- Prior GTEx eQTL sample size used in locus coloc (~175) is **not** the
  tissue expression sample count (277) and must not be conflated.

## Other human reference branches

- SH3PXD2A: top detection 37.74% in oligodendrocytes; large microglia
  and macrophage group exists (33,630 cells) but with relatively low
  expression. A cell type ranking is not cell-specific regulatory causality.
- COL4A2: top detection 64.69% in smooth muscle cells; structural
  basement-membrane genes may be broadly expressed.
- Other detected genes: ALDH2, COL4A1, CALHM2, NEURL1, INA.
- C4orf22 likewise not represented as a matched RNA feature.

## Next precommitted gates

1. Freeze R3_6 source/result hashes and mark FGF5 cell-localization
   NOT_ASSESSABLE_IN_THIS_REFERENCE in evidence dashboard.
2. Prioritize GTEx/HPA cerebellum and vascular data **without** treating
   FGF5 GTEx bulk expression as stroke-specific evidence.
3. Re-run donor-level (six patients) rather than cell-level-as-independent
   pseudobulk specificity checks for SH3PXD2A and COL4A1/A2; donor representation
   per cell type, library sizes and marker QC must be reported.
4. FGF5 locus: p12 colocalization sensitivity, harmonization/LD source,
   pQTL availability and MR/BP pleiotropy evaluation. Still exploratory
   pending paper-grade independent validation.

**No raw data or canonical analysis tables were modified by this documentation change.**
