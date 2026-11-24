# IS 2,225-gene public nominal-QTL evidence coverage — source audit

**Audit date:** 2026-10-15  
**Audited input:** [IS_2225_GENES_PUBLIC_NOMINAL_QTL_READINESS.tsv](https://github.com/kimtk94/IS_Analysis_V3/blob/research/is-broad-discovery-20261009/docs/data/is/IS_2225_GENES_PUBLIC_NOMINAL_QTL_READINESS.tsv)  
**Git blob SHA:** `a99539a73431dc208d010f1a09256e9aad414624`

> Direct read of the **entire 2,225-row GitHub TSV**. This is a focused **public nominal QTL source coverage table**, not the full 2026-10-11 IS server evidence matrix, not a GWAS/molecular inference reanalysis and not a fresh causal validation. Its all-FALSE molecular columns **cannot erase or contradict** separately archived ABF calculations, which must be reconciled by dataset lineage.

## 1. Scope and denomination

- 2,225 TSV data rows = **2,225 unique Ensembl IDs**.
- 869 rows have an EAS regional identifier; 1,356 rows have no EAS region in **this file**. These are not 1,356 biological null genes.
- 14 columns present, including `GWAS_and_QTL_alleles_coharmonized`, `matched_QTL_LD_verified`, `coloc_recomputed`, `ischemic_stroke_causal_gene`.

| Public source status | Genes | Proper interpretation |
|---|---:|---|
| RAW_NOMINAL_QTL_RECORDS_AVAILABLE_NO_COLOC | **149** | Nominal molecular test rows exist; molecular allele/LD/coloc validation is not attested in this table |
| PUBLIC_NOMINAL_QTL_REGIONAL_PILOT_NOT_EXECUTED | **629** | The specific targeted pilot has not run; not a negative eQTL result |
| NO_NOMINAL_ROWS_IN_TARGETED_SLICES_NOT_NEGATIVE | **91** | No nominal rows in the selected slices; tissue/sample/selection limitations are not a negative biological inference |
| NOT_IN_EAS_30_GROUPS | **1,356** | Outside this public EAS 30-group targeted assay, not out of the study |
| **Total** | **2,225** | All original eligible gene IDs preserved |

The first three categories total **869** EAS group-associated genes.

Within the 869 EAS-associated rows, the GWAS-priority label counts are:
- 207 `T1_BBJ_AND_EAS_AIS_REGION_ASSOCIATED_NOT_INDEPENDENT_REPLICATION`
- 16 `T2_CES_SUBTYPE_ASSOCIATED_EAS_AIS_NOT_GWS`
- 57 `T3_AS_ASSOCIATED_EAS_AIS_NOT_GWS`
- 589 `T4_SUGGESTIVE_REGION_NO_5E8_SIGNAL_IN_PRIMARY_OUTCOMES`

These are source-file **annotation categories**, not 207 independently replicated GWAS loci, not 207 causal genes, and not a valid independent-variant LD count.

## 2. Critical molecular readiness gate in this source

| Source table field | TRUE rows | FALSE rows |
|---|---:|---:|
| GWAS_and_QTL_alleles_coharmonized | 0 | **2,225** |
| matched_QTL_LD_verified | 0 | **2,225** |
| coloc_recomputed | 0 | **2,225** |

All 2,225 rows have `ischemic_stroke_causal_gene = NOT_ESTABLISHED`.

**Interpretation:** This narrowly scoped **public nominal-QTL import** does not yet contain the GWAS–molecular harmonization and cohort-appropriate LD necessary to promote *its own* source records to causal/colocalized status. It does **not** say 0 valid coloc exist in the broader IS research project. A separate historical audit documented 646 ABF assay attempts, 30 numeric PASS attempts, and 14 unique genes represented by 30 replayed assays; source and methods must be compared before manuscript statements.

## 3. Previously highlighted loci — per-row check

| Gene (source symbol) | Stable Ensembl ID | Source region | Nominal source status | Interpretation |
|---|---|---|---|---|
| FGF5 | ENSG00000138675 | G0007 | Raw nominal QTL records | Existing input, this source coloc not validated |
| ALDH2 | ENSG00000111275 | G0022 | Pilot not executed | Keep coding/protein/metabolism branch; expression input absence not rejection |
| SH3PXD2A | ENSG00000107957 | G0015 | Raw nominal QTL records | Cell-specific mechanism remains exploratory |
| COL4A2 | ENSG00000134871 | G0023 | Raw nominal QTL records | Vascular mural/cell evidence needed |
| COL4A1 | ENSG00000187498 | G0023 | Raw nominal QTL records | Competing/secondary candidate at G0023 |
| CALHM2 | ENSG00000138172 | G0015 | Raw nominal QTL records | Separate source status, not a causal verdict |
| INA | ENSG00000148798 | G0015 | Raw nominal QTL records | Annotation only |
| C4orf22 | ENSG00000197826 | G0007 | No nominal rows in targeted slices | Not an evidence-of-absence finding |
| NEURL (legacy NEURL1 alias) | ENSG00000107954 | Gene alias requires lineage check | Not assumed from symbol NEURL1 | **Use stable Ensembl ID in all joins** |

Source `region_id` values are prefixed `IS_XDATA_`. Table abbreviations above show the G-code suffix only.

## 4. Next action and scope boundary

1. **Freeze this 2,225-row source coverage denominator**. Avoid artificial narrowing to 4 core genes. The 149 rows with raw nominal records form a source-availability queue, not an unconditional request to re-run 149 coloc tests.
2. Reconcile these **public nominal sources** to the archived 14-gene/30-ABF multi-signal readiness and the wider 2,425 gene–region project manifest. Missing cohort-matched QTL LD must remain UNRESOLVED.
3. Use existing allele-harmonized BBJ/GIGASTROKE/EBI plotted loci when legitimately supported; label marginal GWAS/eQTL differently from colocalization and mediation.
4. Produce a candidate evidence matrix with **SOURCE_AVAILABLE / UNVERIFIED / UNRESOLVED / NOT_ASSESSED** labels. Numeric thresholds do not replace ancestry/cohort matching.
5. **Stop** further data-source hunts or algorithm extensions unless a missing input is essential to a specific high-impact workflow manuscript claim.

**Acceptance test:** every reported candidate and manuscript figure can be traced to a source path, version, ancestry, allele/build, exact statistical test, and limitations. Lack of such evidence means 'not established', not 'effect absent'.
