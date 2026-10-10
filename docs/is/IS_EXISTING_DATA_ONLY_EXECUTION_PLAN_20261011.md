# IS existing-data-only research execution plan

Date: 2026-10-11 (KST). Scope: **No external database request, purchase, restricted download, or controlled-access data dependency.** This is a workstream prioritization and reproducible evidence inventory, not a new causal analysis.

## Direct server audit

The evidence inventory was recomputed from existing local files using `scripts/is/audit_existing_data_no_download_backlog.py`, and passed source consistency checks. Its outputs are under:

`/srv/is-analysis/results/is/audits/is_no_new_db_backlog_20261011_v1/`

- `IS_NO_NEW_DB_BACKLOG_AUDIT.json`: row counts, statuses, SHA-256 of nine input files, scientific caveats, task ledger.
- `IS_NO_NEW_DB_WORKSTREAMS.tsv`: ordered seven workstreams.

Verified baseline:
- **80 provisional regions**: 30 EAS/Japanese components + 50 EUR positional clusters. Not 80 independent causal loci.
- **2,225 unique positional gene IDs**, **2,425 gene-region associations**, **2,426 evidence table rows** (includes one legacy-anchor-only row).
- EAS/Japanese 869 region-gene rows, EUR 1,556 rows; QTL labels: **2,383 `NOT_TESTED`** eQTL and 43 `REASSESS_LEGACY`; all 2,426 rows have `sqtl_status=NOT_TESTED` and `pqtl_status=NOT_TESTED`.
- **646 archived GTEx-v8 coloc ABF gene-tissue assays**: latest local per-SNP replay ledger has **30 PASS**, **616 MISSING_INPUT**. The latter is missing *replay source files in this cache*, not a negative association or proven unavailability from original repositories. Don't call 646 assays independently recomputed.
- **Eight selected gene-tissue direct per-SNP replay entries** verified: FGF5, CALHM2, NEURL1, C4orf22, INA, SH3PXD2A, COL4A2, COL4A1. The existing `p12` grid comprises 45 records for selected genes; changed prior and QTL sample-size sensitivity are not independent mechanism validation.
- Adult **healthy** human vascular/perivascular temporal-lobe scRNA reference: 80,515 cells from 6 donors, 585 pseudobulk rows, 10/10 donor-paired descriptive QC PASS. There is **no ischemic-stroke patient-versus-control differential-expression contrast**.
- Externally matched Japanese cohort LD and molecular-QTL matched LD remain unresolved; no previously BLOCKED candidate becomes causal due to presentation or ranking.

## Execute in order

### P0-A. Evidence-universe matrix: retain all 2,225 genes

Make one traceable row per stable gene ID, with nested/long-form links to its 2,425 region-gene associations. Include source region, ancestry, effect direction where valid, window distance, whether GWS or suggestive, prior coloc/QTL `TESTED/NOT_TESTED`, human healthy-reference presence and literature sources. Deduplicate by stable Ensembl ID (not merely symbol) and preserve legacy-anchor-only details without inflating the positional count. A heuristic **triage rank** may be produced but **never a causal rank**, and must not trim the candidate universe. Unknown is not negative.

Deliverables: `IS_ALL_GENES_EVIDENCE_MATRIX.tsv`, `IS_GENES_BY_WORKSTREAM_READINESS.tsv`, `IS_GENE_EVIDENCE_PROVENANCE.json`.

### P0-B. Existing SNP-level coloc source-input recovery

For the 616 `MISSING_INPUT` rows, inventory **existing** server / Git / permitted Drive backup paths for per-SNP GTEx QTL and BBJ GWAS harmonized inputs. Do not fetch a new DB. Categorize `LOCAL_FOUND`, `ARCHIVE_FOUND`, `NOT_FOUND`, `SCHEMA_CONFLICT` with checksums. Re-run SNP-level `coloc.abf` only when source summary statistics, allele alignment and explicit sample-N conventions are reproducible. Preserve original 646-row master; write derivative outputs and per-pair status. Run `p12` sensitivity on the valid subset, not via post hoc strong-hit cherry-picking.

Deliverables: `IS_646_COLOC_INPUT_RECOVERY.tsv`, `IS_646_SNP_REPLAY_V2.tsv`, `IS_COLOC_P12_SENSITIVITY_V2.tsv`, and source hashes.

### P1-A. Human donor-level functional reference QC

Use 6-donor and 80,515-cell existing healthy reference to evaluate per-donor cell type coverage, library-size imbalance, target-feature detection, and conditionally normalized CPM contrasts in pre-specified cell-type pairs. For `FGF5`, `C4orf22`, document absent features as **NOT ASSESSABLE**, not absent protein function. No stroke differential expression or genotype→expression causation is inferred.

Deliverables: per-donor QC plots, locus × cell-type expression/detection heatmap, donor-paired uncertainty/sensitivity summary.

### P1-B. Existing GWAS EAS/EUR direction and evidence gating

Read existing harmonized BBJ/GIGASTROKE 47,710-row dataset, check allele sign flips and ancestry-specific frequency/heterogeneity. Separate GWS from suggestive regions; don't equate BBJ and GIGASTROKE EAS with independent replication because GIGASTROKE contains BBJ. EUR heterogeneity is cross-ancestry context, with LD and phenotype caveats.

### P1-C. Competing mechanism matrices for key genes

Compare FGF5, CALHM2, NEURL1, C4orf22, INA, SH3PXD2A, COL4A1/2 and ALDH2 across coding, regulatory, tissue, cell-specific, alcohol/metabolic mechanisms. Record H3 vs H4, prior sensitivity, source-harmonization confidence, cell-feature presence. ALDH2 rs671 AIS association may be strong while GTEx bulk eQTL colocalization remains weak; distinguish these.

### P2. Publication-ready figures and dashboard

Plot population locus comparison, multi-omics evidence heatmap with explicit grey NOT_TESTED, 646-coloc coverage bar, `p12` sensitivity, donor-level reference expression and sensitivity, evidence-to-inference flowchart. Update site and existing literature benchmarks only from hashed outputs. A completed plot does not upgrade a scientific gate.

## Explicitly paused

- CKB encrypted IS ZIP: waits for permitted institutional decryption key.
- TPMI 433.21 full GWAS: currently inaccessible by permitted normal route.
- Matched Japanese BBJ genotypes/LD, controlled data or purchase requests.
- Causal MR/instrument independence or disease-cell-specific mechanism claims dependent on those data.

## Guardrails

1. **2225-gene scope stays intact** regardless of data quality or score.
2. `MISSING_INPUT`, `NOT_TESTED`, `FEATURE_NOT_IN_MATRIX`, `NO_ASSOCIATION` and `TESTED_NO_SUPPORT` are different statuses.
3. No claiming shared causal gene from a high `PP.H4` without per-SNP input validation, locus conditioning and appropriate LD.
4. **ALDH2 alcohol fine-mapping BLOCKED; ADH1B alcohol fine-mapping EXPLORATORY** pending existing outstanding science gates.
5. No canonicals overwritten or server git checkout altered; keep IS research work isolated from the CKD server branch.
