# IS Phase10/11 — 646-assay direct SNP replay and 1000G EAS LD audit

**2026-10-11 update | FULL_646_NUMERIC_REPLAY VERIFIED (646 PASS / 0 FAIL / 0 MISSING); remaining scientific causality/LD gates BLOCKED.** [Independent full source audit](IS_PHASE10_11_646_FULL_REPLAY_VERIFIED_20261011.md).

## Scope and source files

The historical GTEx ABF screen contains 646 gene–tissue tests from **four original BBJ ischemic-stroke loci**, not the expanded 80 provisional genomic regions.

Original Drive paths below are relative to MyDrive/MASTER_DEGREE/IS_COLAB:

- Input index: coloc_ready_v2/COLOC_ABF_INPUT_INDEX_V2.tsv
- Original ABF master: results_v2/COLOC_ABF_MASTER_V2.tsv
- Actual per-SNP inputs: coloc_ready_v2/abf/
- Archived settings: coloc.abf p1=p2=1e-4, p12=1e-5, GWAS N=174686, cases=22664, QTL first-scalar N. The archived notebook code supports these literal parameters, but does not attest the exact original runtime.

## Full original 646-row replay implementation

- Batch R script: scripts/is/replay_is_legacy_646_batch.R
- One-Python-cell notebook: notebooks/is/IS_PHASE10_11_646_SNP_COLOC_REPLAY_ONECELL.ipynb
- Notebook builder: scripts/is/build_is_646_snp_replay_colab.py
- Pinned coloc 5.2.3.
- Each source row requires correct locus, tissue, gene, unique SNP keys, source SNP denominator, expected QTL N, finite beta/SE/MAF and direct original H0–H4 replication within absolute error 1e-8.
- Input absent -> MISSING_INPUT; mismatched or unreadable -> FAIL with reason; genuine full-H0–H4 SNP replay -> PASS.
- Each job refuses an existing output folder. No original inputs, canonical data, original ABF results or historical figures are overwritten.

**Historical 2026-10-10 server-cache run only (superseded by full Colab run):** 646 status rows, **30 PASS / 0 FAIL / 616 MISSING_INPUT**. This has also passed against the unannotated raw source master, not only its annotated copy. This is a successful **partial-coverage accounting run**, NOT a completed 646-input replay.

## Run full 646 inputs in Google Colab

[Open and run the one-cell Colab notebook](https://colab.research.google.com/github/kimtk94/IS_Analysis_V3/blob/feat/is-p0-evidence-20261010/notebooks/is/IS_PHASE10_11_646_SNP_COLOC_REPLAY_ONECELL.ipynb)

It mounts Drive and reads every source SNP file directly, bypassing the server Drive API rate limit. It checks 646 input filenames, installs R/coloc 5.2.3 if absent, embeds and parses the exact batch R source, and writes to a **new** Drive directory:

MASTER_DEGREE/IS_COLAB/results/IS_PHASE10_11_LEGACY_646_SNP_REPLAY_V1

The output includes full 646-row status, PASS-only table, per-source SHA256, original source index/master SHA256, runtime log and execution manifest. Only an all-646 PASS result can be called FULL_646_NUMERIC_REPLAY. **Verified 2026-10-11:** The one-cell Colab completed with all 646 original source SNP tests PASS; input/source SHA256, six Drive outputs, H0-H4 full matrix and logged exit status are independently checked. See [final full audit](IS_PHASE10_11_646_FULL_REPLAY_VERIFIED_20261011.md).

## 1000G EAS504 reference LD source audit

Three locus-labelled matrices in results/is/stage3_finemap/japan/bbj/ld_v3/ have verified PLINK logs and matching PSAM data with **504 reference participants**. The source is locally labelled **1000G_EAS_504**, *not* BBJ cohort-specific genotype LD and *not* JPT-only LD. BBJ_IS in filenames denotes GWAS locus.

All 15,749 gene–SNP records across 8 selected original source inputs match their 1000G EAS panel's variant IDs by **GRCh37 GWAS variant_id** (100% source SNP coverage). None match by GRCh38 match_key. Variant order maps and matrix shape/bytes are verified; this does not establish orientation or numerical matrix LD quality.

- LD source audit code: scripts/is/audit_is_1000g_eas_ld_priority8.py
- Separate detailed output: results/is/audits/is_priority8_1000g_eas_ld_alignment_20261010_v1/
- Full direct-SNP batch trial output: results/is/audits/legacy_coloc_646_batch_original_master_20261010_v1/

**Critical blocker:** GTEx molecular cohort matching LD is not established in this audit. GWAS ancestry-matched 1000G EAS LD alone is insufficient for valid joint GWAS–QTL SuSiE; multi-signal evidence stays BLOCKED.

## Scientific interpretation and remaining gates

1. Complete original 638 source-level H0–H4 tests without collapsing failed tests into PASS, preserve all file hashes and software provenance.
2. Verify GRCh37↔GRCh38 variant conversion, effect allele orientation, missing variant coverage and tissue/gene selection multiplicity.
3. Acquire and validate a suitable molecular eQTL LD reference; use ancestry/cohort matched multi-signal methods.
4. Report the distinct layers of GWAS association, coloc shared causal variant, healthy cell-type localization, mediation and disease-specific mechanism.
5. Preserve expanded 80 provisional GWAS regions and 2,225 positional genes as the broader candidate universe.

**No new IS causal gene established.**

## Subsequent partial-server replay snapshot

- The new cache21 report independently reproduced 21/646 H0–H4 results (0 FAIL, 625 MISSING_INPUT) with maximum absolute H0–H4 difference < 1.7e-15. No inference promoted.
- The 21 cached gene–tissue inputs consist of 37,771 source gene–SNP rows and 3 1000G EAS504 LD reference panels. Numerical LD tests verify finite, symmetric, diag=1, bounds [-1,+1], 128 sampled eigenvalues nonnegative within numerical tolerance, and GRCh37 REF/ALT ID agreement; this is NOT cohort-matched GTEx LD.
- Details and source-hash provenance: [21-input replay and numeric LD audit](IS_PHASE10_11_646_SOURCE_REPLAY_CACHED21_20261010.md).

## Latest partial-source snapshot: 30 inputs

Direct SNP replays: 30/646 PASS, 0 FAIL, 616 server inputs not cached. Four 1000G EAS504 panels pass numerical correlation QC (53,979 repeated gene–SNP rows). See [cache30 source audit](IS_PHASE10_11_646_SOURCE_REPLAY_CACHED30_20261010.md). Cache21 remains a preserved historical snapshot; a subsequent 2026-10-11 Colab execution independently verified all 646 source tests (see full report).
