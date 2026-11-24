# Ocular P5c independent TM RNA-seq replication metadata gate — 2026-10-31

## Real execution evidence
[GitHub Actions 41086796126](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41086796126) [aggregate JSON artifact 12992741024](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41086796126/artifacts/12992741024).

Source is the ENA public read_run filereport on PRJNA1243143, **not** a downloaded FASTQ. It exposed eleven available metadata fields including sample_alias, sample_title, sample_description, experimental_factor, experiment_title, library_layout and FASTQ file endpoints. Actual run: **46 read runs, 46 unique BioSamples, 46 FASTQ endpoints**.

The sample-specific fields contained **no explicit DEX versus vehicle labels, no classifiable 1h/6h/2d timepoints, and no TM cell-strain mapping**. All 46 BioSamples classified `UNKNOWN` for treatment and timepoint under a conservative, predeclared parser. Study-wide/experiment titles may mention the word dexamethasone in every record, but **must not be used** as per-run treatment labels.

### Verdict
**HOLD_EXTERNAL_SAMPLE_DESIGN**. Actual existence of 46 read runs alone cannot establish 46 independent donors, treatment pairs, timepoint matching or true independent DEG replication. Do not claim external validation is completed. We need the original paper's sample-to-run submission key, authors' design map or repository sample attributes updated with treatment/time/biological strain, then re-check before any FASTQ download.

### Existing paper design and alternative
Published IOVS 2026 Mehrotra et al. work reports 9 QC-passing human TM strains and 4-5 SCE strains across varying times, not 46 independent biological donors. Literature-only direction crosswalk exists for FKBP5 and FAM107A and counts as biological face validity only. Distinguish 100 nM/48h in the independent study vs nominal GSE313736 500 nM/48h before comparing gene-level effect magnitude.

### Research status
- GSE313736 in-vitro BeamA primary count-based edgeR: PASS, 57 biological pairs, 16,050 tested, 13,717 BH-FDR<0.05.
- Published gene crosswalk: see `P5B_LITERATURE_GENE_CROSSWALK_20261031.md`.
- P5b formal 9-pool leave-one-out and MDS: workflow `.github/workflows/ocular-p5b-pool-sensitivity.yml`.
- Clinical steroid-induced ΔIOP response prediction: **STOP for current dataset**, no matched clinical outcomes. Distinct controlled-access dbGaP phs004037 required for any clinical genomics work.
