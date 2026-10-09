# IS broad discovery: GWAS allele-overlap audit / fine-mapping readiness

Date: 2026-10-09. Working branch: `research/is-broad-discovery-20261009`.

## Why the expanded analysis needs another gate

Existing work has 30 provisional coordinate-overlap components (seven containing GWS signals), and 869 positional gene-region associations. New EAS 1KG GRCh37 genotype reference was acquired for G0004, G0008, G0021; four anchor references already exist (G0007, G0015, G0022, G0023).

**Genome-wide significant coordinate-overlap groups are NOT yet statistically independent GWAS signals.** The reference LD proxy inventory cannot substitute for a complete study-specific signed LD matrix or multi-signal SuSiE fine-mapping.

## New allele-overlap stage

Command:
```bash
python3 scripts/is/audit_gws_gwas_reference_overlap.py
python3 scripts/is/build_gws_study_readiness.py
python3 -m unittest discover -s tests -p 'test_is_*' -q
```

Inputs:
- Five GIGASTROKE EAS GRCh37 canonical GWAS datasets (AS/AIS/CES/LAS/SVS), one BBJ ischemic stroke dataset.
- Seven local 504-sample EAS PVAR indexes with stable `chr:pos:REF:ALT` identifiers.

The audit scans *all* eligible SNPs in seven GWS windows, checks per-position variant pairs against local reference, classifies effect-allele orientation, and separates palindromic variants requiring human review. It never guesses strand flips. For BBJ sources that contain REF/ALT, mismatched reference orientation is explicitly flagged. No source is modified.

Artifacts:
- `IS_GWS_GWAS_REFERENCE_OVERLAP.tsv` — study/group coverage and unresolved allele types.
- `IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz` — reference-matched records with optional ALT-relative beta/EAF, but NOT clinical or causal inferences.
- `IS_GWS_GWAS_REFERENCE_OVERLAP_SUMMARY.json`.
- `IS_GWS_STUDY_INPUT_READINESS.tsv` — source/group/lead status; always blocks full fine-mapping until true LD and signal QC.
- `IS_GWS_STUDY_INPUT_READINESS_SUMMARY.json`.

Scientific limitations: source GWAS may contain study overlap, missing or inconsistent INFO, and no joint LD covariance matrix. Apparent REF/ALT concordance does not prove allele frequency matching, independent replication, sufficient QTL power, or causal proteins/genes.

## EUR acquisition

The source-verified EUR AIS dataset GCST90104540 (`GRCh37`, full original ~283 MB) is retrieved into an isolated noncanonical directory. The download must finish and match the official YAML MD5 before running:

```bash
python3 scripts/is/normalize_gigastroke_verified.py --accession GCST90104540
```

Normalization preserves *effect_allele*, *other_allele*, effect beta/se and reported EAF. The original source lacks reliable REF/ALT; it therefore only creates an unordered allele-pair variant key until genome-reference orientation is independently validated.

## Next calculations

1. Complete EUR raw-file MD5 and normalization QC; audit variants shared with EAS without treating overlapping cohorts as independent replication.
2. Compute ancestry-specific LD matrices for region-wise GWAS summary overlap, with allele/EAF correspondence and numerical PSD checks.
3. Run multi-signal fine mapping, credible-set size/PIP calibration and sensitivity.
4. Apply molecular-QTL coloc only after aligned GWAS LD and complete source provenance. Retain all 869 candidate gene-region records and 9 legacy anchors, with missing evidence recorded as NOT_TESTED.
