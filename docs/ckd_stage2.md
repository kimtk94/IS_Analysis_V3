# CKD Stage 2: candidate planning and full pQTL acquisition

Stage 2 narrows the Stage 1 EUR eGFR discovery screen before locus-level colocalization.

## Follow-up gate

A protein proceeds when all are true:

- EUR eGFR BH-FDR < 0.05
- at least 3 supporting renal phenotypes show the same kidney-risk direction
- at least 2 supporting phenotypes are nominally significant in that same direction
- the pre-specified strongest-F anchor SNP is available in EAS eGFR
- EAS eGFR direction is concordant

This is a follow-up gate, not a causal classification.

The current Stage 1 results yield 9 assays:
UMOD, INHBC, SDCCAG8, ACP1, GSTA3, GSTA1, HLA-E, F12, CPVL.

The exact Olink assay is matched by gene + UniProt accession + OID + assay version.
This prevents accidentally downloading another assay for the same gene.

## Server plan only

```bash
cd /srv/is-analysis/IS_Analysis_V3
git pull --ff-only origin main

CKD_STAGE1_ROOT=/srv/is-analysis/results/ckd/stage1 \
CKD_STAGE2_ROOT=/srv/is-analysis/results/ckd/stage2 \
SYNC_DRIVE=1 \
  server/ckd_run_stage2.sh
```

## Download EUR full protein summary archives

Synapse authentication is required once:

```bash
cd /srv/is-analysis/IS_Analysis_V3
export PATH="$PWD/.venv-ckd/bin:$PATH"
python3 -m venv .venv-ckd  # only if the CKD venv does not already exist
.venv-ckd/bin/python -m pip install synapseclient
synapse login
```

Then:

```bash
CKD_STAGE1_ROOT=/srv/is-analysis/results/ckd/stage1 \
CKD_STAGE2_ROOT=/srv/is-analysis/results/ckd/stage2 \
CKD_STAGE2_PQTL_ROOT=/srv/is-analysis/data/ckd/stage2_pqtl \
STAGE2_DOWNLOAD=1 \
SYNC_DRIVE=1 \
RCLONE_REMOTE=gdrive \
  server/ckd_run_stage2.sh
```

Only the 9 EUR archives are downloaded at this stage. The 9 EAS archive rows are
recorded but marked optional until EUR locus-level colocalization is evaluated.

## Outputs

```text
results/ckd/stage2/
├── stage2_candidates.tsv
├── stage2_anchor_regions_hg19.tsv
├── stage2_pqtl_download_manifest.tsv
├── stage2_pqtl_download_progress.tsv     # after download
├── stage2_pqtl_schema.tsv                # after download
├── STAGE2_PLAN_SUMMARY.json
└── SHA256SUMS.txt
```

Raw candidate archives:

```text
/srv/is-analysis/data/ckd/stage2_pqtl/EUR/<GENE>/<archive>.tar
```

After archive headers are inspected, the next step is to standardize the full
protein summary fields, intersect the cis locus with the full EUR eGFR GWAS by
rsID/allele, and run coloc using the matched locus variants.
