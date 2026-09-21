# CKD Stage 3A: kidney-specific QTL evidence

Stage 3A attaches kidney-tissue genetic evidence to the nine Stage 2 candidates.
It is deliberately evidence-preserving rather than score-based.

## Public data sources

1. Hirohama et al., Nature Medicine 2025, human kidney proteogenomics
   - kidney proteomics/pQTL: n=337
   - same-study kidney eQTL: n=315
   - Supplementary Tables 1-30 are retrieved as the published PMC workbook
     `NIHMS2102371-supplement-Supplementary_Tables.xlsx`.
   - The downloader validates that the response is a real XLSX container before parsing; short HTML/XML error responses are rejected.
2. Susztak kidney eQTL meta-analysis: n=686, significant SNP-gene pairs at q<0.01.
3. Tubule eQTL: n=356, significant SNP-gene pairs at FDR<0.05.
4. Glomerulus eQTL: n=303, significant SNP-gene pairs at FDR<0.05.

The eQTL bulk downloads contain significant pairs only. A zero candidate hit must
not be interpreted as proof that the gene is not expressed or has no cis genetic
regulation.

## Data-use boundary

The Susztak Kidney Biobank publishes a user agreement. The runner does not accept
that agreement on the user's behalf. It requires `ACCEPT_SUSZTAK_TERMS=1` after
the user has reviewed the current agreement at:

https://susztaklab.com/agree.php

Raw kidney QTL files are stored only under
`/srv/is-analysis/data/ckd/stage3a_kidney` and are not copied to GitHub or Google
Drive by this runner. Only candidate-level derived evidence tables are eligible for
Drive sync.

## Run

```bash
cd /srv/is-analysis/IS_Analysis_V3
git pull --ff-only origin main

ACCEPT_SUSZTAK_TERMS=1 \
CKD_STAGE2C_ROOT=/srv/is-analysis/results/ckd/stage2c_susie \
CKD_STAGE3A_RAW_ROOT=/srv/is-analysis/data/ckd/stage3a_kidney \
CKD_STAGE3A_ROOT=/srv/is-analysis/results/ckd/stage3a_kidney \
SYNC_DRIVE=1 \
RCLONE_REMOTE=gdrive \
  server/ckd_run_stage3a.sh
```

## Outputs

```text
/srv/is-analysis/results/ckd/stage3a_kidney/
├── STAGE3A_KIDNEY_EVIDENCE.tsv
├── STAGE3A_HIROHAMA_SUPPLEMENT_HITS.tsv
├── STAGE3A_EQTL_META686_HITS.tsv
├── STAGE3A_EQTL_TUBULE356_HITS.tsv
├── STAGE3A_EQTL_GLOMERULUS303_HITS.tsv
├── STAGE3A_PROVENANCE.json
└── SHA256SUMS.txt
```

`STAGE3A_KIDNEY_EVIDENCE.tsv` is the primary candidate-level table. It combines
Stage 2 ABF/SuSiE status with kidney pQTL, same-study kidney eQTL, published eGFR
colocalization supplement hits, and independent kidney meta/tubule/glomerulus eQTL
hit counts.

## Interpretation order

The primary Stage 2 shared-signal candidates (`SDCCAG8`, `GSTA3`, `ACP1`) are
carried forward first for tissue replication. `INHBC` and `GSTA1` remain explicit
discordance/sensitivity candidates rather than being silently discarded. `UMOD`
remains a special multi-signal/LD diagnostic locus.


## Supplementary workbook fallback

Publisher and PMC supplementary-file endpoints may return anti-bot HTML or a
small metadata response to scripted `curl` requests. The Stage 3A runner
therefore treats the Supplementary Tables workbook as optional:

- it tries PMC and bounded Nature/Springer media-object URLs;
- every candidate file must validate as a real XLSX with >=20 worksheets;
- if all scripted downloads fail, Stage 3A continues instead of aborting;
- candidate-level values directly printed in Hirohama et al. Table 1/main text
  are preserved for ACP1, GSTA1, and INHBC;
- independent kidney eQTL meta/tubule/glomerulus downloads still run normally;
- provenance records whether the workbook was available.

The fallback must not be interpreted as a negative result for candidates absent
from the hard-coded main-text rows; it means the supplementary workbook could
not be programmatically retrieved in that run.
