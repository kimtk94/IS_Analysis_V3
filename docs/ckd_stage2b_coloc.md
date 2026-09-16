# CKD Stage 2B: locus matching and colocalization

## Why this stage is build-aware

The downloaded UKB-PPP per-protein summaries use GRCh38 coordinates, whereas the
Stanzick 2021 CKDGen eGFR summary is GRCh37/hg19. Stage 2B does not directly join
their raw positions.

For each candidate:

1. use the reviewed GRCh38 gene coordinate file;
2. stream only the candidate's chromosome member from its UKB-PPP tar archive;
3. extract the gene interval ±1 Mb;
4. retain additive biallelic SNVs with INFO >= 0.3;
5. lift GRCh38 SNP positions to GRCh37 with the UCSC hg38ToHg19 chain;
6. reverse-complement alleles when the chain maps to the reverse strand;
7. scan the Stanzick eGFR file once and match by GRCh37 chr:pos + allele pair;
8. align outcome beta to UKB-PPP ALLELE1;
9. run coloc.abf.

BOLT/REGENIE-style UKB-PPP fields are interpreted as:
- ALLELE1: effect allele
- A1FREQ: ALLELE1 frequency
- BETA: ALLELE1 effect

## Initial coloc model

The first pass uses coloc.abf with:
- exposure type = quantitative
- exposure sdY = 1 because UKB-PPP protein phenotypes were inverse-rank normalized
- outcome type = quantitative
- outcome sdY estimated by coloc from beta variance, MAF, and N
- p1 = p2 = 1e-4
- default p12 = 1e-5
- p12 sensitivity = 1e-6, 1e-5, 1e-4

This is a single-causal-variant screen. Complex/multi-signal loci, especially
HLA-E/MHC, should not be finalized from coloc.abf alone. LD-aware conditional or
SuSiE colocalization is a later sensitivity step.

## Run

```bash
cd /srv/is-analysis/IS_Analysis_V3
git pull --ff-only origin main

CKD_STAGE2_ROOT=/srv/is-analysis/results/ckd/stage2 \
CKD_STAGE2_PQTL_ROOT=/srv/is-analysis/data/ckd/stage2_pqtl \
CKD_STAGE2B_ROOT=/srv/is-analysis/results/ckd/stage2b_coloc \
SYNC_DRIVE=1 \
RCLONE_REMOTE=gdrive \
  server/ckd_run_stage2b.sh
```

Main results:

```text
/srv/is-analysis/results/ckd/stage2b_coloc/
├── pqtl_locus/*.tsv.gz
├── coloc_input/*.tsv.gz
├── STAGE2B_MATCH_QC.tsv
├── STAGE2B_PREP_SUMMARY.json
├── coloc_results/
│   ├── STAGE2B_COLOC_DEFAULT.tsv
│   ├── STAGE2B_COLOC_SUMMARY.tsv
│   └── <GENE>_coloc_snps.tsv
└── SHA256SUMS.txt
```
