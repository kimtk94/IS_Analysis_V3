# CKD Stage 2C: multi-signal SuSiE colocalization

Stage 2B uses `coloc.abf`, which assumes at most one causal variant per trait in a
region. Stage 2C re-runs all 9 candidate loci with `coloc.susie`.

## LD design

- Summary-statistic build: GRCh37/hg19 after Stage 2B harmonization.
- LD reference: 1000 Genomes Project Phase 3 European-ancestry samples.
- Reference genotypes: chromosome-split GRCh37 PLINK2 resources.
- Close relatives: PLINK resource `deg1_phase3.king.cutoff.out.id` is removed.
- Default region: the published Berisa-Pickrell EUR hg19 approximately-independent
  LD block containing the pre-specified Stage 1 anchor SNP.
- If a block cannot be resolved, fallback is anchor +/-500 kb.
- The selected block is intersected with the dense Stage 2B locus; truncation is
  explicitly recorded.
- Variants must be common in the reference (MAF >= 0.01), biallelic A/C/G/T SNPs,
  and match Stage 2B by GRCh37 position + allele pair.

PLINK2 computes signed unphased correlation with `ref-based`. The matrix is then
sign-flipped variant-by-variant so LD is oriented to the same pQTL effect allele
used by both harmonized beta vectors.

The same 1KG EUR reference LD is used for plasma pQTL and EUR eGFR. This is a
reference-LD sensitivity analysis, not study-specific LD.

## SuSiE / coloc

For each gene:

1. build coloc-style pQTL and eGFR datasets with the same SNP order and signed LD;
2. run `runsusie()` once per trait;
3. compare every pair of detected signals with `coloc.susie()`;
4. repeat colocalization for p12 = 1e-6, 1e-5, and 1e-4;
5. compare the default p12=1e-5 result with the Stage 2B single-variant ABF result.

A signal-pair PP.H4 >= 0.8 is flagged for follow-up; PP.H4 >= 0.9 is separately
flagged as stronger evidence. These are analysis flags, not causal verdicts.

## Server run

Install PLINK2 once if needed:

```bash
sudo apt update
sudo apt install -y plink2
```

Then:

```bash
cd /srv/is-analysis/IS_Analysis_V3
git pull --ff-only origin main

CKD_STAGE2_ROOT=/srv/is-analysis/results/ckd/stage2 \
CKD_STAGE2B_ROOT=/srv/is-analysis/results/ckd/stage2b_coloc \
CKD_STAGE2C_ROOT=/srv/is-analysis/results/ckd/stage2c_susie \
CKD_STAGE2C_WORK_ROOT=/srv/is-analysis/data/ckd/stage2c_ld \
STAGE2C_THREADS=2 \
STAGE2C_MEMORY_MB=2500 \
SYNC_DRIVE=1 \
RCLONE_REMOTE=gdrive \
  server/ckd_run_stage2c_susie.sh
```

By default the large chromosome-level PLINK reference downloads are deleted after
all candidate loci on that chromosome have been converted to regional LD data.
Set `STAGE2C_CLEANUP_CHR_CACHE=0` to retain them.

## Outputs

```text
/srv/is-analysis/results/ckd/stage2c_susie/
├── STAGE2C_REGIONS.tsv
├── STAGE2C_LD_QC.tsv
├── STAGE2C_LD_PROVENANCE.json
├── susie_input/<GENE>.tsv.gz
├── susie_results/
│   ├── STAGE2C_SUSIE_DEFAULT.tsv
│   ├── STAGE2C_SUSIE_ALL_PRIORS.tsv
│   ├── <GENE>_pqtl_susie.rds
│   └── <GENE>_egfr_susie.rds
└── SHA256SUMS.txt
```

Large derived LD matrices stay in
`/srv/is-analysis/data/ckd/stage2c_ld/ld` and are not copied to Drive by this
runner. They are reproducible from the public 1000 Genomes reference.
