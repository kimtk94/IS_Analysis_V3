# CKD thesis data acquisition

Data acquisition precedes downstream MR/colocalization workflow development.

## Acquisition order

1. EUR outcomes: CKDGen Stanzick 2021 eGFRcrea and Wuttke 2019 binary CKD.
2. Renal support traits: BUN, eGFRcys, then UACR.
3. EAS outcome/support: Chen 2024 Taiwan + Japan eGFR and BUN meta-GWAS.
4. Plasma pQTL screening instruments: UKB-PPP cis/independent pQTL signals.
5. Candidate-only UKB-PPP full cis-region summary statistics for colocalization.
6. KoGES individual-level data only after CODA approval; never commit or publicly upload restricted data.

## Pinned Chen 2024 EAS files

- eGFR: `TWB2_BBJ_eGFR_hg19_METAL_FUMA_noNA.gz`, Figshare file 42774049, 111,041,990 bytes.
- BUN: `TWB2_BBJ_BUN_hg19_METAL_FUMA_noNA.gz`, Figshare file 43218600, 132,381,963 bytes.

The seed download workflow fetches both files completely, validates gzip integrity,
creates SHA256 sidecars/manifests, and keeps a short-lived GitHub Actions artifact for transfer testing.

## Why UKB-PPP raw data are not bulk-downloaded first

The existing `ukb_ppp_download_manifest.tsv` contains per-protein archives that are commonly around 0.5 GB each. Bulk-downloading thousands of proteins would therefore create a TB-scale staging requirement. The thesis pipeline first uses compact cis-pQTL/independent-signal instruments to screen proteins by MR. Full regional pQTL data are acquired only for candidates that proceed to colocalization.

## Commands

List core public datasets:

```bash
python3 scripts/download_ckd_public_data.py \
  --dest "$WORK_ROOT/data/rawdata/ckd" --public-core --list
```

Download EUR eGFR and binary CKD:

```bash
python3 scripts/download_ckd_public_data.py \
  --dest "$WORK_ROOT/data/rawdata/ckd" \
  --data-id ckdgen_stanzick2021_egfr_eur \
  --data-id ckdgen_wuttke2019_ckd_eur
```

Download the pinned EAS eGFR and BUN files:

```bash
python3 scripts/download_ckd_public_data.py \
  --dest "$WORK_ROOT/data/rawdata/ckd" \
  --data-id eas_chen2024_egfr_meta \
  --data-id eas_chen2024_bun_meta
```

Every downloaded file is gzip-tested when applicable and gets a sidecar `.download.json` containing source URL, byte size, SHA256, phenotype, ancestry, genome build and acquisition timestamp.
