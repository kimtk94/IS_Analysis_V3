# CKD server + Google Drive execution

Recommended storage model:

1. **Git repository / code**: local server filesystem
   - default: `/srv/is-analysis/IS_Analysis_V3`
2. **Download staging + compute**: local server filesystem
   - default: `/srv/is-analysis/data/ckd`
3. **Long-term public-data storage**: Google Drive via `rclone`
   - `IS_Analysis_V3/data/rawdata/ckd`
   - `IS_Analysis_V3/data/analysis_ready/ckd`
4. **Restricted KoGES individual-level data**: do not sync to this public-data Drive tree.

## 1. Clone/update repository

```bash
curl -fsSL https://raw.githubusercontent.com/kimtk94/IS_Analysis_V3/main/server/ckd_server_bootstrap.sh -o /tmp/ckd_server_bootstrap.sh
chmod +x /tmp/ckd_server_bootstrap.sh
/tmp/ckd_server_bootstrap.sh
```

Or manually:

```bash
sudo mkdir -p /srv/is-analysis
sudo chown -R "$USER:$USER" /srv/is-analysis
git clone https://github.com/kimtk94/IS_Analysis_V3.git /srv/is-analysis/IS_Analysis_V3
cd /srv/is-analysis/IS_Analysis_V3
```

For subsequent updates:

```bash
cd /srv/is-analysis/IS_Analysis_V3
git pull --ff-only origin main
```

## 2. Acquire and verify public CKD data locally

```bash
cd /srv/is-analysis/IS_Analysis_V3
CKD_DATA_ROOT=/srv/is-analysis/data/ckd/rawdata \
CKD_READY_ROOT=/srv/is-analysis/data/ckd/analysis_ready \
  server/ckd_acquire_public_data.sh
```

The acquisition script downloads:
- UKB-PPP/Sun 2023 cis instruments
- EUR eGFR, CKD, BUN, eGFRcys and UACR
- EAS eGFR and BUN
- matched analysis-ready datasets

It performs gzip and SHA256/inventory checks before data are considered ready.

## 3. Configure Google Drive with rclone

Check whether a Drive remote already exists:

```bash
rclone listremotes
```

If not:

```bash
rclone config
```

Create a Google Drive remote named `gdrive`.

Test it:

```bash
rclone lsd gdrive:
```

## 4. Sync verified data to Google Drive

```bash
cd /srv/is-analysis/IS_Analysis_V3
RCLONE_REMOTE=gdrive \
LOCAL_CKD_ROOT=/srv/is-analysis/data/ckd \
DRIVE_BASE='IS_Analysis_V3/data' \
  server/ckd_sync_to_drive.sh
```

Resulting Drive layout:

```text
IS_Analysis_V3/
└── data/
    ├── rawdata/
    │   └── ckd/
    │       ├── outcome/
    │       ├── support/
    │       ├── eas/
    │       └── instruments/
    └── analysis_ready/
        └── ckd/
            ├── instruments/
            ├── EUR/
            ├── EAS/
            ├── MATCH_SUMMARY.json
            └── SHA256SUMS.txt
```

## Why stage locally first?

Writing multi-hundred-MB GWAS files directly to a Google Drive FUSE mount can
be slower and more failure-prone. The preferred sequence is:

```text
Internet
  ↓
server local scratch
  ↓ gzip/SHA256 validation
analysis-ready extraction
  ↓
rclone checksum sync
Google Drive
```

This keeps downloads resumable/auditable while Drive remains the persistent
data store.
