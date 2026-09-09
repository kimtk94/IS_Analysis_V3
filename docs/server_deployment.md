# Personal server deployment (v1 bootstrap)

This deployment converts the Colab/Google Drive execution model into a persistent
Ubuntu server model while preserving the V3 rule that code and runtime data are
separate.

## Capacity policy for the current server

Observed root filesystem at bootstrap planning time:

- root LV: about 231 GiB
- used: about 17 GiB
- available: about 202 GiB

The server is therefore suitable for batched UKB-PPP/GIGASTROKE processing, but
not for keeping a permanent mirror of all source archives.

Hard defaults:

- reserve at least 40 GiB free space
- cap `data/staging` at 45 GiB
- only one raw download at a time
- delete raw archives only after processed output, QC, and manifest/checksum are verified
- keep failed/unverified raw files for diagnosis
- never back up raw staging archives to Google Drive

These limits live in `config/server.yaml` and can be changed after measuring real
batch sizes.

## 1. Confirm machine resources

Run these before tuning worker counts:

```bash
free -h
nproc
lscpu | egrep 'Model name|CPU\(s\)|Thread|Core|Socket'
lsblk -o NAME,SIZE,FSTYPE,MOUNTPOINTS
sudo pvs
sudo vgs
sudo lvs
```

The current configuration intentionally starts conservatively with one download
and at most two analysis jobs.

## 2. Clone the bootstrap branch

```bash
cd ~
git clone https://github.com/kimtk94/IS_Analysis_V3.git
cd IS_Analysis_V3
git switch server-v1-bootstrap
```

After the bootstrap branch is merged, use `main` instead.

## 3. Create the persistent runtime layout

```bash
sudo bash server/bootstrap.sh
```

Expected root:

```text
/srv/is-analysis/
  app -> repository checkout
  data/
    staging/
    rawdata/
    reference/
    standardized/
  results/
  reports/
  manifests/
  logs/
  state/
  backup/
```

Re-check the storage guard at any time:

```bash
python3 scripts/storage_guard.py \
  --path /srv/is-analysis \
  --staging /srv/is-analysis/data/staging \
  --min-free-gib 40 \
  --max-staging-gib 45
```

No download job should start when this command fails.

## 4. Server-only secrets

```bash
cp server/.env.example server/.env
chmod 600 server/.env
nano server/.env
```

Set a long random PostgreSQL password and configure the Synapse token only on the
server. Do not commit `.env`, Synapse credentials, rclone credentials, source
archives, full summary statistics, or analysis results.

## 5. Docker services

Install Docker Engine/Compose using the Ubuntu-supported installation method,
then run:

```bash
cd server
docker compose build
docker compose up -d postgres pipeline
docker compose ps
```

PostgreSQL is intentionally not published to a host TCP port. It is reachable
only on the Compose network unless the configuration is explicitly changed.

## 6. Migrate reviewed V2/V3 inputs

Before production download, populate at minimum:

```text
/srv/is-analysis/data/metadata/ukb_ppp_download_manifest.tsv
/srv/is-analysis/data/reference/gene_coordinates_hg38.tsv
/srv/is-analysis/data/rawdata/outcome/... GIGASTROKE inputs or catalog metadata
```

Use the reviewed Google Drive V2/V3 files as the source. Do not substitute the
synthetic test fixtures for production data.

## 7. Validation order

Do not schedule unattended production runs until these gates pass in order:

1. storage guard
2. repository smoke test
3. Synapse authentication/download of one known source
4. IDO1 bounded EUR/EAS test
5. canonical schema and row-count QC
6. raw-delete-after-verified test
7. Google Drive report/manifest backup test
8. 10-15 gene smoke batch
9. only then enable recurring execution

## Storage layout rationale

The 231 GiB root volume should be treated as compute workspace, not archival
storage. Persistent high-value products are standardized datasets, manifests,
checksums, MR/coloc outputs, reports, and the reference data required for
reproducibility. Re-downloadable source archives are staging material.

If future standardized data plus LD references approach the 40 GiB reserve,
attach another SSD/HDD or expand the LVM before increasing batch concurrency.
