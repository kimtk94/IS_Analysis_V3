# Contributor guidance

- Keep source code independent from Google Drive paths; accept roots as CLI options.
- Never commit credentials, downloaded archives, full summary statistics, or results.
- Python code uses the standard library where practical and must expose a CLI `main`.
- Tests must use only the small files under `tests/fixtures`.
- Run `bash scripts/codex_smoke_test.sh` before committing (or report Windows-specific missing tools instead of claiming PASS).

## CKM x stroke research gates
- Keep original phenotype accession mapping AIS=GCST90104540, LAA=GCST90104542, CE=GCST90104541, SVD=GCST90104543; use publisher GWAS Catalog metadata, not conflicting secondary labels.
- Classify original-source accessibility, raw-byte integrity, allele/schema QC, ancestry-matched LD, LDSC, MR/coloc and cell-type annotation as **separate** gates.
- A working URL, successful checksum, or CI SUCCESS is never a scientific result; evidence must cite actual artifacts and statistics.
- GIANT BMI 2018 HRS allele frequency is not cohort-wide EAF; 44 exported BMI P values were numerically zero in the scanned original file, so handle as underflow explicitly.
- Keep SBP OpenGWAS JWT and T2DGGI download-term access behind authenticated user approval. Do not bypass access controls or paste tokens into chats/logs.
- Prefer small, offline tests and bounded GitHub Actions jobs; report costs and avoid duplicate large downloads. Draft PRs require explicit approval to merge.
- For Windows Orca, use the project-local terminal. Never install other dependencies or change PowerShell execution policy without a documented need.
