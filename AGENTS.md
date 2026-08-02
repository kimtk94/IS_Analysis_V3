# Contributor guidance

- Keep source code independent from Google Drive paths; accept roots as CLI options.
- Never commit credentials, downloaded archives, full summary statistics, or results.
- Python code uses the standard library where practical and must expose a CLI `main`.
- Tests must use only the small files under `tests/fixtures`.
- Run `bash scripts/codex_smoke_test.sh` before committing.
