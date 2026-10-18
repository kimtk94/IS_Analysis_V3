# IS Harness: GitHub Actions → SSH server runbook

Prepared 2026-10-18. The GitHub Actions workflow `.github/workflows/is-harness-ssh-p0.yml` is **committed**, but there has been **no successful SSH run or research server P0 execution** confirmed yet.

## What it does
`workflow_dispatch` offers two modes:
- `probe` (default): checks four Actions secrets, validates the SSH private key without displaying it, connects via `BatchMode=yes` and `StrictHostKeyChecking=yes`, and checks server data sources and numpy/pandas. **Read-only research server inspection.**
- `run`: after passing probe, stages exactly four pinned source scripts from the workflow-triggering Git commit into an isolated server `/tmp` directory, then executes P0.1, P0.2, P0.3 and P0.4 in order. Every run writes new directories suffixed by the GitHub Actions run ID and attempt number. Failure immediately stops later stages. It never touches production LIVE settings or the current research worktree.

## Prerequisites
The repository Settings → Secrets and variables → Actions must contain `SSH_HOST`, `SSH_USER`, `SSH_KEY` (a readable SSH private key with real newlines), and `SSH_KNOWN_HOSTS` (pin the SSH host key fingerprint independently, **never** collect it and trust it blindly inside the same job). Provision only a narrowly scoped dedicated server SSH identity, preferably limited to research output directories, not an unrestricted personal administrator key. Server egress and inbound SSH firewall must allow the selected GitHub runner.

GitHub normally requires a manually dispatchable workflow file to be present on the **default branch**. The committed file currently sits on `research/is-broad-discovery-20261009`. Review and merge/cherry-pick this workflow to the default branch, or use an existing approved default-branch workflow; do not assume the Run workflow button is visible for a feature-branch-only workflow. After the workflow is on the default branch, select the **research branch** in the Actions UI before initiating `probe` so code is sourced from the intended commit. Restrict untrusted PR branches from secrets.

## Operator procedure
1. GitHub repository → Actions → `IS Harness SSH P0 (manual, non-LIVE)`.
2. Choose Run workflow → Branch: `research/is-broad-discovery-20261009`, mode: `probe`. Check log lines `SSH_AUTH=PASS`, `SOURCE_INPUTS=PASS`, `PYTHON_DEPENDENCIES=PASS`. Never print private keys.
3. Only when probe is successful, rerun with mode `run`. The job prints `P0_AUDIT_STATUS=COMPLETE` and `P0_SCORECARD=...` if all four scripts finish.
4. Do **not** mark P0 complete until the resulting manifest/TSV have been reviewed for expected three QTL SNPs, pairwise LD matrix, four gene comparators, and scientific boundary flags. Test and CI should also be verified. The workflow currently validates script execution but **does not run the full test_is_*.py regression suite on the server**.

## Operational cautions
- Original 2,225 gene discovery universe, 646 historical coloc, and 120 ABF sensitivity results must remain unchanged.
- Research outputs are written to new versioned paths, not to original canonical data, but `run` does write research result directories. Only `probe` is fully read-only.
- Server scripts use the run-triggering repository commit; failures should be inspected from Actions logs and fixed via Git commit, not patched manually in place.
- This is GitHub Actions-based remote execution, **not proof that the Codex cloud execution environment itself has an SSH private key or can SSH**. Both are potential execution paths and have separate credentials/network rules.
- An exhausted Remote Desktop Commander plugin quota does not necessarily restrict the independent GitHub Actions job.

Current status: `WORKFLOW_COMMITTED / SECRETS_UNVERIFIED / NO_SSH_RUN_CONFIRMED / P0_RUN_NOT_COMPLETED`.
