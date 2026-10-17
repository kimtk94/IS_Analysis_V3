# Metabolic Resilience pipeline scripts

This directory is reserved for the reproducible server-side scripts supporting
the `METABOLIC_RESILIENCE_MASTER` workflow.

## CI contract

GitHub Actions runs:

1. repository smoke tests,
2. `bash -n` syntax validation,
3. syntax parsing of embedded Python heredocs,
4. ShellCheck,
5. project-policy checks.

The CI intentionally rejects:

- `set -e` / `errexit`,
- obvious hard-coded credentials,
- a fixed `ROOT=/srv/is-analysis` assignment.

Scripts should accept a configurable analysis root:

```bash
ROOT="${IS_ANALYSIS_ROOT:-/srv/is-analysis}"
```

## Planned Stage 3 sequence

```text
C0R2 -> C1 -> C2 -> D0 -> D1A -> D1B -> D2 -> D3 -> E0 -> E1
```
