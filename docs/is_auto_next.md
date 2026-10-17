# IS AUTO_NEXT v1 — controlled execution

AUTO_NEXT reduces manual "next" prompts by executing **one approved research audit per invocation**. Existing `github_issue_queue_worker.py` and its scheduled refresh must remain installed and running independently.

## Review-first rollout

```bash
cd /srv/is-analysis/IS_Analysis_V3
python3 scripts/is/auto_next.py
python3 scripts/is/auto_next.py --state-dir /tmp/is-auto-next-smoke
python3 scripts/is/auto_next.py --run --state-dir /tmp/is-auto-next-smoke
```

The final command genuinely executes one audit: use only after reviewing source code and confirming access to data. Do **not** run the `--run` command merely to inspect status.

## Safety and semantics

- No shell commands from GitHub Issue descriptions, prompts or agent output.
- Only operations from the existing worker's fixed allowlist, with a conservative P1–P5 dependency order.
- `--run` executes **at most one** step; use a systemd timer only after staging verification. Default is dry run.
- State is local atomic JSON. Linux advisory flock prevents concurrent invocation of this controller; it does not lock the separate existing worker, so coordinate maintenance windows before enabling.
- A crashed `RUNNING` task and any `REVIEW_REQUIRED` task stop automatic advancement. Do not silently reset.
- `PASS` means the operation exited successfully, **not** that a scientific claim is validated.
- Tasks are one-shot. Adding new analysis goals or resetting state requires an explicit reviewed change.
- No automated Git merge, server deployment, privileged access, database deletion, or LIVE trading.
- For recurring automation, independently verify input freshness and create a *new reviewed task series* rather than resetting this state file.

## OpenClaw integration (phase 2)

OpenClaw, if installed in an isolated low-privilege environment, may read `state.json` and ask for approval when `REVIEW_REQUIRED` occurs. It should **not** have direct production SSH / credentials or unrestricted command execution. Operate the approved controller using a dedicated least-privilege service account.
