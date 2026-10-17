# AUTO_NEXT staged server rollout — no automatic activation

The existing `github_issue_queue_worker.py` performs a five-minute scheduled refresh. **Do not schedule AUTO_NEXT concurrently** until an operator reviews process and resource concurrency. The controller's flock protects only its own state; it is not a shared lock with the existing worker.

## Stage A: inspect without changing running services

From a reviewed checkout containing this PR:

```bash
bash scripts/is/auto_next_preflight.sh
python3 -m unittest discover -s tests -p 'test_is_auto_next.py' -v
python3 scripts/is/auto_next.py --state-dir /tmp/is-autonext-staging
```

Expected: preflight static PASS, 7 unit tests pass, dry-run says P1 inventory; no `state.json` created in dry run. Inspect all active user/system units and timer triggers (including GitHub Actions) separately.

## Stage B: staged controlled operation

After confirming there is no concurrent operation with the GitHub Issue worker, run **one** manual step using a dedicated staging state directory:

```bash
python3 scripts/is/auto_next.py --run --state-dir /tmp/is-autonext-staging
cat /tmp/is-autonext-staging/state.json
```

Review the recorded outputs. Continue only when technical PASS plus scientific interpretation are both accepted. No current script grants scientific approval automatically.

## Stage C: 5-minute orchestration (not yet enabled)

First migrate the existing refresh and AUTO_NEXT through a **shared execution lock** or a single orchestrator. Configure unit scheduling only after proving no overlapping `research-refresh` executions, and provide a rollback to the current worker. Do not disable the existing worker or use a parallel timer as a shortcut.

## Exit gates

- Verify real server working tree and active timers.
- Unit tests run against the deployed revision.
- Confirm shared lock / one orchestrator across both controllers.
- Record first real run start and finish, plus clean rollback.
- Keep live trading, secrets, git merging, and privileged actions outside scope.
