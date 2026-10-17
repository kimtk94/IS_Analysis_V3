# Five-minute AUTO_NEXT timer — staged, NOT ENABLED

The unit templates in `deploy/systemd/` are inert repository files.
They have not been installed, merged, or activated on the user's server.

## Activation gates

1. PR reviewed and merged, remote checkout updated, CI green.
2. Inspect the existing queue worker service and ensure it uses the same
   shared lock implementation/path on the same filesystem.
3. Check ownership of `/srv/is-analysis/results/is/github_actions_runs`
   and that a shared lock can be acquired by both services.
4. Run `bash scripts/is/auto_next_preflight.sh` and unit tests.
5. Stage one manual `--run` audit and verify state, historical evidence,
   and data completeness.
6. Confirm lock-busy operation is `DEFERRED_LOCK_BUSY` and reattempted on
   the next tick without duplicate audit execution.
7. Only after these checks, install user unit and timer and enable it.
   Do not disable or replace the existing five-minute worker without a rollback.

The timer intentionally invokes only **one** eligible phase per tick.
Technical PASS is not independent scientific validation. Five phases complete
once, and no new research experiments are invented by this timer.

## Rollback

Stop/disable the newly introduced timer only and keep the existing worker
unmodified. Retain state evidence for postmortem. Never delete the state
file to force an automatic retry of uncertain operations.
