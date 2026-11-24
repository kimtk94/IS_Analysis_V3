# IS Harness automatic polling (user-level systemd)

The queue worker reports technical execution separately from scientific validity.
An old successful Issue does not establish that the current worker is running.

The queue worker must be installed at `~/is-harness/github_issue_queue_worker.py`, with a GitHub fine-grained token available in `~/.config/is-harness/github.env` as `IS_QUEUE_GITHUB_TOKEN=...` (0600). The worker is authorized only for predefined IS research tasks.

## One-time installation (run as the dedicated authorized **non-root** research account)

From the reviewed local checkout, select an immutable commit. The installer reads
all four Python modules and both units from that same commit, validates the imports,
and writes a SHA256 receipt. It holds installation while the queue service lock is
busy and does not start services, enable timers, or access the GitHub API.

```bash
REVISION="$(git rev-parse HEAD)"
python3 deploy/is-harness/install_worker.py --revision "$REVISION"
cat ~/is-harness/installed_revision.json
```

The installed modules are `github_issue_queue_worker.py`, `harness_research_jobs.py`,
`harness_science_pipeline.py`, and `harness_archived_evidence.py`. Local uncommitted
changes are deliberately excluded from installation. Configure the token file
separately before activating the service:

```bash
chmod 600 ~/.config/is-harness/github.env
systemctl --user daemon-reload
systemctl --user start is-harness-queue.service
journalctl --user -u is-harness-queue.service -n 30 --no-pager
systemctl --user enable --now is-harness-queue.timer
systemctl --user list-timers --all | grep is-harness
```

If `systemctl --user` reports `Failed to connect to bus`, log in through a normal user session and retry rather than running this as root. For unattended operation after logout or reboot, an administrator must explicitly approve `sudo loginctl enable-linger <research-user>` (not configured by this repository).

Stop using `systemctl --user disable --now is-harness-queue.timer`. Check logs via `journalctl --user -u is-harness-queue.service -n 100 --no-pager`. Do not run a second long-running worker: the queue claim operation is not transactional. The systemd service uses `flock` to limit concurrent invocation of this service, but cannot guard an unrelated manually launched worker. If a job is stranded in `is-harness-running`, investigate and restore manually only after confirming it is not active.

## Stale requests and validation

A queued request older than one hour, or a running request without an update for
one hour, produces `HOLD_STALE_REQUEST`. The scheduled GitHub workflow fails
visibly for a stale research refresh; the worker's `--once` exits with status 2.
Neither path cancels, deletes, or requeues the old request. Inspect the worker
logs, installed revision, service status, and any active research process before
an operator decides how to recover it. These checks do not prove that an old job
has stopped. Completed results must be reconciled before a manual retry.
Missing, invalid, or future request timestamps also hold execution. Queue
inventory is capped at 100 items per pending label; hitting that limit holds
instead of silently ignoring older requests. Pending-state checks run even when
the local input fingerprint is unchanged.

Each service activation executes at most one queued operation. `operation=tests`
runs the maintained `bash scripts/codex_smoke_test.sh` entrypoint in the configured
research checkout; its fixture-only baseline also works on a clean checkout.

All jobs report limited outputs to GitHub Issues. Never include API tokens or
sample-level private metadata in an Issue; this repository may be public.
