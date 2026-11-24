# IS Harness automatic polling (user-level systemd)

The queue worker reports technical execution separately from scientific validity.
An old successful Issue does not establish that the current worker is running.

The queue worker must be installed at `~/is-harness/github_issue_queue_worker.py`, with a GitHub fine-grained token available in `~/.config/is-harness/github.env` as `IS_QUEUE_GITHUB_TOKEN=...` (0600). The worker is authorized only for predefined IS research tasks.

## One-time installation (run as the dedicated authorized **non-root** research account)

From the reviewed local checkout, select an immutable commit. The installer reads
all four Python modules and both units from that same commit, validates the imports,
and writes a SHA256 receipt. It holds installation while the queue service or
shared research execution lock is busy. It does not start services, enable timers,
or access the GitHub API. The CLI expects the existing execution lock at
`/srv/is-analysis/results/is/github_actions_runs/is_research_execution.lock`;
use `--execution-lock` for the actual configured path when it differs.

Inspect the runtime first using the read-only audit below. Stop the queue and
AUTO_NEXT timers if present, and let active work finish before installation.
Stopping a timer does not establish that its service or research process has
finished. Do not kill an active job to make installation proceed.

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

An enabled timer that is currently stopped can run after a reboot or a later
operator start. Installation alone does not certify queue recovery or authorize
a scientific analysis.

## GitHub deployment workflow

`IS Harness Safe Research Worker Update` now validates code on pushes and pull
requests. **A merge does not install or start the worker.** Use its manual
`workflow_dispatch` with `operation=install` to install the selected revision.
The default operation is `validate`.

The installation job exports a small bundle from the exact workflow commit:
four modules, two units, the installer, and a SHA256 manifest. It verifies the
transferred installer before execution, checks the expected commit and every
payload checksum, validates imports and unit contracts, and installs under both
execution locks. A failed replacement restores the prior files. It does not
fetch, reset, or change the server's research checkout.

The remote job holds before replacement if a queue/AUTO_NEXT service or timer
is active, its state cannot be read, or a related process is still present. It
uses the canonical shared execution-lock path above. For a custom lock layout,
use the local installer with the matching explicit path. Successful installation
reloads the unit definitions but does not start services or enable timers;
installation or reload errors fail the workflow.

For offline transport, export from the reviewed checkout:

```bash
python3 deploy/is-harness/install_worker.py \
  --repo-root "$PWD" --revision "$REVISION" --export-bundle /tmp/is-harness-bundle
```

The export directory must not already exist. Transfer it over the established
verified channel, verify the installer digest from that channel, and install with
`--bundle-dir <directory> --revision <expected-full-SHA>`. The expected SHA must
come from the reviewed revision, rather than from an unverified bundle alone.

## Read-only worker diagnosis

Run `IS Server Read-Only SSH Audit` with `audit=worker-status` on the revision
being reviewed. It streams the diagnostic through SSH and does not install
anything on the server. It reports only bounded operational metadata: installed
file hashes and receipt consistency, checkout revision, service/timer states,
lock occupancy, process counts, heartbeat/result freshness, and journal error
categories. It does not read the token environment file or print raw logs,
process command lines, or research rows.

This diagnostic does not import or run the worker, start/stop a service, change
an Issue, or create a missing lock file. An unknown or inconsistent runtime is a
visible hold, not a successful health check. A oneshot service can be inactive
between timer invocations; assess it together with the timer and heartbeat.

Its JSON separates `collection_status` (`COMPLETE` or `PARTIAL`) from runtime
`status` (`OPERATIONAL_OK` or `ATTENTION_REQUIRED`). Exit status 2 still emits the
sanitized report so an unavailable systemd manager, stale heartbeat, or
inconsistent installation can be diagnosed. Even an operationally healthy
snapshot is not authorization to retry an old request.

For an existing authorized server shell, the same check can be streamed without
copying a script to the host:

```bash
ssh <existing-verified-research-target> 'python3 -I -B -' \
  < scripts/is/audit_queue_runtime.py
```

Use the existing pinned-host SSH configuration. For nonstandard layouts, pass
the diagnostic's explicit `--repo-root`, `--install-root`, `--unit-root`,
`--runtime-dir`, and `--results-root` arguments. It observes the configured
default execution lock under the results root; it does not infer overrides from
the credential environment file. A custom execution-lock path must be
reconciled separately before recovery or installation.

The server research checkout is a separate deployment target. Compare its
reported full SHA with the intended revision before running `operation=tests`;
that operation tests the checkout, not the installed home-directory bundle.

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
