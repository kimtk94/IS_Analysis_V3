# IS Harness automatic polling (user-level systemd)

Validated before this file: issue #35 status PASS; #36 inventory PASS; #37 mkdir-output PASS. This does not validate automatic polling.

The queue worker must be installed at `~/is-harness/github_issue_queue_worker.py`, with a GitHub fine-grained token available in `~/.config/is-harness/github.env` as `IS_QUEUE_GITHUB_TOKEN=...` (0600). The worker is authorized only for predefined IS research tasks.

## One-time installation (run as the dedicated authorized **non-root** research account)

```bash
mkdir -p ~/.config/systemd/user ~/is-harness
curl -fL https://raw.githubusercontent.com/kimtk94/IS_Analysis_V3/main/scripts/is/github_issue_queue_worker.py -o ~/is-harness/github_issue_queue_worker.py
curl -fL https://raw.githubusercontent.com/kimtk94/IS_Analysis_V3/main/deploy/is-harness/is-harness-queue.service -o ~/.config/systemd/user/is-harness-queue.service
curl -fL https://raw.githubusercontent.com/kimtk94/IS_Analysis_V3/main/deploy/is-harness/is-harness-queue.timer -o ~/.config/systemd/user/is-harness-queue.timer
python3 -m py_compile ~/is-harness/github_issue_queue_worker.py
chmod 600 ~/.config/is-harness/github.env
systemctl --user daemon-reload
systemctl --user start is-harness-queue.service
journalctl --user -u is-harness-queue.service -n 30 --no-pager
systemctl --user enable --now is-harness-queue.timer
systemctl --user list-timers --all | grep is-harness
```

If `systemctl --user` reports `Failed to connect to bus`, log in through a normal user session and retry rather than running this as root. For unattended operation after logout or reboot, an administrator must explicitly approve `sudo loginctl enable-linger <research-user>` (not configured by this repository).

Stop using `systemctl --user disable --now is-harness-queue.timer`. Check logs via `journalctl --user -u is-harness-queue.service -n 100 --no-pager`. Do not run a second long-running worker: the queue claim operation is not transactional. The systemd service uses `flock` to limit concurrent invocation of this service, but cannot guard an unrelated manually launched worker. If a job is stranded in `is-harness-running`, investigate and restore manually only after confirming it is not active.

**Important:** All jobs report limited outputs to GitHub Issues. Keep the repository private for sensitive data, and never put API tokens or sample-level private metadata in an Issue.
