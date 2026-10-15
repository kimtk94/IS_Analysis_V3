# IS Research GitHub Issue Queue (restricted executor)

ChatGPT can create GitHub Issues with label `is-harness-queued` and exact JSON body such as `{"operation":"status"}`. A server-side polling worker executes only `status`, `inventory`, `tests`, or `mkdir-output`, and writes a bounded result comment back to the Issue. No arbitrary shell, sudo, checkout, deletion, or LIVE trading operations.

## One-time server setup (manual, NOT performed yet)

1. Review `scripts/is/github_issue_queue_worker.py` and install it to a protected server path. Use a dedicated unprivileged account with permissions limited to IS research files; no sudo, Kalman or credentials access.
2. Create a fine-grained GitHub token scoped ONLY to `kimtk94/IS_Analysis_V3` with **Issues: Read and write**; no repository administration or code write permission. Store it in a mode-0600 server-only environment file. **Never paste the token in chat or commit it to Git.**
3. Create labels `is-harness-queued`, `is-harness-running`, `is-harness-done`, `is-harness-failed` in the repository.
4. Run a single poll: `IS_QUEUE_GITHUB_TOKEN=... python3 scripts/is/github_issue_queue_worker.py --once` (load the token from the protected environment, not inline shell history). After verification configure a systemd service with `EnvironmentFile=`, `User=` and `ExecStart=/usr/bin/python3 ... --interval 60`.
5. Open a GitHub Issue authored by `kimtk94` with the queued label and body exactly `{"operation":"status"}`. Verify status comment, then inventory/tests. For write check use `{"operation":"mkdir-output"}`.

## Security and limitations

This is a prototype, **not a hardened distributed queue**. Removing a label is not an atomic claim across workers. Run exactly ONE worker instance; don't enable multiple machines. A crash after label removal leaves a stranded job for manual review. Trusted author verification is not a defense against an account compromise; use branch/issue access controls and token rotation. Issue comments can contain filenames, so use a private repository if results are sensitive. No scientific validity is implied by execution success. This prototype intentionally does not execute arbitrary requested commands or auto-sync branches. Server installation and live operation have not been verified.
