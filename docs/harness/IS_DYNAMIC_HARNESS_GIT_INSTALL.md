# IS-only dynamic Harness extension

This branch adds an IS-only stage generator for the *standalone* Harness at
`/home/taehoon/agent-harness`. It does not modify Kalman or canonical IS files.

## Installation without Remote Desktop Commander

The remote Ubuntu host must have git access to `kimtk94/IS_Analysis_V3`.
From a server shell, retrieve files **without switching the research worktree branch**:

```bash
git -C /srv/is-analysis/IS_Analysis_V3 fetch origin feat/is-harness-dynamic-loop-20261013
mkdir -p /home/taehoon/agent-harness/is-git-install
git -C /srv/is-analysis/IS_Analysis_V3 show FETCH_HEAD:scripts/harness/is_dynamic_loop.py > /home/taehoon/agent-harness/is-git-install/is_dynamic_loop.py
git -C /srv/is-analysis/IS_Analysis_V3 show FETCH_HEAD:scripts/harness/install_is_dynamic.py > /home/taehoon/agent-harness/is-git-install/install_is_dynamic.py
python3 -m py_compile /home/taehoon/agent-harness/is-git-install/*.py
python3 /home/taehoon/agent-harness/is-git-install/install_is_dynamic.py --apply
```

Check both status and operation, explicitly **scoped to IS**:

```bash
python3 /home/taehoon/agent-harness/experiment_runner.py --job is-research --status
python3 /home/taehoon/agent-harness/experiment_runner.py --job is-research --steps 2 --run
```

## Contract

- Starts after existing stage 24 and uses review/checkpoint pairs
  stage 25/26, stage 27/28, etc.
- Reads only the handoff narrative and previous stage review, copying both
  into a fresh isolated workspace.
- Uses `gpt-5.6-luna` to write a *new* experiment script; does not ask the model
  to execute it.
- Statically checks generated imports and unsafe calls; executes in
  `codex sandbox` workspace-write (never deploy or promote).
- Requires nonempty metrics, limitations, provenance, and at least one new
  metric key to advance.
- Records the first failure and stops; no silent quota retries, no unbounded
  loops, no claim of scientific completion.
- `--steps` limits the number of new stages, max 20 per invocation.

## Known limitations

- Static inspection and `codex sandbox` are **not a formal security boundary**.
  The execution context must be audited before using sensitive datasets.
- Numeric contract validation does **not** establish statistical or causal
  correctness, and the handoff narrative is not primary genotype data.
- Requires a working Codex subscription, available model capacity, and the
  Ubuntu sandbox configuration.
- On AI quota failure, checkpoint remains `WAITING_AI_QUOTA` and requires
  deliberate, audited retry (not automatic deletion).
- Do not use `--job all` for this research-only workflow.

## Rollback

The installer creates
`/home/taehoon/agent-harness/experiment_runner.py.pre_is_dynamic_20261013.bak`.
Copy it back after confirming no Harness process is executing.

This is intentionally a reviewable GitHub-side implementation. Merely merging
this change will NOT install it into the standalone server Harness.
