#!/usr/bin/env python3
"""IS-only autonomous Harness research loop; never edits canonical IS or Kalman data.

Install into /home/taehoon/agent-harness/is_dynamic_loop.py.
Usage: python3 is_dynamic_loop.py --steps 5 --run
Historical stages 18/20/22/24 are not modified. New reviews are 26,28,... .
Generated code and data remain confined to versioned temporary workspaces.
"""
from __future__ import annotations
import argparse
import ast
import datetime as dt
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HARNESS = Path("/home/taehoon/agent-harness")
STAGES = HARNESS / "tasks/is-research/stages"
SOURCE = Path("/srv/is-analysis/worktrees/is-broad-discovery-20261009/docs/harness/IS_RESEARCH_HANDOFF_20261012.md")
WORK = HARNESS / "workspaces"
FIRST_CYCLE = 5
MAX_STEPS = 20
MODEL = "gpt-5.6-luna"
ALLOWED_IMPORTS = {"json", "math", "re", "csv", "statistics", "collections", "itertools", "functools", "pathlib", "datetime", "decimal", "typing"}
BLOCKED_CALLS = {"eval", "exec", "compile", "__import__", "input", "breakpoint"}
BLOCKED_TOKENS = ("/opt/kalman", "/srv/is-analysis", "/home/taehoon/.config", "/home/taehoon/agent-harness/tasks")

def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def atomic_json(path, obj):
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    os.replace(temp, path)

def read_json(path):
    return json.loads(path.read_text())

def review_path(cycle):
    return STAGES / ("stage_%03d_autonomous_review.json" % (16 + 2 * cycle))

def work_record_path(cycle):
    return STAGES / ("stage_%03d_autonomous.json" % (15 + 2 * cycle))

def subprocess_result(command, cwd, timeout):
    try:
        proc = subprocess.run(command, cwd=cwd, stdin=subprocess.DEVNULL,
                              text=True, capture_output=True, timeout=timeout,
                              env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
        return {"exit": proc.returncode, "stdout_tail": proc.stdout[-2500:],
                "stderr_tail": proc.stderr[-2000:]}
    except subprocess.TimeoutExpired:
        return {"exit": 124, "stderr_tail": "TIMEOUT"}

def source_gate(script):
    if not script.is_file() or script.is_symlink() or script.stat().st_size > 35000:
        return "MISSING_OR_OVERSIZE_SCRIPT"
    code = script.read_text()
    try:
        tree = ast.parse(code)
    except SyntaxError:
        return "INVALID_PYTHON"
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            names = ([x.name.split(".")[0] for x in node.names]
                     if isinstance(node, ast.Import)
                     else [node.module.split(".")[0] if node.module else ""])
            if any(name not in ALLOWED_IMPORTS for name in names):
                return "DISALLOWED_IMPORT"
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            if node.func.id in BLOCKED_CALLS:
                return "DISALLOWED_CALL"
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            if any(token in node.value for token in BLOCKED_TOKENS):
                return "EXTERNAL_SOURCE_REFERENCE"
    return None

def validate_result(result, previous):
    if not isinstance(result, dict) or result.get("status") != "ANALYSIS_COMPLETE":
        return "INVALID_RESULT_STATUS"
    metrics = result.get("metrics")
    if not isinstance(metrics, dict) or not metrics:
        return "EMPTY_METRICS"
    if not isinstance(result.get("limitations"), list) or not result["limitations"]:
        return "LIMITATIONS_MISSING"
    if not isinstance(result.get("source_files"), list) or not result["source_files"]:
        return "PROVENANCE_MISSING"
    current = set(metrics)
    prior = set(previous.get("metric_keys", []))
    if prior and current and len(current - prior) == 0:
        return "NO_NEW_METRIC_KEYS"
    try:
        json.dumps(result, allow_nan=False)
    except (ValueError, TypeError):
        return "NONFINITE_OR_NONJSON_RESULT"
    return None

def next_cycle():
    for cycle in range(FIRST_CYCLE, FIRST_CYCLE + 10000):
        if not review_path(cycle).exists():
            return cycle
    raise RuntimeError("cycle index exhausted")

def run_one(cycle):
    previous_file = review_path(cycle - 1)
    if not previous_file.exists():
        return {"cycle": cycle, "status": "BLOCKED_PREVIOUS_REVIEW_MISSING"}
    prior = read_json(previous_file)
    if prior.get("status") != "AUTONOMOUS_EXPERIMENT_REVIEW_REQUIRED":
        return {"cycle": cycle, "status": "BLOCKED_PREVIOUS_REVIEW", "prior_status": prior.get("status")}
    if not SOURCE.is_file() or SOURCE.is_symlink() or SOURCE.stat().st_size > 3_000_000:
        return {"cycle": cycle, "status": "BLOCKED_HANDOFF_UNAVAILABLE"}
    workspace = Path(tempfile.mkdtemp(prefix="is-cycle-%03d-" % cycle, dir=WORK))
    data = workspace / "data"
    data.mkdir()
    shutil.copyfile(SOURCE, data / "handoff.md")
    (data / "handoff.md").chmod(0o444)
    feedback = {
        "status": prior["status"],
        "metric_keys": prior.get("metric_keys", []),
        "metrics": prior.get("result", {}).get("metrics", {}),
        "limitations": prior.get("result", {}).get("limitations", []),
    }
    atomic_json(data / "previous_findings.json", feedback)
    (data / "previous_findings.json").chmod(0o444)
    prompt = (
        "Create exactly experiment.py in this disposable workspace. Do not run it. "
        "Use ONLY ./data/handoff.md and ./data/previous_findings.json. "
        "Implement a deterministic NEW quantitative IS evidence diagnostic with Python standard library. "
        "Write result.json containing status='ANALYSIS_COMPLETE', nonempty metrics dict, "
        "nonempty limitations list, and source_files list. "
        "Read previous_findings.json, do NOT repeat all prior metric keys. "
        "Clearly distinguish summary-derived exploration from genuine genotype-level analysis. "
        "No internet, subprocess, imports beyond json/math/re/csv/statistics/collections/"
        "itertools/functools/pathlib/datetime/decimal/typing, no secret access, no mutation outside workspace. "
        "Do not claim causal inference, formal coloc, or validation without primary data. "
        "If the narrative lacks required values, record missingness instead of inventing numbers."
    )
    agent = subprocess_result(
        ["codex", "exec", "--ephemeral", "--skip-git-repo-check",
         "-m", MODEL, "-s", "workspace-write", "-C", str(workspace), prompt],
        workspace, 300
    )
    script = workspace / "experiment.py"
    failure = "AI_GENERATION_FAILED" if agent["exit"] != 0 else source_gate(script)
    execution = {"exit": None, "stderr_tail": "NOT_RUN"}
    result = {}
    if failure is None:
        execution = subprocess_result(
            ["codex", "sandbox", "-c", 'sandbox_mode="workspace-write"',
             "--", sys.executable, "-I", "-S", "experiment.py"], workspace, 90
        )
        if execution["exit"] == 0:
            try:
                result = read_json(workspace / "result.json")
            except (OSError, ValueError):
                failure = "MISSING_OR_INVALID_OUTPUT"
        else:
            failure = "SANDBOX_EXECUTION_FAILED"
    if failure is None:
        failure = validate_result(result, prior)
    quota = "usage limit" in agent.get("stderr_tail", "").lower()
    state = ("AUTONOMOUS_EXPERIMENT_REVIEW_REQUIRED" if failure is None
             else "WAITING_AI_QUOTA" if quota else "AUTONOMOUS_EXPERIMENT_BLOCKED")
    record = {
        "job": "is-research", "cycle": cycle, "status": state,
        "error": failure, "at_utc": utc_now(), "workspace": str(workspace),
        "source_snapshot_sha256": hashlib.sha256((data / "handoff.md").read_bytes()).hexdigest(),
        "agent": agent, "execution": execution, "result": result if failure is None else {},
        "metric_keys": list(result.get("metrics", {})) if failure is None else [],
        "scientific_completion": False, "live_modified": False,
        "source_modified": False, "promoted": False,
    }
    atomic_json(work_record_path(cycle), {"cycle": cycle, "workspace": str(workspace),
               "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest() if script.is_file() else None,
               "source_sha256": record["source_snapshot_sha256"], "agent_exit": agent["exit"]})
    atomic_json(review_path(cycle), record)
    return {"cycle": cycle, "status": state, "error": failure, "metrics": record["metric_keys"][:12]}

def main():
    p = argparse.ArgumentParser(description="IS-only bounded autonomous experimental loop")
    p.add_argument("--steps", type=int, default=1)
    p.add_argument("--run", action="store_true")
    p.add_argument("--status", action="store_true")
    args = p.parse_args()
    if not 1 <= args.steps <= MAX_STEPS:
        p.error("--steps must be 1..20")
    STAGES.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    if args.status or not args.run:
        c = next_cycle()
        print(json.dumps({"job": "is-research", "next_cycle": c,
                          "last_review": str(review_path(c - 1)), "running": False,
                          "ready": review_path(c - 1).exists()}, ensure_ascii=False))
        return
    with (STAGES / "is_dynamic_loop.lock").open("a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print(json.dumps({"state": "SKIPPED_RUNNING", "job": "is-research"}))
            return
        results = []
        for _ in range(args.steps):
            r = run_one(next_cycle())
            results.append(r)
            if r["status"] != "AUTONOMOUS_EXPERIMENT_REVIEW_REQUIRED":
                break
        print(json.dumps({"job": "is-research", "new_steps": len(results), "results": results,
                          "research_complete": False}, ensure_ascii=False))

if __name__ == "__main__":
    main()
