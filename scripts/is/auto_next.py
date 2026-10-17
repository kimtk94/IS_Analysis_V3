#!/usr/bin/env python3
"""Conservative AUTO_NEXT controller for the existing IS research harness.

One bounded task per invocation. No LLM-originated shell, live trading, or auto merge.
Requires explicit --run; otherwise produces a read-only preview.
"""
import argparse
import datetime as dt
import fcntl
import json
import os
from pathlib import Path
import tempfile
import time

from github_issue_queue_worker import run_job

DEFAULT_ROOT = Path("/srv/is-analysis/results/is/github_actions_runs/auto_next")
ALLOWED = {"status", "inventory", "tests", "headers", "eas5-audit",
           "candidate-2225-audit", "figure-manifest", "research-refresh"}
STEPS = [
    {"id": "P1-inventory", "phase": "P1", "operation": "inventory", "requires": []},
    {"id": "P2-headers", "phase": "P2", "operation": "headers", "requires": ["P1-inventory"]},
    {"id": "P3-eas5", "phase": "P3", "operation": "eas5-audit", "requires": ["P2-headers"]},
    {"id": "P4-figure", "phase": "P4", "operation": "figure-manifest", "requires": ["P3-eas5"]},
    {"id": "P5-refresh", "phase": "P5", "operation": "research-refresh", "requires": ["P4-figure"]},
]

def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def write_atomic(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".autonext-", suffix=".json", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as out:
            json.dump(data, out, indent=2, ensure_ascii=False)
            out.write("\n")
            out.flush()
            os.fsync(out.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)

def load(path):
    if not path.is_file():
        return {"schema": "is-auto-next-v1", "tasks": {}, "history": []}
    obj = json.loads(path.read_text(encoding="utf-8"))
    if obj.get("schema") != "is-auto-next-v1" or not isinstance(obj.get("tasks"), dict):
        raise ValueError("Unknown or invalid state schema; refusing to overwrite")
    return obj

def choose(state):
    statuses = state["tasks"]
    for step in STEPS:
        old = statuses.get(step["id"], {})
        status = old.get("status")
        if status == "PASS":
            continue
        if status in {"REVIEW_REQUIRED", "RUNNING", "BLOCKED"}:
            return step, "BLOCKED_" + status
        if any(statuses.get(dep, {}).get("status") != "PASS" for dep in step["requires"]):
            return step, "BLOCKED_DEPENDENCY"
        if old.get("attempts", 0) >= 2:
            return step, "BLOCKED_RETRY_LIMIT"
        return step, "READY"
    return None, "ALL_TASKS_COMPLETED"

def validate_output(operation, output):
    """Minimal evidence checks; successful subprocess exit alone is insufficient.

    This only confirms expected audit artifacts appear in output. It does not
    establish scientific validity.
    """
    if not isinstance(output, str):
        raise ValueError("Non-text audit result")
    if operation == "headers":
        rows = output.splitlines()
        if any(row.startswith("MISSING\t") for row in rows):
            raise ValueError("Missing required header input")
        if sum(row.startswith("FILE\t") for row in rows) < 4:
            raise ValueError("Incomplete required header inventory")
    elif operation == "inventory":
        matches = [line for line in output.splitlines()
                   if line.strip().endswith(("IS_EAS5", ".tsv", ".tsv.gz", ".json"))]
        if not matches:
            raise ValueError("Empty EAS5 input inventory")
    return True


def execute(state_path, run=False):
    lock_path = state_path.with_suffix(".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return {"status": "ALREADY_RUNNING"}
        state = load(state_path)
        step, decision = choose(state)
        if not step:
            return {"status": decision, "tasks": state["tasks"]}
        if decision != "READY":
            return {"status": decision, "task": step["id"]}
        if step["operation"] not in ALLOWED:
            return {"status": "BLOCKED_NOT_ALLOWLISTED", "task": step["id"]}
        if not run:
            return {"status": "DRY_RUN", "next": step, "state_path": str(state_path)}
        previous = dict(state["tasks"].get(step["id"], {}))
        record = {"status": "RUNNING", "operation": step["operation"],
                  "phase": step["phase"], "attempts": previous.get("attempts", 0) + 1,
                  "started_at": now(), "ended_at": None}
        state["tasks"][step["id"]] = record
        write_atomic(state_path, state)  # crash leaves RUNNING for human review
        started = time.monotonic()
        try:
            output = run_job(step["operation"])
            validate_output(step["operation"], output)
            record["status"] = "PASS"  # technical execution only, not scientific validation
            record["output_excerpt"] = str(output)[-3000:]
        except RuntimeError as exc:
            if str(exc) == "SHARED_EXECUTION_LOCK_BUSY":
                # A collision is not an audit attempt. Try again next invocation.
                if previous:
                    state["tasks"][step["id"]] = previous
                else:
                    state["tasks"].pop(step["id"], None)
                state["history"].append({"task": step["id"], "status": "DEFERRED_LOCK_BUSY",
                                         "observed_at": now()})
                state["history"] = state["history"][-100:]
                write_atomic(state_path, state)
                return {"task": step["id"], "status": "DEFERRED_LOCK_BUSY",
                        "next_action": "RETRY_NEXT_CYCLE"}
            record["status"] = "REVIEW_REQUIRED"
            record["error"] = f"{type(exc).__name__}: {exc}"[-1000:]
        except Exception as exc:
            record["status"] = "REVIEW_REQUIRED"
            record["error"] = f"{type(exc).__name__}: {exc}"[-1000:]
        record["ended_at"] = now()
        record["duration_seconds"] = round(time.monotonic() - started, 3)
        state["history"].append({"task": step["id"], **record})
        state["history"] = state["history"][-100:]
        write_atomic(state_path, state)
        return {"task": step["id"], **record,
                "next_action": "REVIEW" if record["status"] != "PASS" else "RUN_NEXT_CYCLE"}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", action="store_true", help="Actually run one allowlisted operation")
    p.add_argument("--state-dir", type=Path, default=DEFAULT_ROOT)
    args = p.parse_args()
    print(json.dumps(execute(args.state_dir / "state.json", args.run), ensure_ascii=False, indent=2))

if __name__ == "__main__":
    main()
