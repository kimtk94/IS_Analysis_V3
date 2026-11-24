#!/usr/bin/env python3
"""Bounded queue diagnostics: metadata only, no installation or worker execution."""
from __future__ import annotations

import argparse
import datetime as dt
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import selectors
import stat
import subprocess
import time

MODULES = (
    "github_issue_queue_worker.py", "harness_research_jobs.py",
    "harness_science_pipeline.py", "harness_archived_evidence.py",
)
UNITS = ("is-harness-queue.service", "is-harness-queue.timer")
PROCESS_NAMES = {
    "github_issue_queue_worker.py": "queue_worker",
    "harness_research_jobs.py": "research_helper",
    "auto_next.py": "auto_next",
}
MAX_FILE_BYTES = 1024 * 1024
MAX_COMMAND_BYTES = 256 * 1024
MAX_PROCESSES = 4096
EXPECTED_SERVICE_TIMEOUT_SECONDS = 600
RUNNING_SERVICE_STATES = {
    ("activating", "start-pre"), ("activating", "start"),
    ("activating", "start-post"),
}
ENUMS = {
    "LoadState": {"loaded", "not-found", "error", "bad-setting", "masked", "stub", "merged"},
    "ActiveState": {"active", "inactive", "failed", "activating", "deactivating", "reloading", "maintenance", "refreshing"},
    "SubState": {"dead", "running", "exited", "failed", "waiting", "elapsed", "start-pre", "start", "start-post", "stop", "stop-post", "stop-sigterm", "stop-sigkill", "auto-restart", "condition", "cleaning"},
    "Result": {"success", "resources", "timeout", "exit-code", "signal", "core-dump", "watchdog", "start-limit-hit", "protocol", "oom-kill", "condition", "assert"},
    "UnitFileState": {"enabled", "enabled-runtime", "disabled", "static", "indirect", "masked", "masked-runtime", "linked", "linked-runtime", "alias", "generated", "transient", "bad"},
}
JOURNAL_CATEGORIES = {
    "stale_request_hold": r"HOLD_STALE_REQUEST",
    "unknown_request_age": r"HOLD_REQUEST_TIME_UNKNOWN",
    "inventory_hold": r"HOLD_QUEUE_INVENTORY_LIMIT",
    "uncertain_queue_state": r"HOLD_QUEUE_STATE_UNCERTAIN",
    "token_missing": r"IS_QUEUE_GITHUB_TOKEN missing",
    "authorization_error": r"HTTP(?:Error)?[^\n]{0,60}\b(?:401|403)\b|unauthorized|forbidden",
    "network_error": r"URLError|TimeoutError|timed out|Connection refused|Name or service not known",
    "missing_dependency": r"helper missing|module (?:missing|not installed)|ModuleNotFoundError|FileNotFoundError",
    "execution_failed": r"issue #\d+: FAIL|Research audit failed|Command failed",
    "shared_execution_busy": r"SHARED_EXECUTION_LOCK_BUSY",
}


def read_file(path: Path, limit: int = MAX_FILE_BYTES) -> dict:
    """Read a regular, non-symlink file without creating or modifying it."""
    descriptor = None
    try:
        descriptor = os.open(path, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW | os.O_NONBLOCK)
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            return {"status": "not_regular"}
        if metadata.st_size > limit:
            return {"status": "oversized"}
        with os.fdopen(descriptor, "rb") as stream:
            descriptor = None
            data = stream.read(limit + 1)
        if len(data) > limit:
            return {"status": "oversized"}
        return {"status": "ok", "data": data, "mtime": metadata.st_mtime}
    except FileNotFoundError:
        return {"status": "missing"}
    except OSError as error:
        return {"status": "symlink" if error.errno == errno.ELOOP else "unreadable"}
    finally:
        if descriptor is not None:
            os.close(descriptor)


def read_json(path: Path) -> dict:
    result = read_file(path)
    if result["status"] != "ok":
        return result
    try:
        value = json.loads(result.pop("data"))
        if not isinstance(value, dict):
            return {"status": "invalid_json_shape"}
        result["value"] = value
        return result
    except (ValueError, UnicodeError, RecursionError):
        return {"status": "invalid_json"}


def command(arguments: list[str], timeout: float = 5, limit: int = MAX_COMMAND_BYTES) -> dict:
    """Bound fixed read-only command output and discard stderr, never echo errors."""
    environment = dict(os.environ, LC_ALL="C", SYSTEMD_PAGER="cat", SYSTEMD_COLORS="0",
                       GIT_OPTIONAL_LOCKS="0", GIT_TERMINAL_PROMPT="0")
    try:
        process = subprocess.Popen(arguments, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                   stderr=subprocess.DEVNULL, env=environment)
    except OSError:
        return {"status": "unavailable", "returncode": None}
    chunks = bytearray()
    deadline = time.monotonic() + timeout
    outcome = "ok"
    try:
        with selectors.DefaultSelector() as selector:
            selector.register(process.stdout, selectors.EVENT_READ)
            while True:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    outcome = "timeout"
                    break
                ready = selector.select(remaining)
                if not ready:
                    outcome = "timeout"
                    break
                data = os.read(process.stdout.fileno(), min(65536, limit - len(chunks) + 1))
                if not data:
                    break
                chunks.extend(data)
                if len(chunks) > limit:
                    outcome = "oversized"
                    break
        if outcome != "ok":
            process.kill()
        try:
            code = process.wait(timeout=max(0.01, deadline - time.monotonic()))
        except subprocess.TimeoutExpired:
            outcome = "timeout"
            process.kill()
            code = process.wait()
    finally:
        process.stdout.close()
        if process.poll() is None:
            process.kill()
            process.wait()
    if outcome != "ok":
        return {"status": outcome, "returncode": code}
    if code != 0:
        return {"status": "command_failed", "returncode": code}
    return {"status": "ok", "returncode": code,
            "output": bytes(chunks).decode("utf-8", errors="replace")}


def installation(install_root: Path, unit_root: Path) -> dict:
    receipt = read_json(install_root / "installed_revision.json")
    value = receipt.get("value", {})
    revision = value.get("source_commit")
    valid_revision = isinstance(revision, str) and bool(re.fullmatch(r"[0-9a-f]{40}", revision))
    expected = value.get("files_sha256", {})
    if not isinstance(expected, dict):
        expected = {}
    files = {}
    for name in MODULES + UNITS:
        source = read_file((install_root if name in MODULES else unit_root) / name)
        digest = hashlib.sha256(source["data"]).hexdigest() if source["status"] == "ok" else None
        expected_hash = expected.get(name)
        if not isinstance(expected_hash, str) or not re.fullmatch(r"[0-9a-f]{64}", expected_hash):
            expected_hash = None
        files[name] = {"status": source["status"], "sha256": digest,
                       "receipt_sha256": expected_hash,
                       "matches_receipt": bool(digest and expected_hash and digest == expected_hash)}
    verified = receipt["status"] == "ok" and valid_revision and all(row["matches_receipt"] for row in files.values())
    return {"status": "verified" if verified else "unverified", "receipt_status": receipt["status"],
            "receipt_commit": revision if valid_revision else None, "files": files}


def repository_head(repo_root: Path) -> dict:
    result = command(["git", "--no-optional-locks", "-C", str(repo_root), "rev-parse", "--verify", "HEAD"])
    head = result.get("output", "").strip()
    valid = result["status"] == "ok" and bool(re.fullmatch(r"[0-9a-f]{40}", head))
    return {"status": "ok" if valid else "unknown", "head": head if valid else None,
            "query_status": result["status"]}


def nonnegative_integer(value):
    if isinstance(value, str) and re.fullmatch(r"\d{1,20}", value):
        return int(value)
    return None


def duration_seconds(value):
    if not isinstance(value, str):
        return None
    if re.fullmatch(r"\d{1,20}", value):
        return int(value) / 1000000
    units = {"us": 0.000001, "ms": 0.001, "s": 1, "min": 60, "h": 3600, "d": 86400}
    parts = re.findall(r"(\d+(?:\.\d+)?)(us|ms|min|s|h|d)", value)
    if not parts or "".join(number + unit for number, unit in parts) != value.replace(" ", ""):
        return None
    return sum(float(number) * units[unit] for number, unit in parts)


def systemd_unit(unit: str, timer: bool = False) -> dict:
    names = ["LoadState", "ActiveState", "SubState"]
    names += (["UnitFileState", "NextElapseUSecMonotonic", "LastTriggerUSecMonotonic"] if timer else
              ["Result", "ExecMainStatus", "ExecMainStartTimestampMonotonic", "ExecMainExitTimestampMonotonic", "TimeoutStartUSec"])
    result = command(["systemctl", "--user", "show", unit, "--no-pager", "--property=" + ",".join(names)])
    if result["status"] != "ok":
        return {"status": "unknown", "query_status": result["status"], "properties": {}}
    raw = dict(line.split("=", 1) for line in result.get("output", "").splitlines() if "=" in line)
    properties = {name: raw.get(name) if raw.get(name) in ENUMS[name] else "unknown" for name in names if name in ENUMS}
    for name in names:
        if name not in ENUMS:
            if name == "TimeoutStartUSec":
                properties["TimeoutStartSeconds"] = duration_seconds(raw.get(name))
            else:
                properties[name] = nonnegative_integer(raw.get(name))
    required = ["LoadState", "ActiveState", "SubState"] + (["UnitFileState"] if timer else ["Result"])
    unknown = any(properties[name] == "unknown" for name in required)
    if not timer:
        unknown |= (properties["ExecMainStatus"] is None or properties["ExecMainStartTimestampMonotonic"] is None
                    or properties["TimeoutStartSeconds"] is None)
    return {"status": "unknown" if unknown else "observed", "query_status": "ok", "properties": properties}


def existing_lock(path: Path) -> dict:
    descriptor = None
    try:
        descriptor = os.open(path, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW | os.O_NONBLOCK)
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            return {"status": "unknown", "reason": "not_regular"}
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return {"status": "held"}
        fcntl.flock(descriptor, fcntl.LOCK_UN)
        return {"status": "free"}
    except FileNotFoundError:
        return {"status": "missing"}
    except OSError:
        return {"status": "unknown", "reason": "unreadable"}
    finally:
        if descriptor is not None:
            os.close(descriptor)


def process_counts(proc_root: Path) -> dict:
    counts = dict.fromkeys(PROCESS_NAMES.values(), 0)
    inspected = unreadable = 0
    try:
        with os.scandir(proc_root) as entries:
            for entry in entries:
                if not entry.name.isdigit():
                    continue
                inspected += 1
                if inspected > MAX_PROCESSES:
                    return {"status": "truncated", "counts": counts}
                try:
                    if entry.stat(follow_symlinks=False).st_uid != os.getuid():
                        continue
                except FileNotFoundError:
                    continue
                result = read_file(Path(entry.path) / "cmdline", 8192)
                if result["status"] == "missing":
                    continue  # A process may exit during this bounded snapshot.
                if result["status"] != "ok":
                    unreadable += 1
                    continue
                args = result["data"].split(b"\0")
                names = [os.path.basename(arg.decode("utf-8", "replace")) for arg in args if arg]
                if not names:
                    continue
                if not (re.fullmatch(r"python(?:\d+(?:\.\d+)*)?", names[0]) or names[0] in PROCESS_NAMES):
                    continue
                for name, category in PROCESS_NAMES.items():
                    if name in names:
                        counts[category] += 1
    except OSError:
        return {"status": "unknown", "counts": counts}
    return {"status": "partial" if unreadable else "observed", "counts": counts,
            "unreadable_processes": unreadable, "snapshot_only": True}


def timestamp(value):
    if not isinstance(value, str) or len(value) > 50:
        return None
    try:
        parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
        return parsed.astimezone(dt.timezone.utc) if parsed.tzinfo is not None else None
    except ValueError:
        return None


def evidence_freshness(path: Path, schema: str, now: dt.datetime, max_age: int, heartbeat=False) -> dict:
    result = read_json(path)
    if result["status"] != "ok":
        return {"status": result["status"], "age_seconds": None}
    value = result["value"]
    observed = timestamp(value.get("observed_at" if heartbeat else "generated_at"))
    elapsed = (now - observed).total_seconds() if observed else None
    age = int(elapsed) if elapsed is not None else None
    valid = value.get("schema") == schema and observed is not None and elapsed >= 0
    summary = {"status": "fresh" if valid and age <= max_age else "stale" if valid else "unknown",
               "age_seconds": age if valid else None,
               "recorded_at": observed.isoformat() if observed else None}
    if heartbeat:
        state = value.get("watcher_status", "")
        if not isinstance(state, str):
            state = ""
        exact = {"UNCHANGED_SOURCES_FAST_CHECK", "PENDING_REFRESH_EXISTS", "SLOT_ALREADY_ENQUEUED"}
        prefixes = ("HOLD_STALE_REQUEST", "HOLD_REQUEST_TIME_UNKNOWN", "HOLD_QUEUE_INVENTORY_LIMIT")
        summary["queue_decision"] = (state if state in exact else
            "FALLBACK_ENQUEUED" if re.fullmatch(r"FALLBACK_ENQUEUED_ISSUE_[1-9]\d*", state) else
            next((prefix for prefix in prefixes if state.startswith(prefix)), "unknown"))
    else:
        integrity = value.get("integrity")
        summary["integrity"] = integrity if isinstance(integrity, str) and integrity in {"PASS", "FAIL"} else "unknown"
    return summary


def journal_summary(line_limit: int = 100) -> dict:
    result = command(["journalctl", "--user", "-u", "is-harness-queue.service", "--no-pager",
                      "--since=-24h", "-n", str(line_limit), "--output=json",
                      "--output-fields=MESSAGE,PRIORITY,__REALTIME_TIMESTAMP"])
    if result["status"] != "ok":
        return {"status": "unknown", "query_status": result["status"], "error_counts": {}}
    counts = dict.fromkeys(JOURNAL_CATEGORIES, 0)
    counts["other_error"] = 0
    rows = invalid = successful = 0
    for line in result.get("output", "").splitlines()[:line_limit]:
        try:
            row = json.loads(line)
            message = row.get("MESSAGE")
            if not isinstance(message, str):
                raise ValueError("Not a textual journal entry")
        except (ValueError, AttributeError):
            invalid += 1
            continue
        rows += 1
        matched = False
        for category, pattern in JOURNAL_CATEGORIES.items():
            if re.search(pattern, message, re.IGNORECASE):
                counts[category] += 1
                matched = True
        if not matched and (str(row.get("PRIORITY")) in {"0", "1", "2", "3"} or "queue poll error:" in message):
            counts["other_error"] += 1
        if re.search(r"issue #\d+: PASS\b", message):
            successful += 1
    return {"status": "unknown" if invalid or not rows else "observed", "query_status": "ok",
            "entries_examined": rows, "invalid_entries": invalid, "window_hours": 24,
            "line_limit": line_limit, "completed_pass_count": successful, "error_counts": counts}


def audit(args, now=None) -> dict:
    now = now or dt.datetime.now(dt.timezone.utc)
    installed = installation(args.install_root, args.unit_root)
    repository = repository_head(args.repo_root)
    service = systemd_unit(UNITS[0])
    timer = systemd_unit(UNITS[1], timer=True)
    locks = {"queue": existing_lock(args.runtime_dir / "is-harness-queue.lock"),
             "research_execution": existing_lock(args.results_root / "is_research_execution.lock")}
    processes = process_counts(args.proc_root)
    heartbeat = evidence_freshness(args.results_root / "evidence/harness_heartbeat.json",
                                   "is-harness-heartbeat-v1", now, args.freshness_seconds, heartbeat=True)
    results = evidence_freshness(args.results_root / "evidence/research_manifest.json",
                                 "is-research-evidence-v1", now, args.freshness_seconds)
    journal = journal_summary(args.journal_lines)
    problems = []
    if os.geteuid() == 0:
        problems.append("UNEXPECTED_ROOT_ACCOUNT")
    if installed["status"] != "verified":
        problems.append("INSTALLATION_UNVERIFIED")
    if repository["status"] != "ok":
        problems.append("REPOSITORY_HEAD_UNKNOWN")
    elif installed["receipt_commit"] and repository["head"] != installed["receipt_commit"]:
        problems.append("REPOSITORY_INSTALL_REVISION_MISMATCH")
    state = service["properties"]
    service_state = (state.get("ActiveState"), state.get("SubState"))
    service_idle = service_state == ("inactive", "dead")
    service_running = service_state in RUNNING_SERVICE_STATES
    if service["status"] != "observed":
        problems.append("SERVICE_STATE_UNKNOWN")
    elif state["LoadState"] != "loaded" or state["ActiveState"] == "failed" or state["Result"] != "success" or state["ExecMainStatus"] != 0:
        problems.append("SERVICE_BLOCKED_OR_FAILED")
    elif not state["ExecMainStartTimestampMonotonic"]:
        problems.append("SERVICE_EXECUTION_NOT_OBSERVED")
    if service["status"] == "observed":
        if state["TimeoutStartSeconds"] != EXPECTED_SERVICE_TIMEOUT_SECONDS:
            problems.append("SERVICE_TIMEOUT_MISMATCH")
        if not (service_idle or service_running):
            problems.append("SERVICE_STATE_UNSUPPORTED")
    timer_state = timer["properties"]
    if timer["status"] != "observed":
        problems.append("TIMER_STATE_UNKNOWN")
    elif (timer_state["LoadState"] != "loaded" or timer_state["ActiveState"] != "active" or
          timer_state["UnitFileState"] not in {"enabled", "enabled-runtime"}):
        problems.append("TIMER_NOT_ACTIVE_AND_ENABLED")
    elif (timer_state["SubState"] not in {"waiting", "running"} or
          (timer_state["SubState"] == "running" and not service_running)):
        problems.append("TIMER_STATE_INCOMPATIBLE")
    elif timer_state["NextElapseUSecMonotonic"] is None:
        problems.append("TIMER_NEXT_ACTIVATION_UNKNOWN")
    elif service_idle and timer_state["NextElapseUSecMonotonic"] == 0:
        problems.append("TIMER_NO_PENDING_ACTIVATION")
    if any(item["status"] not in {"held", "free"} for item in locks.values()):
        problems.append("LOCK_STATE_UNKNOWN_OR_MISSING")
    if processes["status"] != "observed":
        problems.append("PROCESS_SNAPSHOT_INCOMPLETE")
    if processes["counts"]["queue_worker"] > 1:
        problems.append("MULTIPLE_QUEUE_WORKERS")
    if heartbeat["status"] != "fresh" or heartbeat.get("queue_decision") == "unknown":
        problems.append("HEARTBEAT_STALE_OR_UNKNOWN")
    elif heartbeat.get("queue_decision", "").startswith("HOLD_"):
        problems.append("RECORDED_QUEUE_HOLD")
    if results["status"] != "fresh" or results.get("integrity") != "PASS":
        problems.append("RESULT_METADATA_STALE_OR_UNVERIFIED")
    if journal["status"] != "observed":
        problems.append("JOURNAL_UNAVAILABLE_OR_EMPTY")
    elif any(value for name, value in journal["error_counts"].items() if name != "shared_execution_busy"):
        problems.append("ERRORS_IN_BOUNDED_JOURNAL_WINDOW")
    complete = (
        repository["status"] == "ok" and service["status"] == "observed" and timer["status"] == "observed"
        and processes["status"] == "observed" and journal["query_status"] == "ok"
        and not journal.get("invalid_entries", 0)
        and installed["receipt_status"] in {"ok", "missing"}
        and all(item["status"] in {"ok", "missing"} for item in installed["files"].values())
        and all(item["status"] in {"held", "free", "missing"} for item in locks.values())
        and heartbeat["status"] in {"fresh", "stale", "missing"}
        and results["status"] in {"fresh", "stale", "missing"}
    )
    return {"schema": "is-queue-runtime-audit-v1", "observed_at": now.isoformat(), "read_only": True,
            "collection_status": "COMPLETE" if complete else "PARTIAL",
            "status": "OPERATIONAL_OK" if not problems else "ATTENTION_REQUIRED",
            "scientific_validity": "NOT_EVALUATED", "uid": os.getuid(), "problems": problems,
            "configuration_scope": "CLI paths only; token environment is not read",
            "recovery_authorized": False,
            "installation": installed, "repository": repository, "service": service, "timer": timer,
            "locks": locks, "processes": processes, "heartbeat": heartbeat, "result_metadata": results,
            "journal": journal}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("/srv/is-analysis/IS_Analysis_V3"))
    parser.add_argument("--install-root", type=Path, default=Path.home() / "is-harness")
    parser.add_argument("--unit-root", type=Path, default=Path.home() / ".config/systemd/user")
    parser.add_argument("--runtime-dir", type=Path, default=Path(os.environ.get("XDG_RUNTIME_DIR", f"/run/user/{os.getuid()}")))
    parser.add_argument("--results-root", type=Path, default=Path("/srv/is-analysis/results/is/github_actions_runs"))
    parser.add_argument("--proc-root", type=Path, default=Path("/proc"))
    parser.add_argument("--freshness-seconds", type=int, default=3600)
    parser.add_argument("--journal-lines", type=int, default=100)
    args = parser.parse_args(argv)
    if not 1 <= args.journal_lines <= 200 or not 60 <= args.freshness_seconds <= 86400:
        parser.error("journal-lines must be 1..200 and freshness-seconds must be 60..86400")
    try:
        report = audit(args)
    except Exception:
        # Unknown state must fail visibly without disclosing exception payloads.
        report = {"schema": "is-queue-runtime-audit-v1", "read_only": True,
                  "collection_status": "PARTIAL",
                  "status": "ATTENTION_REQUIRED", "problems": ["DIAGNOSTIC_INTERNAL_ERROR"],
                  "scientific_validity": "NOT_EVALUATED"}
    print(json.dumps(report, sort_keys=True, indent=2))
    return 0 if report["status"] == "OPERATIONAL_OK" else 2


if __name__ == "__main__":
    raise SystemExit(main())
