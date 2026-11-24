#!/usr/bin/env python3
"""Install the queue bundle from one local, immutable Git commit; never start it."""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


MODULES = (
    "github_issue_queue_worker.py",
    "harness_research_jobs.py",
    "harness_science_pipeline.py",
    "harness_archived_evidence.py",
)
UNITS = ("is-harness-queue.service", "is-harness-queue.timer")


def git_file(repo: Path, revision: str, path: str) -> bytes:
    return subprocess.run(
        ["git", "-C", str(repo), "show", f"{revision}:{path}"],
        check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
    ).stdout


def atomic_write(path: Path, data: bytes, mode: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".install-", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.chmod(temporary, mode)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def install_worker(repo: Path, revision: str, home: Path, runtime_dir: Path) -> Path:
    """Validate the entire bundle before replacing files under the service lock."""
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise ValueError("--revision must be a full lowercase 40-character commit SHA")
    actual = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "--verify", f"{revision}^{{commit}}"],
        check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=30,
    ).stdout.strip()
    if actual != revision:
        raise ValueError("Revision did not resolve to the requested commit")
    if not runtime_dir.is_dir():
        raise ValueError("User runtime directory is unavailable; use a normal user login session")

    bundle = {name: git_file(repo, revision, f"scripts/is/{name}") for name in MODULES}
    bundle.update({name: git_file(repo, revision, f"deploy/is-harness/{name}") for name in UNITS})
    with tempfile.TemporaryDirectory(prefix="is-harness-install-") as directory:
        stage = Path(directory)
        for name in MODULES:
            (stage / name).write_bytes(bundle[name])
        # Imports verify companion modules as well as Python syntax. No CLI main
        # is invoked and no research operation, API request, or service is started.
        validation = """
import importlib.util, pathlib, sys
base = pathlib.Path(sys.argv[1])
sys.path.insert(0, str(base))
for name in sys.argv[2:]:
    spec = importlib.util.spec_from_file_location(pathlib.Path(name).stem, base / name)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
"""
        subprocess.run([sys.executable, "-I", "-c", validation, str(stage), *MODULES],
                       check=True, capture_output=True, text=True, timeout=30)

    destination = home / "is-harness"
    units_dir = home / ".config/systemd/user"
    with (runtime_dir / "is-harness-queue.lock").open("a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError("Queue worker is active; installation held without changing files") from error
        for name, data in bundle.items():
            target = destination / name if name in MODULES else units_dir / name
            atomic_write(target, data, 0o600 if name in MODULES else 0o644)
        receipt = {
            "source_commit": revision,
            "files_sha256": {name: hashlib.sha256(data).hexdigest() for name, data in bundle.items()},
            "service_started": False,
        }
        receipt_path = destination / "installed_revision.json"
        atomic_write(receipt_path, (json.dumps(receipt, indent=2) + "\n").encode(), 0o600)
    return receipt_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--revision", required=True, help="Reviewed, full Git commit SHA")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    if os.geteuid() == 0:
        parser.error("Install as the dedicated non-root research account")
    runtime_dir = Path(os.environ.get("XDG_RUNTIME_DIR", f"/run/user/{os.getuid()}"))
    try:
        receipt = install_worker(args.repo_root.resolve(), args.revision, Path.home(), runtime_dir)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        parser.exit(1, f"Installation held: {error}\n")
    print(f"Installed validated queue bundle from {args.revision}; receipt={receipt}")
    print("Service and timer remain unchanged. Review configuration before activation.")


if __name__ == "__main__":
    main()
