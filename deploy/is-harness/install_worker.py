#!/usr/bin/env python3
"""Export or install a pinned queue bundle; never start or enable user services."""
from __future__ import annotations

import argparse
from contextlib import ExitStack, contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import stat
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
INSTALLER = "install_worker.py"
MANIFEST = "manifest.json"
BUNDLE_FILES = MODULES + UNITS + (INSTALLER,)
BUNDLE_SCHEMA = "is-harness-install-bundle-v1"
MAX_PAYLOAD_BYTES = 1024 * 1024
MAX_MANIFEST_BYTES = 16 * 1024
DEFAULT_EXECUTION_LOCK = Path(
    "/srv/is-analysis/results/is/github_actions_runs/is_research_execution.lock"
)
RUNTIME_UNITS = UNITS + ("is-auto-next.service", "is-auto-next.timer")
PROCESS_PATTERN = (
    r"[g]ithub_issue_queue_worker\.py|[h]arness_research_jobs\.py|"
    r"[h]arness_science_pipeline\.py|[h]arness_archived_evidence\.py|[a]uto_next\.py"
)
# Unknown or duplicate directives are rejected, including additional ExecStart*
# commands. Description text is harmless; every execution directive is pinned.
UNIT_CONTRACTS = {
    UNITS[0]: {
        ("Unit", "Description"): None,
        ("Unit", "Wants"): "network-online.target",
        ("Unit", "After"): "network-online.target",
        ("Service", "Type"): "oneshot",
        ("Service", "EnvironmentFile"): "%h/.config/is-harness/github.env",
        ("Service", "WorkingDirectory"): "%h/is-harness",
        ("Service", "ExecStart"):
            "/usr/bin/flock -n %t/is-harness-queue.lock /usr/bin/python3 "
            "%h/is-harness/github_issue_queue_worker.py --once",
        ("Service", "TimeoutStartSec"): "600",
        ("Service", "UMask"): "0077",
        ("Service", "NoNewPrivileges"): "yes",
        ("Service", "PrivateTmp"): "yes",
        ("Service", "ProtectSystem"): "full",
        ("Service", "ProtectHome"): "no",
    },
    UNITS[1]: {
        ("Unit", "Description"): None,
        ("Timer", "OnBootSec"): "30s",
        ("Timer", "OnUnitInactiveSec"): "60s",
        ("Timer", "AccuracySec"): "10s",
        ("Timer", "Unit"): UNITS[0],
        ("Install", "WantedBy"): "timers.target",
    },
}


def validate_revision(revision: str) -> None:
    if not isinstance(revision, str) or not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise ValueError("--revision must be a full lowercase 40-character commit SHA")


def git_file(repo: Path, revision: str, path: str) -> bytes:
    """Read only a bounded regular blob, never a Git symlink or working-tree file."""
    prefix = ["git", "-C", str(repo)]
    size = int(subprocess.run(
        [*prefix, "cat-file", "-s", f"{revision}:{path}"],
        check=True, capture_output=True, text=True, timeout=30,
    ).stdout.strip())
    if size < 1 or size > MAX_PAYLOAD_BYTES:
        raise ValueError(f"Pinned file size is outside the allowed bound: {path}")
    entry = subprocess.run(
        [*prefix, "ls-tree", revision, "--", path],
        check=True, capture_output=True, text=True, timeout=30,
    ).stdout.rstrip("\n").split("\t")
    if len(entry) != 2 or entry[1] != path:
        raise ValueError(f"Pinned file is missing or ambiguous: {path}")
    metadata = entry[0].split()
    if len(metadata) != 3 or metadata[0] not in {"100644", "100755"} or metadata[1] != "blob":
        raise ValueError(f"Pinned file must be regular, not a link: {path}")
    data = subprocess.run(
        [*prefix, "cat-file", "blob", metadata[2]],
        check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
    ).stdout
    if len(data) != size:
        raise ValueError(f"Pinned file changed size while reading: {path}")
    return data


def git_payloads(repo: Path, revision: str, *, include_installer: bool) -> dict[str, bytes]:
    validate_revision(revision)
    actual = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "--verify", f"{revision}^{{commit}}"],
        check=True, capture_output=True, text=True, timeout=30,
    ).stdout.strip()
    if actual != revision:
        raise ValueError("Revision did not resolve to the requested commit")
    payloads = {name: git_file(repo, revision, f"scripts/is/{name}") for name in MODULES}
    names = UNITS + ((INSTALLER,) if include_installer else ())
    payloads.update({name: git_file(repo, revision, f"deploy/is-harness/{name}") for name in names})
    return payloads


def validate_units(payloads: dict[str, bytes]) -> None:
    """Validate unit text without loading systemd or reading EnvironmentFile."""
    for name in UNITS:
        parsed = {}
        section = None
        sections = set()
        for raw_line in payloads[name].decode("utf-8").splitlines():
            line = raw_line.strip()
            if not line or line.startswith(("#", ";")):
                continue
            if line.startswith("[") and line.endswith("]"):
                section = line[1:-1]
                if section in sections:
                    raise ValueError(f"Duplicate unit section in {name}")
                sections.add(section)
                continue
            if section is None or "=" not in line or line.endswith("\\"):
                raise ValueError(f"Malformed unit directive in {name}")
            key, value = (part.strip() for part in line.split("=", 1))
            item = (section, key)
            if item in parsed or not value:
                raise ValueError(f"Duplicate or empty unit directive in {name}")
            parsed[item] = value
        contract = UNIT_CONTRACTS[name]
        if set(parsed) != set(contract):
            raise ValueError(f"Unexpected or missing unit directives in {name}")
        if sections != {section for section, _ in contract}:
            raise ValueError(f"Unexpected unit section in {name}")
        for item, expected in contract.items():
            if expected is not None and parsed[item] != expected:
                raise ValueError(f"Unsafe unit contract in {name}: {item[0]}.{item[1]}")


def _strict_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate manifest key: {key}")
        result[key] = value
    return result


def _read_regular(name, *, dir_fd=None, limit=MAX_PAYLOAD_BYTES,
                  allow_empty=False) -> tuple[bytes, int]:
    """Use a non-following descriptor so links and special files cannot be read."""
    descriptor = os.open(name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=dir_fd)
    with os.fdopen(descriptor, "rb") as stream:
        info = os.fstat(stream.fileno())
        if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
            raise ValueError(f"Bundle/install file must be regular and not linked: {name}")
        if info.st_size < (0 if allow_empty else 1) or info.st_size > limit:
            raise ValueError(f"File size is outside the allowed bound: {name}")
        data = stream.read(limit + 1)
        if len(data) != info.st_size or len(data) > limit:
            raise ValueError(f"File changed size while reading: {name}")
        return data, stat.S_IMODE(info.st_mode)


def read_bundle(bundle_dir: Path, revision: str) -> tuple[dict[str, bytes], bytes]:
    validate_revision(revision)
    descriptor = os.open(bundle_dir, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        expected = set(BUNDLE_FILES) | {MANIFEST}
        if set(os.listdir(descriptor)) != expected:
            raise ValueError("Bundle has missing, extra, or unexpected files")
        manifest_bytes, _ = _read_regular(MANIFEST, dir_fd=descriptor, limit=MAX_MANIFEST_BYTES)
        manifest = json.loads(manifest_bytes.decode("utf-8"), object_pairs_hook=_strict_object)
        if not isinstance(manifest, dict) or set(manifest) != {"schema", "source_commit", "files_sha256"}:
            raise ValueError("Malformed or incomplete bundle manifest")
        if manifest["schema"] != BUNDLE_SCHEMA:
            raise ValueError("Unsupported bundle manifest schema")
        if manifest["source_commit"] != revision:
            raise ValueError("Bundle source_commit does not match the expected revision")
        hashes = manifest["files_sha256"]
        if not isinstance(hashes, dict) or set(hashes) != set(BUNDLE_FILES):
            raise ValueError("Manifest must cover exactly the seven payload files")
        if any(not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value)
               for value in hashes.values()):
            raise ValueError("Manifest SHA256 values must be lowercase hexadecimal digests")
        payloads = {}
        for name in BUNDLE_FILES:
            data, _ = _read_regular(name, dir_fd=descriptor)
            if hashlib.sha256(data).hexdigest() != hashes[name]:
                raise ValueError(f"Bundle hash mismatch: {name}")
            payloads[name] = data
        if set(os.listdir(descriptor)) != expected:
            raise ValueError("Bundle directory changed during validation")
    finally:
        os.close(descriptor)
    validate_units(payloads)
    compile(payloads[INSTALLER], INSTALLER, "exec")
    return payloads, manifest_bytes


def export_bundle(repo: Path, revision: str, bundle_dir: Path) -> Path:
    """Export seven Git-pinned payloads and their manifest without importing them."""
    if bundle_dir.exists() or bundle_dir.is_symlink():
        raise ValueError("Export directory must be new and must not already exist")
    payloads = git_payloads(repo, revision, include_installer=True)
    validate_units(payloads)
    # Syntax checks execute no working-tree or pinned module code.
    for name in MODULES + (INSTALLER,):
        compile(payloads[name], name, "exec")
    manifest = {
        "schema": BUNDLE_SCHEMA,
        "source_commit": revision,
        "files_sha256": {name: hashlib.sha256(data).hexdigest() for name, data in payloads.items()},
    }
    payloads[MANIFEST] = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
    bundle_dir.mkdir(mode=0o700)  # Exclusive creation; do not create missing parent trees.
    created = []
    try:
        for name, data in payloads.items():
            path = bundle_dir / name
            descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
            created.append(path)
            with os.fdopen(descriptor, "wb") as stream:
                stream.write(data)
                stream.flush()
                os.fsync(stream.fileno())
    except BaseException:
        for path in reversed(created):
            path.unlink(missing_ok=True)
        bundle_dir.rmdir()
        raise
    return bundle_dir / MANIFEST


def validate_imports(payloads: dict[str, bytes]) -> None:
    with tempfile.TemporaryDirectory(prefix="is-harness-validate-") as directory:
        stage = Path(directory)
        for name in MODULES:
            (stage / name).write_bytes(payloads[name])
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
        subprocess.run([sys.executable, "-I", "-B", "-c", validation, str(stage), *MODULES],
                       check=True, capture_output=True, text=True, timeout=30)


def _reject_linked_parents(path: Path) -> None:
    for current in (path, *path.parents):
        try:
            info = current.lstat()
        except FileNotFoundError:
            continue
        if stat.S_ISLNK(info.st_mode):
            raise ValueError(f"Installation path must not be a symlink: {current}")
        if current != path and not stat.S_ISDIR(info.st_mode):
            raise ValueError(f"Installation parent must be a directory: {current}")


@contextmanager
def installation_lock(path: Path, *, create: bool, busy_message: str):
    _reject_linked_parents(path)
    flags = os.O_RDWR | os.O_NOFOLLOW | os.O_NONBLOCK
    if create:
        flags |= os.O_CREAT
    descriptor = os.open(path, flags, 0o600)
    with os.fdopen(descriptor, "rb") as stream:
        info = os.fstat(stream.fileno())
        if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
            raise ValueError(f"Execution lock must be a regular unlinked file: {path}")
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(busy_message) from error
        yield


def require_inactive_runtime() -> None:
    """Read only fixed user-unit properties and process counts while locks are held."""
    for unit in RUNTIME_UNITS:
        result = subprocess.run(
            ["systemctl", "--user", "show", unit, "--property=LoadState,ActiveState", "--no-pager"],
            capture_output=True, text=True, timeout=15,
        )
        fields = {}
        for line in result.stdout.splitlines():
            key, separator, value = line.partition("=")
            if not separator or key in fields:
                raise RuntimeError("Unable to verify inactive user services")
            fields[key] = value
        if (result.returncode not in {0, 4} or set(fields) != {"LoadState", "ActiveState"} or
                fields["LoadState"] not in {"loaded", "not-found", "masked"} or
                fields["ActiveState"] not in {"inactive", "failed"}):
            raise RuntimeError(f"Runtime unit is active or unverifiable: {unit}")
    result = subprocess.run(["pgrep", "-fc", PROCESS_PATTERN],
                            capture_output=True, text=True, timeout=15)
    if result.returncode != 1 or result.stdout.strip() != "0":
        raise RuntimeError("Queue/AUTO_NEXT worker processes are active or unverifiable")


def _stage_file(parent: Path, data: bytes, mode: int, *, prefix: str) -> Path:
    descriptor, temporary = tempfile.mkstemp(prefix=prefix, dir=parent)
    path = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.chmod(path, mode)
    except BaseException:
        path.unlink(missing_ok=True)
        raise
    return path


def replace_installation(files: dict[Path, tuple[bytes, int]]) -> None:
    """Stage all destinations first and restore every attempted replacement on error."""
    snapshots = {}
    for target in files:
        _reject_linked_parents(target)
        try:
            snapshots[target] = _read_regular(target, allow_empty=True)
        except FileNotFoundError:
            snapshots[target] = None
    created_dirs = []
    staged = {}
    backups = {}
    temporary_files = []
    attempted = []
    keep_backups = set()
    success = False
    try:
        for target, (data, mode) in files.items():
            missing = []
            parent = target.parent
            while not parent.exists():
                missing.append(parent)
                parent = parent.parent
            for directory in reversed(missing):
                directory.mkdir(mode=0o700)
                created_dirs.append(directory)
            _reject_linked_parents(target)
            staged[target] = _stage_file(target.parent, data, mode, prefix=".install-new-")
            temporary_files.append(staged[target])
            previous = snapshots[target]
            if previous is not None:
                backups[target] = _stage_file(target.parent, *previous, prefix=".install-rollback-")
                temporary_files.append(backups[target])
        for target in files:
            attempted.append(target)
            os.replace(staged[target], target)
        success = True
    except BaseException as error:
        rollback_errors = []
        for target in reversed(attempted):
            try:
                if snapshots[target] is None:
                    target.unlink(missing_ok=True)
                else:
                    os.replace(backups[target], target)
            except OSError as rollback_error:
                rollback_errors.append(f"{target}: {rollback_error}")
                if target in backups:
                    keep_backups.add(backups[target])
        if rollback_errors:
            raise RuntimeError("Installation rollback incomplete; protected backup files retained: "
                               + "; ".join(rollback_errors)) from error
        raise
    finally:
        for path in temporary_files:
            if path not in keep_backups:
                path.unlink(missing_ok=True)
        if not success:
            for directory in reversed(created_dirs):
                try:
                    directory.rmdir()
                except OSError:
                    pass


def _install(payloads: dict[str, bytes], revision: str, home: Path, runtime_dir: Path,
             *, execution_lock: Path | None, require_inactive: bool, receipt_extra=None) -> Path:
    validate_revision(revision)
    validate_units(payloads)
    if not runtime_dir.is_dir():
        raise ValueError("User runtime directory is unavailable; use a normal user login session")
    if require_inactive and execution_lock is None:
        raise ValueError("Inactive-runtime verification requires the configured shared execution lock")
    service_lock = runtime_dir / "is-harness-queue.lock"
    if execution_lock is not None:
        _reject_linked_parents(execution_lock)
        info = execution_lock.lstat()  # Never create the production result tree or its lock.
        if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
            raise ValueError("Configured shared execution lock must already be a regular file")
        if os.path.abspath(execution_lock) == os.path.abspath(service_lock):
            raise ValueError("Service and shared execution locks must be different")
    validate_imports(payloads)
    destination = home / "is-harness"
    units_dir = home / ".config/systemd/user"
    receipt = {
        "source_commit": revision,
        "files_sha256": {name: hashlib.sha256(payloads[name]).hexdigest() for name in MODULES + UNITS},
        "service_started": False,
        **(receipt_extra or {}),
    }
    receipt_path = destination / "installed_revision.json"
    files = {
        (destination / name if name in MODULES else units_dir / name):
            (payloads[name], 0o600 if name in MODULES else 0o644)
        for name in MODULES + UNITS
    }
    files[receipt_path] = ((json.dumps(receipt, indent=2, sort_keys=True) + "\n").encode(), 0o600)
    with ExitStack() as locks:
        locks.enter_context(installation_lock(
            service_lock, create=True,
            busy_message="Queue worker is active; installation held without changing files",
        ))
        if execution_lock is not None:
            locks.enter_context(installation_lock(
                execution_lock, create=False,
                busy_message="Shared research execution lock is busy; installation held without changing files",
            ))
        if require_inactive:
            require_inactive_runtime()
        replace_installation(files)
    return receipt_path


def install_worker(repo: Path, revision: str, home: Path, runtime_dir: Path,
                   *, execution_lock: Path | None = None, require_inactive: bool = False) -> Path:
    """Keep the existing local-Git API; CLI installations require both runtime locks."""
    payloads = git_payloads(repo, revision, include_installer=False)
    return _install(payloads, revision, home, runtime_dir,
                    execution_lock=execution_lock, require_inactive=require_inactive)


def install_bundle(bundle_dir: Path, revision: str, home: Path, runtime_dir: Path,
                   *, execution_lock: Path | None = None, require_inactive: bool = False) -> Path:
    payloads, manifest_bytes = read_bundle(bundle_dir, revision)
    return _install(payloads, revision, home, runtime_dir,
                    execution_lock=execution_lock, require_inactive=require_inactive,
                    receipt_extra={
                        "bundle_schema": BUNDLE_SCHEMA,
                        "bundle_manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
                        "installer_sha256": hashlib.sha256(payloads[INSTALLER]).hexdigest(),
                    })


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--revision", required=True, help="Reviewed, full Git commit SHA")
    sources = parser.add_mutually_exclusive_group()
    sources.add_argument("--repo-root", type=Path)
    sources.add_argument("--bundle-dir", type=Path)
    parser.add_argument("--export-bundle", type=Path, help="Create a new, flat pinned bundle directory")
    parser.add_argument("--execution-lock", type=Path, default=DEFAULT_EXECUTION_LOCK,
                        help="Existing lock shared by queue and AUTO_NEXT; never created by the installer")
    args = parser.parse_args(argv)
    if args.export_bundle is not None and args.bundle_dir is not None:
        parser.error("--export-bundle cannot be combined with --bundle-dir")
    repo = args.repo_root or Path(__file__).resolve().parents[2]
    try:
        if args.export_bundle is not None:
            manifest = export_bundle(repo, args.revision, args.export_bundle)
            print(f"Exported pinned queue bundle from {args.revision}; manifest={manifest}")
            return
        if os.geteuid() == 0:
            parser.error("Install as the dedicated non-root research account")
        configured_lock = os.environ.get("IS_HARNESS_EXECUTION_LOCK")
        if configured_lock and os.path.abspath(configured_lock) != os.path.abspath(args.execution_lock):
            raise ValueError("--execution-lock does not match IS_HARNESS_EXECUTION_LOCK")
        runtime_dir = Path(os.environ.get("XDG_RUNTIME_DIR", f"/run/user/{os.getuid()}"))
        operation = install_bundle if args.bundle_dir is not None else install_worker
        source = args.bundle_dir if args.bundle_dir is not None else repo
        receipt = operation(source, args.revision, Path.home(), runtime_dir,
                            execution_lock=args.execution_lock, require_inactive=True)
    except (OSError, ValueError, RuntimeError, SyntaxError, subprocess.SubprocessError) as error:
        parser.exit(1, f"Installation/export held: {error}\n")
    print(f"Installed validated queue bundle from {args.revision}; receipt={receipt}")
    print("No service or timer was started or enabled.")


if __name__ == "__main__":
    main()
