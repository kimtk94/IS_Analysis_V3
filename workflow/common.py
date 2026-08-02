"""Shared, dependency-free provenance helpers for workflow entry points."""
from __future__ import annotations
import argparse, hashlib, json
from datetime import datetime, timezone
from pathlib import Path

def existing_file(value: str) -> Path:
    path = Path(value)
    if not path.is_file(): raise argparse.ArgumentTypeError(f"file not found: {value}")
    return path

def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""): h.update(block)
    return h.hexdigest()

def checkpoint(stage: str, output: Path, inputs: list[Path], parameters: dict) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {"stage": stage, "created_at": datetime.now(timezone.utc).isoformat(),
               "inputs": [{"path": str(p), "sha256": digest(p)} for p in inputs],
               "parameters": parameters, "status": "validated"}
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
