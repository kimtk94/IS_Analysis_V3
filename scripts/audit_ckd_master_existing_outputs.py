#!/usr/bin/env python3
"""Audit existing CKD outputs against the MASTER Stage 0-19 map.

Read-only. Verifies mapped artifacts and separates conceptual completion from
mere file existence.
"""
from __future__ import annotations
import argparse, csv, gzip, json, glob
from pathlib import Path

def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as f:
        yield from csv.DictReader(f, delimiter="\t")

def open_text(path):
    p = Path(path)
    return gzip.open(p, "rt", encoding="utf-8", errors="replace", newline="") if p.suffix == ".gz" else p.open("r", encoding="utf-8", errors="replace", newline="")

def header(path):
    p = Path(path)
    if p.suffix.lower() not in {".tsv", ".gz", ".txt", ".csv"}:
        return ""
    try:
        with open_text(p) as f:
            line = f.readline().rstrip("\r\n")
        return line[:4000]
    except Exception:
        return ""

def resolve(pattern):
    hits = sorted(glob.glob(pattern))
    if hits:
        return [Path(x) for x in hits]
    p = Path(pattern)
    return [p] if p.exists() else []

def next_action(status, exists):
    if status == "EXISTING" and exists:
        return "bridge_or_reuse"
    if status == "PARTIAL" and exists:
        return "materialize_generic_stage_output"
    if status == "PLANNED":
        return "run_or_acquire_required_inputs"
    if not exists:
        return "locate_or_regenerate_expected_artifact"
    return "review"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True)
    a = ap.parse_args()

    rows = []
    for m in read_tsv(a.map):
        hits = resolve(m["path_pattern"])
        chosen = hits[0] if hits else None
        row = dict(m)
        row.update({
            "artifact_exists": int(bool(hits)),
            "resolved_path": str(chosen) if chosen else "",
            "match_count": len(hits),
            "size_bytes": chosen.stat().st_size if chosen and chosen.is_file() else "",
            "header": header(chosen) if chosen and chosen.is_file() else "",
            "next_action": next_action(m["current_status"], bool(hits)),
        })
        rows.append(row)

    fields = list(rows[0]) if rows else []
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with a.output.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    summary = {
        "stage_count": len(rows),
        "artifacts_found": sum(r["artifact_exists"] for r in rows),
        "existing_stage_artifacts_found": [
            int(r["stage"]) for r in rows
            if r["current_status"] == "EXISTING" and r["artifact_exists"]
        ],
        "partial_stage_artifacts_found": [
            int(r["stage"]) for r in rows
            if r["current_status"] == "PARTIAL" and r["artifact_exists"]
        ],
        "planned_stages": [
            int(r["stage"]) for r in rows if r["current_status"] == "PLANNED"
        ],
        "missing_expected_existing": [
            int(r["stage"]) for r in rows
            if r["current_status"] == "EXISTING" and not r["artifact_exists"]
        ],
        "ready_for_bridge": [
            int(r["stage"]) for r in rows
            if r["artifact_exists"] and r["current_status"] in {"EXISTING","PARTIAL"}
        ],
    }
    a.summary.parent.mkdir(parents=True, exist_ok=True)
    a.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("CKD_MASTER_EXISTING_OUTPUT_AUDIT_PASS")

if __name__ == "__main__":
    main()
