#!/usr/bin/env python3
"""Combine independent HCC P0 Actions artifacts with a fail-closed CI gate.

PASS_EXECUTION means the probe ran, NOT that biological data, patient-level
RNA/protein pairs, prognostic endpoint validity or novel claims were verified.
Do not load original datasets here: only bounded aggregate per-track metadata.
"""
import argparse
import datetime as dt
import json
import pathlib

from run_parallel_p0_track import TRACKS, tasks_for

MAX_TRACK_STATUS_BYTES = 300_000
SUMMARY_SCHEMA = "hcc-p0-parallel-summary-v1"


def aggregate(artifact_root):
    artifact_root = pathlib.Path(artifact_root)
    rows = []
    failure_count = 0
    for name in TRACKS:
        candidates = [
            artifact_root / f"hcc-p0-{name}" / "track_status.json",
            artifact_root / name / "track_status.json",
        ]
        path = next((p for p in candidates if p.is_file()), None)
        expected_names = [task for task, _, _ in tasks_for(name, pathlib.Path("ignored"))]
        entry = {
            "track": name, "expected_tasks": expected_names,
            "status": "MISSING",
            "executed_task_count": 0, "successful_task_count": 0,
            "failed_tasks": expected_names.copy(),
        }
        if path is not None:
            try:
                if path.stat().st_size > MAX_TRACK_STATUS_BYTES:
                    raise ValueError("oversized status file")
                data = json.loads(path.read_text(encoding="utf-8"))
                items = data["tasks"]
                actual = [x["task"] for x in items]
                if (
                    data.get("schema_version") != "hcc-p0-parallel-v1"
                    or data.get("track") != name
                    or actual != expected_names
                    or data.get("task_count") != len(items)
                    or data.get("executed_task_count") != len(items)
                ):
                    raise ValueError("track status schema or task identity differs")
                successes = sum(x.get("result") == "PASS_EXECUTION" for x in items)
                if successes != data.get("successful_task_count"):
                    raise ValueError("successful task count mismatch")
                # Every script's successful execution must be backed by its
                # declared artifact if the lane task requires one.
                required = {t: expected is not None
                            for t, _, expected in tasks_for(name, pathlib.Path("ignored"))}
                valid = all(
                    x.get("result") in ("PASS_EXECUTION", "FAIL_EXECUTION")
                    and isinstance(x.get("exit_code"), int)
                    and (
                        x["result"] != "PASS_EXECUTION"
                        or (x["exit_code"] == 0 and
                            (not required[x["task"]] or x.get("artifact_present") is True))
                    )
                    for x in items
                )
                if not valid:
                    raise ValueError("inconsistent execution status or missing output")
                failed = [x["task"] for x in items if x["result"] != "PASS_EXECUTION"]
                entry.update({
                    "status": "PASS_EXECUTION" if not failed else "PARTIAL_EXECUTION_FAILURE",
                    "executed_task_count": len(items),
                    "successful_task_count": successes,
                    "failed_tasks": failed,
                })
            except (OSError, UnicodeDecodeError, ValueError, KeyError, TypeError):
                entry["status"] = "INVALID_STATUS_FILE"
        if entry["status"] != "PASS_EXECUTION":
            failure_count += 1
        rows.append(entry)

    total_expected = sum(len(x["expected_tasks"]) for x in rows)
    total_pass = sum(x["successful_task_count"] for x in rows)
    summary = {
        "schema_version": SUMMARY_SCHEMA,
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "track_count": len(rows),
        "tracks_passed_execution": len(rows) - failure_count,
        "tracks_failed_or_missing": failure_count,
        "tasks_expected": total_expected,
        "tasks_passed_execution": total_pass,
        "all_execution_tracks_passed": failure_count == 0,
        "V02_matched_rna_protein_quantitative_join": "BLOCKED",
        "V04_etiology_group_inference": "PARTIAL",
        "V12_prognostic_event_codebook": "PARTIAL",
        "clinical_event_counts_claimed": False,
        "biological_gate_pass_claimed": False,
        "notes": "CI PASS only checks execution and aggregate artifacts. Source discovery or parser availability never establishes clinical validity, RNA/protein patient matching, etiology cases or independent prognosis.",
        "tracks": rows,
    }
    return summary


def to_markdown(report):
    lines = [
        "# HCC parallel P0 execution audit", "",
        f"- Timestamp UTC: {report['generated_utc']}",
        f"- Tracks: {report['tracks_passed_execution']}/{report['track_count']} executed successfully",
        f"- Tasks: {report['tasks_passed_execution']}/{report['tasks_expected']} executed successfully",
        f"- Overall: {'PASS_EXECUTION' if report['all_execution_tracks_passed'] else 'FAIL_OR_INCOMPLETE_EXECUTION'}",
        "- Scientific evidence gates: V02 BLOCKED; V04 PARTIAL; V12 PARTIAL (unchanged)", "",
        "| Track | Execution | Successful tasks | Failures |",
        "| --- | --- | ---: | --- |"
    ]
    for row in report["tracks"]:
        failed = ", ".join(row["failed_tasks"]) if row["failed_tasks"] else "None"
        lines.append(
            f"| {row['track']} | {row['status']} | "
            f"{row['successful_task_count']}/{len(row['expected_tasks'])} | {failed} |"
        )
    lines.extend([
        "", "**Guardrail:** successful public metadata probes do not establish a",
        "patient-matched RNA/protein quantitative dataset, true HBV/HCV cases,",
        "validated RFS/OS events, or scientific significance.", ""
    ])
    return "\n".join(lines)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--artifact-root", required=True)
    p.add_argument("--out-dir", default="results/hcc/parallel/batch")
    args = p.parse_args()
    report = aggregate(args.artifact_root)
    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "batch_status.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (out / "batch_report.md").write_text(to_markdown(report), encoding="utf-8")
    print(to_markdown(report))
    if not report["all_execution_tracks_passed"]:
        raise SystemExit(1)

if __name__ == "__main__":
    main()
