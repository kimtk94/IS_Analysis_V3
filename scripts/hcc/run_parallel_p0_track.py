#!/usr/bin/env python3
"""Run one independent HCC P0 audit lane in GitHub Actions.

Each task runs even if a previous task in the same lane fails.
Only aggregate public audit outputs and status metadata are uploaded.
GitHub Actions green = execution success, never biological gate PASS.
"""
import argparse
import datetime as dt
import json
import pathlib
import subprocess
import sys
import time

TRACKS = {
    "ng-public": [
        ("ng_gao_metadata", "probe_ng_gao_public.py", "--output"),
        ("quant_matrix_manifest", "probe_quant_matrix_manifest.py", "--output"),
        ("quant_candidate_index", "probe_quant_candidates.py", "--out"),
    ],
    "pdc-mapping": [
        ("pdc_quant_headers", "probe_pdc_quant_headers.py", "--output"),
        ("pdc_sample_mapping", "probe_pdc_sample_mapping.py", "--output"),
        ("pdc_biospecimen_lineage", "probe_pdc_biospecimen_lineage.py", "--output"),
        ("pdc_experimental_design", "probe_pdc_experimental_design.py", "--out"),
    ],
    "rna-access": [
        ("node_oep000321_access", "probe_node_oep000321.py", "--out"),
        ("hccdb25_access", "probe_hccdb25_access.py", "--out"),
        ("hccdb25_frontend_routes", "probe_hccdb25_assets.py", "--out"),
    ],
    "geo-clinical": [
        ("geo_metadata", "probe_public_metadata.py", "--output"),
    ],
    "xing-clinical": [
        ("xing_clinical_completeness", "probe_xing2023_clinical_completeness.py", "--out"),
        ("xing_clinical_complete_case", "probe_xing2023_clinical_complete_case.py", "--out"),
        ("xing_etiology_codebook_preflight", "probe_xing2023_etiology_codebook_preflight.py", "--out"),
    ],
    "regression": [
        ("feasibility_registry_report", "feasibility_audit.py", "--out-dir"),
    ],
}
ALL_TRACKS = tuple(TRACKS)
OUTPUT_SCOPE = "public metadata, source-checked aggregate-only, no patient identifiers"


def tasks_for(track, out_dir):
    if track not in TRACKS:
        raise ValueError("unknown HCC P0 lane")
    if track == "regression":
        return [
            ("all_hcc_synthetic_tests", [
                sys.executable, "-m", "unittest", "discover", "-s", "tests",
                "-p", "test_hcc_*.py", "-q"
            ], None),
            ("feasibility_registry_report", [
                sys.executable, "scripts/hcc/feasibility_audit.py",
                "--out-dir", str(out_dir)
            ], None),
        ]
    return [
        (name, [sys.executable, f"scripts/hcc/{script}", option,
                str(out_dir / f"{name}.json")], out_dir / f"{name}.json")
        for name, script, option in TRACKS[track]
    ]


def run_track(track, out_dir, runner=subprocess.run):
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    entries = []
    for task, command, expected in tasks_for(track, out_dir):
        started = time.monotonic()
        try:
            # The called scripts own their logs, but their stdout/stderr can
            # include raw response data. Keep output private to this process.
            run = runner(command, capture_output=True, text=True, timeout=570,
                         check=False)
            code = int(run.returncode)
            failure = "nonzero_exit" if code != 0 else None
        except subprocess.TimeoutExpired:
            code, failure = 124, "time_limit_exceeded"
        except OSError:
            code, failure = 127, "runner_or_source_unavailable"
        if code == 0 and expected is not None:
            if not expected.is_file():
                code, failure = 65, "missing_audit_output"
            else:
                try:
                    parsed = json.loads(expected.read_text(encoding="utf-8"))
                    if not isinstance(parsed, dict):
                        raise ValueError("audit JSON must be object")
                except (OSError, ValueError):
                    code, failure = 65, "malformed_audit_json"
        entry = {
            "task": task, "result": "PASS_EXECUTION" if code == 0 else "FAIL_EXECUTION",
            "exit_code": code, "error_class": failure,
            "seconds_elapsed": round(time.monotonic() - started, 2),
            "artifact_present": expected is not None and expected.is_file(),
        }
        entries.append(entry)
        print(f"{track}/{task}: {entry['result']} ({entry['seconds_elapsed']}s)",
              flush=True)
    all_pass = all(x["result"] == "PASS_EXECUTION" for x in entries)
    summary = {
        "schema_version": "hcc-p0-parallel-v1",
        "track": track,
        "checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": "PASS_EXECUTION" if all_pass else "FAIL_EXECUTION",
        "tasks": entries,
        "task_count": len(entries),
        "executed_task_count": len(entries),
        "successful_task_count": sum(x["result"] == "PASS_EXECUTION" for x in entries),
        "biological_claim_supported": False,
        "rna_protein_same_patient_quant_join_verified": False,
        "candidate_recurrence_event_codebook_verified": False,
        "interpretation": ("CI success checks accessibility/parsing/metadata and "
                           "data integrity ONLY; cannot upgrade V02/V04/V12 to "
                           "biological PASS or infer new paired sample counts."),
        "artifact_policy": OUTPUT_SCOPE,
    }
    (out_dir / "track_status.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8"
    )
    return summary


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--track", choices=ALL_TRACKS, required=True)
    p.add_argument("--out-dir")
    args = p.parse_args()
    target = args.out_dir or f"results/hcc/parallel/{args.track}"
    result = run_track(args.track, target)
    if result["status"] != "PASS_EXECUTION":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
