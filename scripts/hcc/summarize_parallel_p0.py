#!/usr/bin/env python3
"""Summarize multi-lane HCC P0 Actions artifacts; never promote scientific gates."""
import argparse
import datetime as dt
import json
import os
import pathlib

from run_parallel_p0_track import ALL_TRACKS


def summarize(root, expected=ALL_TRACKS):
    base = pathlib.Path(root)
    tracks = []
    for lane in expected:
        source = base / f"hcc-p0-{lane}" / "track_status.json"
        status = "MISSING_ARTIFACT"
        tasks = []
        if source.is_file():
            try:
                doc = json.loads(source.read_text(encoding="utf-8"))
                if (doc.get("schema_version") == "hcc-p0-parallel-v1"
                    and doc.get("track") == lane
                    and isinstance(doc.get("tasks"), list)
                    and doc.get("task_count") == len(doc["tasks"])
                    and doc.get("task_count") >= 1
                    and doc.get("biological_claim_supported") is False):
                    tasks = doc["tasks"]
                    all_good = all(t.get("result") == "PASS_EXECUTION" for t in tasks)
                    status = ("PASS_EXECUTION" if all_good and
                              doc.get("status") == "PASS_EXECUTION"
                              else "FAIL_EXECUTION")
                else:
                    status = "MALFORMED_ARTIFACT"
            except (OSError, ValueError, TypeError):
                status = "MALFORMED_ARTIFACT"
        tracks.append({
            "track": lane, "status": status, "task_count": len(tasks),
            "task_pass": sum(t.get("result") == "PASS_EXECUTION" for t in tasks),
            "task_fail": sum(t.get("result") != "PASS_EXECUTION" for t in tasks),
        })
    success = all(t["status"] == "PASS_EXECUTION" for t in tracks)
    return {
        "schema_version": "hcc-p0-parallel-summary-v1",
        "evaluated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "overall_ci_status": "PASS_EXECUTION" if success else "FAIL_OR_MISSING",
        "lane_count": len(tracks),
        "lane_success_count": sum(t["status"] == "PASS_EXECUTION" for t in tracks),
        "task_pass_count": sum(t["task_pass"] for t in tracks),
        "task_failed_count": sum(t["task_fail"] for t in tracks),
        "tracks": tracks,
        "V02_patient_quantitative_join": "BLOCKED_UNTIL_PATIENT_LEVEL_RNA_PROTEIN_VALUES_PROVEN",
        "V04_etiology_group_scientific_status": "PARTIAL_SOURCE_CODEBOOK_UNVERIFIED",
        "V12_recurrence_codebook_scientific_status": "PARTIAL_EVENT_MEANING_UNVERIFIED",
        "biological_go_not_inferred_from_ci": True,
        "interpretation": "Execution PASS verifies scripts and artifacts only; independent clinical and biological evidence gates remain unchanged.",
    }


def as_markdown(report):
    lines = [
        "# HCC Parallel P0 — GitHub Actions evidence audit", "",
        f"**Execution:** {report['overall_ci_status']} "
        f"({report['lane_success_count']}/{report['lane_count']} lanes, "
        f"{report['task_pass_count']} tasks passed)", "",
        "| Lane | Execution status | Tasks passed | Tasks failed |",
        "|---|---|---:|---:|",
    ]
    for t in report["tracks"]:
        lines.append(
            f"| {t['track']} | {t['status']} | {t['task_pass']} | {t['task_fail']} |"
        )
    lines += [
        "", "**Research gates remain on evidence:** V02 RNA–protein join BLOCKED; "
        "V04 etiology subgroup source codebook PARTIAL; "
        "V12 RFS/OS source codebook PARTIAL.",
        "A green CI is not an eligible omics discovery, a verified patient N, "
        "an event count, or a clinical prediction result.", ""
    ]
    return "\n".join(lines)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root", default="artifacts")
    p.add_argument("--out", default="results/hcc/parallel/summary.json")
    args=p.parse_args()
    report=summarize(args.root)
    out=pathlib.Path(args.out)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",
                   encoding="utf-8")
    summary=as_markdown(report)
    print(summary)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"],"a",encoding="utf-8") as fh:
            fh.write(summary+"\n")
    if report["overall_ci_status"] != "PASS_EXECUTION":
        raise SystemExit(1)


if __name__=="__main__":
    main()
