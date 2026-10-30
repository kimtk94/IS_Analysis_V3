#!/usr/bin/env python3
"""Validate HCC checkpoint registry and render an evidence-level status report.

This is protocol governance, not a biological test. PASS cannot be inferred from
file URLs, earlier publications, a successful GitHub Action, or an empty p value.
"""
import argparse
import collections
import datetime as dt
import json
import pathlib

ALLOWED = {"NOT_STARTED", "PARTIAL", "PASS", "BLOCKED", "FAIL"}
EXPECTED = [f"V{i:02d}" for i in range(16)]


def validate(reg):
    if reg.get("schema_version") != "hcc-validation-1.0":
        raise ValueError("wrong HCC checkpoint schema")
    rows = reg.get("checkpoints") or []
    ids = [item.get("id") for item in rows]
    if sorted(ids) != EXPECTED:
        raise ValueError("Missing, duplicated or unexpected HCC validation checkpoints")
    refs = {ref["id"] for ref in reg.get("references", [])}
    for row in rows:
        cid = row["id"]
        status = row.get("status")
        if status not in ALLOWED:
            raise ValueError(f"{cid}: unknown status: {status}")
        if not row.get("pass_criteria") or not row.get("next_action"):
            raise ValueError(f"{cid}: missing actionable acceptance criteria")
        if any(ref not in refs for ref in row.get("paper_refs", [])):
            raise ValueError(f"{cid}: unidentified literature benchmark")
        evidence = row.get("evidence", [])
        if any("url" not in x or not x["url"].startswith("https://") for x in evidence):
            raise ValueError(f"{cid}: evidence URLs must be HTTPS")
        if status == "PASS":
            if not row.get("reviewer_approved"):
                raise ValueError(f"{cid}: PASS without reviewed evidence")
            audit = row.get("analysis_evidence") or []
            if not audit:
                raise ValueError(f"{cid}: PASS without analysis artifact")
            if any(not x.get("path") or not x.get("source_sha256") or
                   not x.get("code_sha") or not x.get("checked_at") for x in audit):
                raise ValueError(f"{cid}: PASS requires path, source SHA256, code SHA and timestamp")
            if cid == "V02" and not (
                isinstance(row.get("independent_patient_n"), int)
                and row["independent_patient_n"] > 0
            ):
                raise ValueError("V02: paired matrix PASS requires independent patient N")
            if cid in ("V12", "V13") and not (
                isinstance(row.get("independent_event_n"), int)
                and row["independent_event_n"] > 0
            ):
                raise ValueError(f"{cid}: clinical PASS needs verified events")
        if status == "FAIL" and not row.get("analysis_evidence"):
            raise ValueError(f"{cid}: failure claim must also be evidence-backed")
        if status == "NOT_STARTED" and row.get("analysis_evidence"):
            raise ValueError(f"{cid}: evidence exists but stage marked not started")
    return True


def summarize(reg):
    validate(reg)
    counts = collections.Counter(x["status"] for x in reg["checkpoints"])
    blockers = [x for x in reg["checkpoints"] if x["status"] in ("BLOCKED", "PARTIAL")]
    return {"created_at": reg["created_at"], "evaluated_at_utc":
            dt.datetime.now(dt.timezone.utc).isoformat(),
            "counts": {key: counts[key] for key in sorted(ALLOWED)},
            "primary_claim_allowed": all(
                next(x for x in reg["checkpoints"] if x["id"] == key)["status"] == "PASS"
                for key in ("V02", "V03", "V04", "V05", "V07", "V08", "V09")
            ),
            "clinical_claim_allowed": all(
                next(x for x in reg["checkpoints"] if x["id"] == key)["status"] == "PASS"
                for key in ("V12", "V13")
            ),
            "unresolved": [{"id": x["id"], "title": x["title"],
                            "status": x["status"], "priority": x["priority"],
                            "next_action": x["next_action"]} for x in blockers]
            }


def to_markdown(result):
    lines = ["# HCC 연구 검증 체크포인트 상태", "",
             f"- Generated UTC: {result['evaluated_at_utc']}",
             "- **Note:** This is a guardrail report, NOT a molecular result.",
             "- Matched cross-etiology biology claim: " +
             ("ALLOWED (all gates passed)" if result["primary_claim_allowed"] else "NOT YET ALLOWED"),
             "- Independent clinical prognosis claim: " +
             ("ALLOWED (all gates passed)" if result["clinical_claim_allowed"] else "NOT YET ALLOWED"),
             "", "## Current gate distribution", ""]
    lines.extend(f"- {k}: {v}" for k, v in result["counts"].items())
    lines += ["", "## Next unresolved checkpoints", ""]
    lines.extend(
        f"- **{v['id']} {v['status']}** — {v['title']} ({v['priority']}): {v['next_action']}"
        for v in result["unresolved"])
    lines += ["", "To continue: open research/hcc/validation/HCC_VALIDATION_MASTER.md "
              "and the matching gate in checkpoints.json. Never upgrade PASS without "
              "verified artifacts, reviewer approval, and biological counts.", ""]
    return "\n".join(lines)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--registry", default="research/hcc/validation/checkpoints.json")
    p.add_argument("--out-dir", default="results/hcc/validation")
    args = p.parse_args()
    registry = json.loads(pathlib.Path(args.registry).read_text(encoding="utf-8"))
    summary = summarize(registry)
    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "validation_status.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (out / "validation_status.md").write_text(to_markdown(summary), encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "unresolved"}, indent=2))


if __name__ == "__main__":
    main()
