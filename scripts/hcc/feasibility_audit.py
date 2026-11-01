#!/usr/bin/env python3
"""Metadata-only HCC feasibility audit. Never treats publication-level counts as joined patient data."""
import argparse
import datetime as dt
import json
import pathlib
import urllib.error
import urllib.request


def validate(registry):
    if registry.get("schema_version") != "hcc-feasibility-v1":
        raise ValueError("unsupported schema version")
    required = {
        "id", "name", "modalities", "url", "access", "next_action",
        "verified_pairing", "verified_pair_n", "verified_etiology_counts",
        "verified_recurrence_events", "verified_external_recurrence_events", "evidence_artifact"
    }
    ids = set()
    for item in registry.get("datasets", []):
        missing = required.difference(item)
        if missing:
            raise ValueError(f"{item.get('id', '?')}: missing {sorted(missing)}")
        if item["id"] in ids:
            raise ValueError(f"duplicate cohort: {item['id']}")
        ids.add(item["id"])
        if not item["url"].startswith("https://"):
            raise ValueError(f"{item['id']}: requires HTTPS source")
        n = item["verified_pair_n"]
        events = item["verified_recurrence_events"]
        external = item["verified_external_recurrence_events"]
        if any(not isinstance(x, int) or isinstance(x, bool) or x < 0 for x in (n, events, external)):
            raise ValueError(f"{item['id']}: invalid verified patient/event counts")
        groups = item["verified_etiology_counts"]
        if not isinstance(groups, dict) or any(
            not isinstance(v, int) or isinstance(v, bool) or v < 0 for v in groups.values()
        ):
            raise ValueError(f"{item['id']}: invalid verified etiology counts")
        if (item["verified_pairing"] or n or events or external or any(groups.values())) and not item["evidence_artifact"]:
            raise ValueError(f"{item['id']}: verified flags/counts require a reviewed evidence artifact path")
        if n and not item["verified_pairing"]:
            raise ValueError(f"{item['id']}: verified pair N requires verified_pairing=true")
    if not ids:
        raise ValueError("empty dataset registry")
    return True


def decide(registry):
    validate(registry)
    policy = registry["gate_policy"]
    cohorts = registry["datasets"]
    paired = [x for x in cohorts if x["verified_pairing"]
              and x["verified_pair_n"] >= policy["min_verified_paired_n"]]
    # Etiology subgroup N must be confirmed within patient-linked RNA/protein samples.
    groups = {}
    for cohort in paired:
        for group, count in cohort["verified_etiology_counts"].items():
            groups[group] = groups.get(group, 0) + count
    eligible_groups = sorted(g for g, n in groups.items()
                             if n >= policy["min_verified_etiology_group_n"])
    recurrence = [x for x in cohorts
                  if x["verified_recurrence_events"] >= policy["min_verified_recurrence_events"]]
    external = [x for x in cohorts
                if x["verified_external_recurrence_events"] >= policy["min_verified_external_recurrence_events"]]
    gates = {
        "file_level_rna_protein_pairs": bool(paired),
        "two_sufficient_etiology_groups": len(eligible_groups) >= policy["min_etiology_groups"],
        "event_defined_recurrence_dataset": bool(recurrence),
        "independent_external_recurrence_dataset": bool(external),
    }
    mechanism_go = gates["file_level_rna_protein_pairs"] and gates["two_sufficient_etiology_groups"]
    # External validation must be from a different cohort than the event discovery set.
    independent_external = [x for x in external if any(y["id"] != x["id"] for y in recurrence)]
    gates["independent_external_recurrence_dataset"] = bool(independent_external)
    clinical_go = gates["event_defined_recurrence_dataset"] and gates["independent_external_recurrence_dataset"]
    return {
        "mechanistic_decision": "GO" if mechanism_go else "HOLD",
        "recurrence_prediction_decision": "GO" if mechanism_go and clinical_go else "HOLD",
        "overall_decision": "GO" if mechanism_go and clinical_go else "HOLD",
        "gate_results": gates,
        "eligible_etiology_groups": eligible_groups,
        "verified_paired_cohorts": [x["id"] for x in paired],
        "verified_recurrence_cohorts": [x["id"] for x in recurrence],
        "verified_external_cohorts": [x["id"] for x in external],
        "interpretation": "HOLD means missing file-level evidence, not proof the study is impossible."
    }


def probe_urls(registry, limit=5, timeout=8):
    """Link reachability check only: HTTP success does NOT verify cohort data availability."""
    results = []
    for item in registry["datasets"][:limit]:
        url = item["url"]
        try:
            req = urllib.request.Request(url, method="HEAD",
                headers={"User-Agent": "HCC-Feasibility-Metadata-Audit/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                results.append({"id": item["id"], "url": url, "http_status": resp.status,
                                "link_status": "REACHABLE"})
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as exc:
            # HTTP blocks / network outages are unknown; no access claims.
            status = getattr(exc, "code", None)
            results.append({"id": item["id"], "url": url, "http_status": status,
                            "link_status": "UNVERIFIED", "error": str(exc)[:180]})
    return results


def render(registry, result, probes):
    lines = [
        "# HCC Multi-omics Initial Feasibility Audit",
        "",
        f"Metadata registry reviewed: {registry['reviewed_at']}",
        f"Audit generated (UTC): {dt.datetime.now(dt.timezone.utc).isoformat()}",
        "",
        "**Evidence level:** literature / repository metadata only; patient-level join NOT performed.",
        "**Caution:** N is publication-reported specimens/samples unless explicitly labeled verified patient count.",
        "",
        f"**Mechanistic feasibility: {result['mechanistic_decision']}**",
        f"**Recurrence prediction feasibility: {result['recurrence_prediction_decision']}**",
        f"**Overall: {result['overall_decision']}**",
        "",
        "## Hard gates (require file-level evidence)",
    ]
    for name, ok in result["gate_results"].items():
        lines.append(f"- {'PASS' if ok else 'UNVERIFIED'}: {name}")
    lines += ["", "## Dataset inventory", "",
              "| Cohort | Omics | Published N | RNA/protein matched in paper | Verified pair N | Status |",
              "|---|---|---:|---:|---:|---|"]
    for item in registry["datasets"]:
        n = item.get("subjects_reported")
        reported = str(n) if n is not None else "unknown"
        paper_pairs = item.get("publication_paired_rna_protein_specimens")
        pairs = str(paper_pairs) if paper_pairs is not None else "unreported"
        lines.append(
            f"| [{item['id']}]({item['url']}) | {', '.join(item['modalities'])} | {reported} | "
            f"{pairs} | {item['verified_pair_n']} | {'FILE VERIFIED' if item['verified_pairing'] else 'PUBLICATION ONLY / NOT VERIFIED'} |"
        )
    lines += ["", "## Data access blockers / next actions", ""]
    for item in registry["datasets"]:
        lines.extend([f"### {item['id']}", f"- Access: {item['access']}",
                      f"- Next: {item['next_action']}", ""])
    if probes:
        lines += ["## Metadata URL reachability (NOT file availability)", ""]
        for p in probes:
            lines.append(f"- {p['id']}: {p['link_status']} HTTP={p.get('http_status')}")
        lines.append("")
    lines += [
        "## Do not overclaim",
        "- A paper's multi-omics sample N is not a file-level intersection.",
        "- GSE14520's 488 microarrays are not 488 distinct patients.",
        "- TCGA RPPA is a targeted antibody panel, not an unbiased global MS proteome.",
        "- Microvascular invasion (MVI), overall survival (OS), and treatment response are not interchangeable with recurrence-free survival (RFS).",
        "- HBV/MASH CD45+ single-cell profiles cannot directly localize malignant hepatocyte signals.",
        "- Controlled EGA transcriptomes must not be ingested into public CI.",
        "",
        "## Promotion rules",
        "Populate a verified field ONLY after committing a non-sensitive manifest or audit summary "
        "that documents the patient/sample join, disease etiology, endpoint time/event definitions and counts. "
        "Set `evidence_artifact` to its repository path. Keep protected patient IDs outside the repository.",
    ]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", default="research/hcc/feasibility/datasets.json")
    parser.add_argument("--out-dir", default="results/hcc/feasibility")
    parser.add_argument("--probe-links", action="store_true")
    args = parser.parse_args()
    with open(args.registry, encoding="utf-8") as fh:
        registry = json.load(fh)
    result = decide(registry)
    probes = probe_urls(registry) if args.probe_links else []
    output = pathlib.Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    (output / "audit.json").write_text(json.dumps({
        "reviewed_at": registry["reviewed_at"], "decision": result, "url_probes": probes
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (output / "audit.md").write_text(render(registry, result, probes), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
