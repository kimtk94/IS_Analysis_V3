#!/usr/bin/env python3
"""Read-only ROSMAP brain single-nucleus eQTL provenance/access audit.

Uses public Synapse entity metadata only; never attempts to bypass access
restrictions or download controlled human data. Outcome is provenance/access
status, NOT variant-level eQTL data eligibility.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import urllib.error
import urllib.parse
import urllib.request
import os

SOURCES = {
    "Astrocyte_subtype_6": "syn52363813",
    "Microglia_subtype_2": "syn52363851",
    "Excitatory_neuron_subtype_13": "syn52363818",
    "Inhibitory_neuron_subtype_12": "syn52363834",
}
BASE = "https://repo-prod.prod.sagebase.org/repo/v1/entity/"
CAP = 128000


def parse_entity(payload, entity_id):
    if payload.get("id") != entity_id:
        raise ValueError("Synapse entity ID mismatch")
    return {
        "synapse_id": entity_id,
        "entity_name": payload.get("name"),
        "entity_type": payload.get("concreteType"),
        "file_handle_present": bool(payload.get("dataFileHandleId")),
        "content_size_bytes": payload.get("contentSize"),
        "entity_modified_on": payload.get("modifiedOn"),
        "entity_created_on": payload.get("createdOn"),
        "entity_etag": payload.get("etag"),
        "metadata_source": BASE + entity_id,
    }


def audit_one(name, entity_id):
    item = {
        "cell_subtype": name,
        "synapse_id": entity_id,
        "source_study": "ROSMAP",
        "tissue": "dorsolateral_prefrontal_cortex",
        "assay": "single_nucleus_eQTL_summary_statistics",
        "record": f"https://www.synapse.org/Synapse:{entity_id}",
        "metadata_status": "SOURCE_UNAVAILABLE",
        "download_access": "NOT_TESTED",
        "full_cis_eqtl_coverage": "HOLD",
        "grch37_allele_harmonization": "HOLD",
        "non_AD_stroke_donor_state": "NOT_APPLICABLE",
    }
    req = urllib.request.Request(
        BASE + entity_id,
        headers={"Accept": "application/json", "User-Agent": "IS_Analysis_V3_ROSMAP_EvidenceAudit/1.0"},
    )
    try:
        with urllib.request.urlopen(req, timeout=20) as response:
            final = urllib.parse.urlparse(response.geturl())
            if final.scheme != "https" or final.hostname != "repo-prod.prod.sagebase.org":
                raise ValueError("Unexpected Synapse metadata redirect")
            content = response.read(CAP + 1)
            if len(content) > CAP:
                raise ValueError("Metadata exceeds allowed bytes")
            payload = json.loads(content)
        item.update(parse_entity(payload, entity_id))
        item["metadata_status"] = "PUBLIC_ENTITY_METADATA_ACCESSIBLE"
        item["download_access"] = "NOT_TESTED_MAY_REQUIRE_SYNapse_LOGIN_DUC"
    except urllib.error.HTTPError as exc:
        item["metadata_status"] = "AUTH_OR_ACCESS_RESTRICTION" if exc.code in (401, 403) else "HTTP_UNAVAILABLE"
        item["http_status"] = exc.code
        item["download_access"] = "HOLD_LOGIN_OR_DUC_REQUIRED" if exc.code in (401, 403) else "NOT_TESTED"
    except (urllib.error.URLError, TimeoutError, OSError, ValueError, json.JSONDecodeError) as exc:
        item["error"] = f"{type(exc).__name__}: {str(exc)[:140]}"
    return item


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", default="artifacts/is_rosmap_source_provenance")
    args = p.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    results = [audit_one(k, v) for k, v in SOURCES.items()]
    evidence = {
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "SYNPASE_ENTITY_METADATA_ONLY",
        "variant_level_eQTL_analysis_performed": False,
        "donor_access_attempted": False,
        "consent_bypass_attempted": False,
        "tissue_context": "aged DLPFC, not ischemic-stroke lesion nor patient-matched IS tissue",
        "study_limitation": "Neuropathology/age/composition confounding; brain signal is tissue context rather than EAS stroke causality.",
        "results": results,
    }
    (out / "source_status.json").write_text(json.dumps(evidence, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    lines = ["# ROSMAP subtype eQTL sources (metadata only)", "",
             "| Cell subtype | Synapse | Metadata | File access / eQTL QC |",
             "|---|---|---|---|"]
    for row in results:
        lines.append(f"| {row['cell_subtype']} | {row['synapse_id']} | {row['metadata_status']} | HOLD |")
    lines += ["", "ROSMAP tissues are aged DLPFC and are not stroke lesions; genotype/allele/LD/consent/DUC and full variants remain unverified.",
              "This is not a completed colocalization, ROSMAP eQTL replication, or stroke cell-type-specific effect estimate."]
    summary = "\n".join(lines) + "\n"
    (out / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)
    if os.getenv("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as stream:
            stream.write(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
