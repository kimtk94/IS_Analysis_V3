#!/usr/bin/env python3
"""Offline, fail-closed readiness ledger checks; never a biological PASS."""
from __future__ import annotations
import argparse
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
LEDGER = ROOT / "research/operations/daily_20261101_readiness.json"
IDS = ("is", "glaucoma", "hcc", "radiomics", "ocular", "muscle", "oral", "skin")

def require(ok, message):
    if not ok:
        raise ValueError(message)

def validate_all(data):
    require(data.get("schema") == "master-degree-eight-track-readiness-v1", "schema")
    require(data.get("ckd") == "HOLD_NO_EXECUTION", "CKD must not be run")
    tracks = data.get("tracks")
    require(isinstance(tracks, list) and len(tracks) == 8, "eight tracks")
    require(tuple(sorted(t["id"] for t in tracks)) == tuple(sorted(IDS)), "unique IDs")
    for t in tracks:
        require(t["audit_status"].startswith("PASS_"), "audit type not evidenced")
        require(t["scientific_gate"] in ("HOLD", "STOP"), "metadata cannot be biological PASS")
        require(t.get("task") and t.get("change"), "task and delta required")
        require(len(t.get("known_next_gates", [])) > 0, "next gates required")
        require(len(t.get("evidence_urls", [])) >= 2, "two evidence URLs required")
        require(all(u.startswith("https://") for u in t["evidence_urls"]), "invalid provenance URL")
    return {t["id"]: t for t in tracks}

def validate_track(t):
    o = t["observations"]
    k = t["id"]
    if k == "is":
        require(o["mode"] == "CI_ONLY_NO_SSH", "IS mode")
        code = (ROOT / "scripts/is/harness_ci_only.py").read_text()
        workflow = (ROOT / ".github/workflows/is-harness-ci-only.yml").read_text()
        for token in ("CI_ONLY_NO_SSH", "BLOCKED_REMOTE_INPUTS_NOT_AVAILABLE", "no_live_operations"):
            require(token in code, "IS fail-closed token: " + token)
        require("run: python3 scripts/is/harness_ci_only.py" in workflow, "IS workflow command")
        require(o["remote_reexecution"] == "NOT_VERIFIED", "IS remote execution must not be promoted")
    elif k == "glaucoma":
        low = o["low_expression"]
        adequate = o["adequate_expression"]
        require(set(low) == {"ERFE", "CHAC1", "PARM1"}, "low-expression identities")
        require(all(0 <= float(v) < o["threshold_logcpm"] for v in low.values()), "low-expression threshold")
        require(len(adequate) == len(set(adequate)) == 9, "adequate nine")
        require(o["discordant_total"] == len(low) + len(adequate) + o["remaining_others"] == 96, "96 partition")
    elif k == "hcc":
        require(o["cases"] == 171 and o["aliquots"] == 331, "PDC catalog denominators")
        require(o["tumor"] + o["normal"] + o["other"] == o["aliquots"], "biospecimen categories")
        require(o["unique_nonreference_channels"] + o["repeated_reference_channels"] == o["tmt_channels"] == 363, "TMT channels")
        require(o["paired_metadata_cases"] <= o["cases"], "metadata pairs")
        require(o["rna_join"] == "UNVERIFIED", "RNA join must remain unverified")
    elif k == "radiomics":
        require(o["patients"] * o["readers"] * o["voxel_grids"] * o["binwidths"] == o["extraction_rows"] == 168, "MRI factorial")
        require(o["patients"] * o["readers"] * o["features"] == o["frozen_feature_crosschecks"] == 112, "MRI frozen profiles")
        require(o["clinical_eligibility"] == "NOT_ASSESSED", "clinical gate")
    elif k == "ocular":
        require(o["runs"] == o["library_pairs"] * 2 == 46, "SRA paired libraries")
        require(o["fastq_files"] == o["runs"] * 2 == 92, "paired-end fastqs")
        require(o["alias_pairs"] + o["biosample_resolved_pairs"] == o["library_pairs"], "treatment pairing")
        require(o["tm_strain_labels"] == 9 and o["sce_strain_labels"] == 4, "strain labels")
        require(o["timepoint_mapping"] == "UNVERIFIED", "timepoints")
    elif k == "muscle":
        code = (ROOT / "scripts/muscle_stage2a_gtex.py").read_text()
        for token in ("SKIPPED_NO_GTEX_VARIANT", "b37VariantId", "variantId", "Muscle_Skeletal", "gtex_v10"):
            require(token in code, "GTEx source token: " + token)
        require(o["unresolved"] == "SKIPPED_NO_GTEX_VARIANT", "unresolved variants")
    elif k == "oral":
        require(o["accession"] == "GSE103322" and o["patients"] == 18 and o["cells"] == 5902, "oral cohort source")
        require(o["paired_proteomics"] == "UNVERIFIED", "oral proteome join")
    elif k == "skin":
        require(o["ukb_participants"] == 423992 and o["reported_loci"] == 74, "perceived-age study")
        require(o["full_sumstats_downloaded"] is False, "no claim of raw verification")
        require(o["readme_schema_verified"] is True and o["readme_effect_scale"] == "LOG_OR", "Bristol readme schema")
        require(set(o["readme_allele_columns"]) == {"ALLELE1", "ALLELE0", "A1FREQ"}, "allele fields")
    else:
        raise ValueError("unknown track")
    return {"track": k, "audit_check": "PASS", "scientific_gate": t["scientific_gate"],
            "task": t["task"], "remaining_known_gates": len(t["known_next_gates"]),
            "verification_scope": "offline source and arithmetic checks only; no scientific inference"}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--track", choices=IDS + ("all",), required=True)
    p.add_argument("--out-dir", default="results/research_operations")
    a = p.parse_args()
    data = json.loads(LEDGER.read_text())
    tracks = validate_all(data)
    selected = IDS if a.track == "all" else (a.track,)
    results = [validate_track(tracks[k]) for k in selected]
    dest = pathlib.Path(a.out_dir)
    dest.mkdir(parents=True, exist_ok=True)
    for r in results:
        (dest / (r["track"] + "_audit.json")).write_text(json.dumps(r, indent=2, ensure_ascii=False) + "\n")
        print(json.dumps(r, ensure_ascii=False))
    print("AUDIT_ONLY: scientific PASS not granted")

if __name__ == "__main__":
    main()
