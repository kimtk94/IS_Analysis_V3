#!/usr/bin/env python3
from pathlib import Path
import csv
import json
import os
import sys

ROOT = Path(os.environ.get("IS_ANALYSIS_ROOT", "/srv/is-analysis"))

INST = ROOT / "results/metabolic_resilience/stage3_full_pgwas/instruments/STAGE3C2_FROZEN_INSTRUMENTS_ALL.tsv"
ANN = ROOT / "results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3C4_DBSNP_RSID_ANNOTATION.tsv"

OUTDIR = ROOT / "results/metabolic_resilience/stage3_full_pgwas/instruments"
OUT = OUTDIR / "STAGE3C4_FROZEN_INSTRUMENTS_ANNOTATED.tsv"

AUDIT = ROOT / "results/metabolic_resilience/stage3_full_pgwas/audit"
SUMMARY = AUDIT / "STAGE3C4_INSTRUMENT_OVERLAY_SUMMARY.json"

OUTDIR.mkdir(parents=True, exist_ok=True)
AUDIT.mkdir(parents=True, exist_ok=True)

if not INST.exists() or not ANN.exists():
    payload = {
        "status": "HOLD_MISSING_INPUT",
        "frozen_instruments_exists": INST.exists(),
        "c4_annotation_exists": ANN.exists(),
    }
    SUMMARY.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    raise SystemExit(2)

with INST.open("r", encoding="utf-8", newline="") as f:
    inst = list(csv.DictReader(f, delimiter="\t"))

with ANN.open("r", encoding="utf-8", newline="") as f:
    ann = list(csv.DictReader(f, delimiter="\t"))

ann_by_vid = {}
for r in ann:
    vid = r.get("variant_id", "")
    if not vid:
        continue
    if r.get("status") != "PASS_RSID":
        continue
    rsid = r.get("rsid", "")
    if not rsid.startswith("rs"):
        continue
    ann_by_vid[vid] = r

out = []
missing = []

for r in inst:
    vid = r["variant_id"]
    a = ann_by_vid.get(vid)

    rr = dict(r)
    rr["rsid_c2_original"] = r.get("rsid", "")

    if a is None:
        rr["rsid"] = "."
        rr["rsid_annotation_source"] = ""
        rr["rsid_annotation_status"] = "MISSING_C4_ANNOTATION"
        missing.append(vid)
    else:
        rr["rsid"] = a["rsid"]
        rr["rsid_annotation_source"] = "NCBI_Variation_Services_dbSNP_GRCh37p13"
        rr["rsid_annotation_status"] = "PASS_RSID"

    out.append(rr)

unique_vids = {r["variant_id"] for r in inst}
unique_annotated = {r["variant_id"] for r in out if r["rsid_annotation_status"] == "PASS_RSID"}

if out:
    fields = list(out[0].keys())
else:
    fields = ["status"]

with OUT.open("w", encoding="utf-8", newline="") as f:
    wr = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    wr.writeheader()
    wr.writerows(out)

payload = {
    "status": "PASS" if len(unique_vids) == len(unique_annotated) and not missing else "FAIL",
    "instrument_rows": len(inst),
    "unique_variants": len(unique_vids),
    "unique_variants_rsid_annotated": len(unique_annotated),
    "rsid_annotation_pct": round(100 * len(unique_annotated) / len(unique_vids), 3) if unique_vids else 0,
    "missing_annotation_rows": len(missing),
    "canonical_variant_key": "GRCh37 chr:pos:REF:ALT",
    "instrument_definition_policy": "Stage3C2 frozen IV membership is unchanged; C4 adds rsID metadata only.",
    "output": str(OUT),
}

SUMMARY.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

raise SystemExit(0 if payload["status"] == "PASS" else 3)
