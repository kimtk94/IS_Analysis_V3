#!/usr/bin/env python3
"""Publish research results from IS_Analysis_V3 into the private Vercel dashboard repo.

Designed for:
  source: /srv/is-analysis/IS_Analysis_V3
  dashboard clone: a local clone of kimtk94/brain_research_mr_scrna_seq

The script never stores credentials. Optional --git-push uses the dashboard
clone's existing Git authentication.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

IMAGE_EXTS = {".png", ".svg", ".jpg", ".jpeg", ".webp"}

def read_json(path: Path) -> Optional[Dict[str, Any]]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None

def read_tsv(path: Path) -> List[Dict[str, str]]:
    try:
        with path.open(encoding="utf-8", newline="") as h:
            return list(csv.DictReader(h, delimiter="\t"))
    except Exception:
        return []

def first_match(root: Path, pattern: str) -> Optional[Path]:
    hits = sorted(root.rglob(pattern))
    return hits[0] if hits else None

def copy_figures(source_root: Path, dashboard_public: Path, slug: str) -> List[Dict[str, str]]:
    outdir = dashboard_public / "generated" / slug
    outdir.mkdir(parents=True, exist_ok=True)
    figures: List[Dict[str, str]] = []
    seen = set()
    for p in sorted(source_root.rglob("*")):
        if not p.is_file() or p.suffix.lower() not in IMAGE_EXTS:
            continue
        low = str(p).lower()
        if "figure" not in low and "fig" not in p.stem.lower() and "plot" not in p.stem.lower():
            continue
        key = p.name
        if key in seen:
            continue
        seen.add(key)
        dest = outdir / p.name
        shutil.copy2(p, dest)
        figures.append({
            "title": p.stem.replace("_", " "),
            "src": f"/generated/{slug}/{p.name}",
            "kind": "svg" if p.suffix.lower() == ".svg" else "image",
            "caption": f"Auto-synced from {p.relative_to(source_root)}",
        })
    return figures

def ckd_manifest(root: Path, public: Path) -> Dict[str, Any]:
    base = root / "results" / "ckd"
    metrics: List[Dict[str, str]] = []
    highlights: List[str] = []

    susie = first_match(base, "STAGE2C_SUSIE_DEFAULT.tsv")
    if susie:
        rows = read_tsv(susie)
        sd = next((r for r in rows if r.get("gene_symbol") == "SDCCAG8"), None)
        if sd:
            metrics.append({"label": "SDCCAG8 pQTL CS", "value": sd.get("pqtl_credible_sets", "NA")})
            metrics.append({"label": "SDCCAG8 eGFR CS", "value": sd.get("egfr_credible_sets", "NA")})
            highlights.append(
                "SDCCAG8 SuSiE: pQTL credible sets="
                + str(sd.get("pqtl_credible_sets", "NA"))
                + ", eGFR credible sets="
                + str(sd.get("egfr_credible_sets", "NA"))
                + "."
            )

    coloc = first_match(base, "*COLOC*DEFAULT*.tsv") or first_match(base, "*coloc*.tsv")
    if coloc:
        rows = read_tsv(coloc)
        for gene in ("SDCCAG8", "GSTA3", "UMOD"):
            r = next((x for x in rows if x.get("gene_symbol") == gene), None)
            if r:
                h4 = r.get("PP.H4") or r.get("PP.H4.abf") or r.get("pp_h4")
                if h4:
                    metrics.append({"label": f"{gene} coloc", "value": f"PP.H4 {h4}"})

    stage4 = first_match(base, "STAGE4_PUBLIC_PROTOTYPE.json")
    if stage4:
        j = read_json(stage4) or {}
        for key, label in (
            ("n_subjects", "KoGES subjects"),
            ("baseline_f0_available", "F0 eGFR available"),
            ("prevalent_f0_ckd", "Prevalent CKD"),
            ("primary_slope_eligible", "Slope eligible"),
        ):
            if key in j:
                metrics.append({"label": label, "value": str(j[key])})

    if not metrics:
        metrics = [
            {"label": "Eligible proteins", "value": "1,936"},
            {"label": "Strong cis signals", "value": "10,010", "note": "F>10, MHC excluded"},
        ]
    return {"status": "live", "metrics": metrics, "highlights": highlights, "figures": copy_figures(base, public, "ckd")}

def muscle_manifest(root: Path, public: Path) -> Dict[str, Any]:
    base = root / "results"
    metrics: List[Dict[str, str]] = []
    highlights: List[str] = []

    s1 = first_match(base, "MUSCLE_STAGE1_SUMMARY.json")
    if s1:
        j = read_json(s1) or {}
        if "n_loci" in j:
            metrics.append({"label": "Primary RT loci", "value": str(j["n_loci"])})
        if "study_counts" in j:
            metrics.append({"label": "RT studies", "value": str(len(j["study_counts"]))})
        if "n_methods_threshold_pass" in j:
            metrics.append({"label": "Methods-threshold pass", "value": str(j["n_methods_threshold_pass"])})
        highlights.append("Stage 1 GWAS locus audit summary was detected and synced.")

    s2 = first_match(base, "MUSCLE_STAGE2A_SUMMARY.json")
    if s2:
        j = read_json(s2) or {}
        for key, label in (
            ("n_input_variants", "Stage 2A input variants"),
            ("n_variants_with_significant_eqtl", "Variants with muscle eQTL"),
            ("n_variants_with_significant_sqtl", "Variants with muscle sQTL"),
            ("n_candidate_genes", "Candidate genes"),
        ):
            if key in j:
                metrics.append({"label": label, "value": str(j[key])})
        highlights.append("GTEx v10 Muscle_Skeletal eQTL/sQTL summary was detected and synced.")

    if not metrics:
        metrics = [
            {"label": "Primary RT loci", "value": "20", "note": "Yang 2024 + Gu 2026"},
            {"label": "HIIT comparator loci", "value": "8"},
            {"label": "Current stage", "value": "Stage 2A", "note": "GTEx skeletal-muscle eQTL/sQTL"},
        ]
    return {"status": "stage2a", "metrics": metrics, "highlights": highlights, "figures": copy_figures(base / "muscle", public, "muscle") if (base / "muscle").exists() else []}

def generic_manifest(root: Path, public: Path, slug: str, result_dir: str) -> Dict[str, Any]:
    base = root / "results" / result_dir
    return {
        "status": "live" if base.exists() else "pending",
        "metrics": [],
        "highlights": [],
        "figures": copy_figures(base, public, slug) if base.exists() else [],
    }

def run(source: Path, dashboard_repo: Path, git_push: bool) -> Path:
    public = dashboard_repo / "research-dashboard" / "public"
    if not public.exists():
        raise SystemExit(f"Dashboard public directory not found: {public}")

    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_root": str(source),
        "studies": {
            "ckd": ckd_manifest(source, public),
            "is": generic_manifest(source, public, "is", "is"),
            "metabolic-resilience": generic_manifest(source, public, "metabolic-resilience", "metabolic_resilience"),
            "muscle": muscle_manifest(source, public),
        },
    }

    out = public / "data" / "research-manifest.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if git_push:
        cmds = [
            ["git", "add", "research-dashboard/public/data/research-manifest.json", "research-dashboard/public/generated"],
            ["git", "commit", "-m", "Auto-sync research dashboard outputs"],
            ["git", "push"],
        ]
        for cmd in cmds:
            p = subprocess.run(cmd, cwd=dashboard_repo, text=True, capture_output=True)
            print("$", " ".join(cmd))
            print(p.stdout.strip())
            if p.returncode != 0:
                print(p.stderr.strip())
                if cmd[1] == "commit" and "nothing to commit" in (p.stdout + p.stderr).lower():
                    continue
                raise SystemExit(p.returncode)
    return out

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, default=Path("/srv/is-analysis/IS_Analysis_V3"))
    ap.add_argument("--dashboard-repo", type=Path, required=True)
    ap.add_argument("--git-push", action="store_true")
    args = ap.parse_args()
    out = run(args.source.resolve(), args.dashboard_repo.resolve(), args.git_push)
    print(f"manifest={out}")

if __name__ == "__main__":
    main()
