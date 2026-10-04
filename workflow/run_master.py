#!/usr/bin/env python3
"""MASTER pipeline orchestrator.

The runner is deliberately conservative:
- --dry-run plans commands without executing them.
- --execute runs only stages marked executable in the registry.
- unavailable stages are reported as BLOCKED rather than silently skipped.
"""
from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

STAGE_NAMES = [
    "dataset_audit","exposure","instrument_qc","primary_mr","mr_robustness",
    "coloc","finemap","cross_ancestry","cross_platform","transcriptomics",
    "phenotype_expansion","risk_factor","disease_subtype","individual_validation",
    "tissue","single_cell","spatial","phewas","druggability","evidence_integration",
]

@dataclass(frozen=True)
class StagePlan:
    stage: int
    name: str
    status: str
    command: list[str] | None
    note: str = ""

def load_config(disease: str) -> dict:
    path = ROOT / "config" / f"master_{disease}.json"
    with path.open(encoding="utf-8") as f:
        return json.load(f)

def ckd_registry() -> dict[int, StagePlan]:
    py = sys.executable
    return {
        0: StagePlan(0,"dataset_audit","READY",[py,"scripts/inspect_ckd_stage2_pqtl.py"],
                     "Existing CKD pQTL audit entry point."),
        1: StagePlan(1,"exposure","READY",[py,"scripts/extract_ukbppp_st16_cis.py"],
                     "UKB-PPP ST16 cis exposure extraction."),
        2: StagePlan(2,"instrument_qc","READY",[py,"scripts/run_ckd_stage2_candidates.py"],
                     "Existing candidate/instrument preparation; generic QC wrapper remains a refactor target."),
        3: StagePlan(3,"primary_mr","READY",[py,"scripts/run_ckd_stage1_mr.py"],
                     "Existing CKD primary MR engine."),
        4: StagePlan(4,"mr_robustness","BLOCKED",None,
                     "Dedicated shared WM/Egger/Q/Steiger/LOO runner not yet implemented."),
        5: StagePlan(5,"coloc","READY",["Rscript","scripts/run_ckd_stage2b_coloc.R"],
                     "Requires prepared stage2b inputs."),
        6: StagePlan(6,"finemap","READY",["Rscript","scripts/run_ckd_stage2c_susie.R"],
                     "Requires ancestry-matched LD prepared by prepare_ckd_stage2c_ld.py."),
    }

def is_registry() -> dict[int, StagePlan]:
    return {
        0: StagePlan(0,"dataset_audit","BLOCKED",None,
                     "IS dataset registry/audit runner is not yet standardized."),
        1: StagePlan(1,"exposure","READY",[sys.executable,"scripts/extract_ukbppp_st16_cis.py"],
                     "Shared UKB-PPP exposure extraction."),
        2: StagePlan(2,"instrument_qc","BLOCKED",None,
                     "IS-specific outcome harmonization/QC wrapper required."),
        3: StagePlan(3,"primary_mr","BLOCKED",None,
                     "Generic MR runner required before IS execution."),
        4: StagePlan(4,"mr_robustness","BLOCKED",None,
                     "Shared WM/Egger/Q/Steiger/LOO runner required."),
        5: StagePlan(5,"coloc","BLOCKED",None,
                     "Generic coloc runner required."),
        6: StagePlan(6,"finemap","BLOCKED",None,
                     "Generic SuSiE runner + EAS/EUR LD routing required."),
    }

def registry(disease: str) -> dict[int, StagePlan]:
    return ckd_registry() if disease == "ckd" else is_registry()

def parse_stage_spec(text: str) -> list[int]:
    out: set[int] = set()
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            a,b = token.split("-",1)
            out.update(range(int(a),int(b)+1))
        else:
            out.add(int(token))
    bad = sorted(x for x in out if x < 0 or x >= len(STAGE_NAMES))
    if bad:
        raise ValueError(f"invalid stages: {bad}")
    return sorted(out)

def quote(cmd: list[str] | None) -> str:
    return " ".join(shlex.quote(x) for x in cmd) if cmd else "-"

def build_plan(disease: str, stages: list[int]) -> list[StagePlan]:
    cfg = load_config(disease)
    reg = registry(disease)
    plans: list[StagePlan] = []
    for s in stages:
        name = STAGE_NAMES[s]
        enabled = bool(cfg["stages"].get(name, False))
        if not enabled:
            plans.append(StagePlan(s,name,"DISABLED",None,"Disabled by disease config."))
        elif s in reg:
            plans.append(reg[s])
        else:
            plans.append(StagePlan(s,name,"PLANNED",None,"Stage contract exists; executable runner not yet wired."))
    return plans

def print_plan(disease: str, plans: list[StagePlan]) -> None:
    print(f"MASTER_PIPELINE disease={disease}")
    for p in plans:
        print(f"[{p.stage:02d}] {p.name:22s} {p.status:8s} {quote(p.command)}")
        if p.note:
            print(f"     note: {p.note}")

def execute(plans: list[StagePlan], allow_blocked: bool) -> int:
    blocked = [p for p in plans if p.status in {"BLOCKED","PLANNED"}]
    if blocked and not allow_blocked:
        print("ERROR: requested plan contains non-executable stages:", file=sys.stderr)
        for p in blocked:
            print(f"  stage {p.stage}: {p.name} ({p.status})", file=sys.stderr)
        return 2
    for p in plans:
        if p.status != "READY" or not p.command:
            continue
        print(f"RUN stage={p.stage} name={p.name}: {quote(p.command)}", flush=True)
        completed = subprocess.run(p.command, cwd=ROOT)
        if completed.returncode != 0:
            print(f"FAILED stage={p.stage} rc={completed.returncode}", file=sys.stderr)
            return completed.returncode
    return 0

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--disease", choices=["ckd","is"], required=True)
    ap.add_argument("--stages", default="0-6", help="e.g. 0-6 or 3,5,6")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--dry-run", action="store_true", help="plan only (default)")
    mode.add_argument("--execute", action="store_true", help="execute READY stages")
    ap.add_argument("--allow-blocked", action="store_true",
                    help="with --execute, run READY stages even when requested stages include BLOCKED/PLANNED stages")
    ap.add_argument("--harmonized-input", type=Path,
                    help="harmonized TSV/TSV.GZ for Stage 4 MR robustness")
    ap.add_argument("--robust-output", type=Path,
                    help="output TSV for Stage 4 MR robustness")
    ap.add_argument("--susie-input-dir", type=Path,
                    help="locus-wide summary-statistic directory for Stage 6")
    ap.add_argument("--ld-dir", type=Path,
                    help="ancestry-matched LD directory for Stage 6")
    ap.add_argument("--abf-file", type=Path,
                    help="default ABF coloc summary for Stage 6 comparison")
    ap.add_argument("--susie-output-dir", type=Path,
                    help="Stage 6 SuSiE output directory")
    ap.add_argument("--outcome-type", choices=["quant","cc","auto"], default="auto",
                    help="Stage 6 outcome type")
    ap.add_argument("--ancestry", choices=["EUR","EAS","AFR","SAS","AMR"],
                    help="ancestry used by the Stage 6 LD reference")
    args = ap.parse_args()
    stages = parse_stage_spec(args.stages)
    plans = build_plan(args.disease, stages)
    if 4 in stages and args.harmonized_input and args.robust_output:
        cmd = [sys.executable, "scripts/run_master_mr_robustness.py",
               "--input", str(args.harmonized_input), "--output", str(args.robust_output)]
        plans = [StagePlan(p.stage, p.name, "READY", cmd,
                           "Shared MR robustness: IVW/Q/weighted median/Egger/LOO/optional Steiger.")
                 if p.stage == 4 else p for p in plans]

    susie_args = [args.susie_input_dir, args.ld_dir, args.abf_file, args.susie_output_dir, args.ancestry]
    if 6 in stages and all(x is not None for x in susie_args):
        cmd = ["Rscript", "scripts/run_master_susie.R",
               str(args.susie_input_dir), str(args.ld_dir), str(args.abf_file),
               str(args.susie_output_dir), str(args.ancestry), args.outcome_type]
        plans = [StagePlan(p.stage, p.name, "READY", cmd,
                           "Generic ancestry-aware SuSiE/coloc.susie runner.")
                 if p.stage == 6 else p for p in plans]
    print_plan(args.disease, plans)
    if not args.execute:
        return 0
    return execute(plans, args.allow_blocked)

if __name__ == "__main__":
    raise SystemExit(main())
