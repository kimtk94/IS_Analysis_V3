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
    ap.add_argument("--discovery-mr", type=Path,
                    help="Stage 7 discovery ancestry MR summary")
    ap.add_argument("--replication-mr", type=Path,
                    help="Stage 7 replication ancestry MR summary")
    ap.add_argument("--cross-ancestry-output", type=Path,
                    help="Stage 7 cross-ancestry output TSV")
    ap.add_argument("--discovery-ancestry", choices=["EUR","EAS","AFR","SAS","AMR"])
    ap.add_argument("--replication-ancestry", choices=["EUR","EAS","AFR","SAS","AMR"])
    ap.add_argument("--discovery-coloc", type=Path)
    ap.add_argument("--replication-coloc", type=Path)
    ap.add_argument("--platform-discovery-mr", type=Path,
                    help="Stage 8 discovery-platform MR summary")
    ap.add_argument("--platform-replication-mr", type=Path,
                    help="Stage 8 replication-platform MR summary")
    ap.add_argument("--platform-output", type=Path,
                    help="Stage 8 cross-platform output TSV")
    ap.add_argument("--discovery-platform")
    ap.add_argument("--replication-platform")
    ap.add_argument("--platform-mapping", type=Path)
    ap.add_argument("--platform-discovery-id-column")
    ap.add_argument("--platform-replication-id-column")
    ap.add_argument("--platform-discovery-coloc", type=Path)
    ap.add_argument("--platform-replication-coloc", type=Path)
    ap.add_argument("--transcript-protein-mr", type=Path,
                    help="Stage 9 protein MR summary")
    ap.add_argument("--smr-results", type=Path,
                    help="Stage 9 SMR/HEIDI result table")
    ap.add_argument("--transcript-phenotype")
    ap.add_argument("--transcript-output", type=Path)
    ap.add_argument("--eqtl-coloc", type=Path)
    ap.add_argument("--probe-map", type=Path)
    ap.add_argument("--phenotype-manifest", type=Path,
                    help="Stage 10/12 manifest of MR and optional coloc files")
    ap.add_argument("--phenotype-long-output", type=Path)
    ap.add_argument("--phenotype-wide-output", type=Path)
    ap.add_argument("--risk-disease-mr", type=Path)
    ap.add_argument("--risk-mr", type=Path)
    ap.add_argument("--risk-disease-phenotype")
    ap.add_argument("--risk-output", type=Path)
    ap.add_argument("--risk-sign-map", type=Path)
    ap.add_argument("--individual-subject-file", type=Path)
    ap.add_argument("--individual-predictor")
    ap.add_argument("--individual-output", type=Path)
    ap.add_argument("--baseline-endpoint")
    ap.add_argument("--slope-endpoint")
    ap.add_argument("--incident-endpoint")
    ap.add_argument("--individual-covariates", default="")
    ap.add_argument("--time-to-event")
    ap.add_argument("--event-col")
    ap.add_argument("--longitudinal-file", type=Path)
    ap.add_argument("--id-col", default="participant_id")
    ap.add_argument("--time-col", default="time_years")
    ap.add_argument("--egfr-col", default="egfr")
    ap.add_argument("--localization-candidates", type=Path)
    ap.add_argument("--tissue-manifest", type=Path)
    ap.add_argument("--cell-manifest", type=Path)
    ap.add_argument("--spatial-manifest", type=Path)
    ap.add_argument("--localization-output", type=Path)
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
    cross_args = [args.discovery_mr, args.replication_mr, args.cross_ancestry_output,
                  args.discovery_ancestry, args.replication_ancestry]
    if 7 in stages and all(x is not None for x in cross_args):
        cmd = [sys.executable, "scripts/run_master_cross_ancestry.py",
               "--discovery-mr", str(args.discovery_mr),
               "--replication-mr", str(args.replication_mr),
               "--output", str(args.cross_ancestry_output),
               "--discovery-ancestry", str(args.discovery_ancestry),
               "--replication-ancestry", str(args.replication_ancestry)]
        if args.discovery_coloc:
            cmd += ["--discovery-coloc", str(args.discovery_coloc)]
        if args.replication_coloc:
            cmd += ["--replication-coloc", str(args.replication_coloc)]
        plans = [StagePlan(p.stage, p.name, "READY", cmd,
                           "Cross-ancestry direction/effect/coloc replication summary.")
                 if p.stage == 7 else p for p in plans]
    platform_args = [args.platform_discovery_mr, args.platform_replication_mr,
                     args.platform_output, args.discovery_platform, args.replication_platform]
    if 8 in stages and all(x is not None for x in platform_args):
        cmd = [sys.executable, "scripts/run_master_cross_platform.py",
               "--discovery-mr", str(args.platform_discovery_mr),
               "--replication-mr", str(args.platform_replication_mr),
               "--output", str(args.platform_output),
               "--discovery-platform", str(args.discovery_platform),
               "--replication-platform", str(args.replication_platform)]
        if args.platform_mapping:
            cmd += ["--mapping", str(args.platform_mapping)]
        if args.platform_discovery_id_column:
            cmd += ["--discovery-id-column", str(args.platform_discovery_id_column)]
        if args.platform_replication_id_column:
            cmd += ["--replication-id-column", str(args.platform_replication_id_column)]
        if args.platform_discovery_coloc:
            cmd += ["--discovery-coloc", str(args.platform_discovery_coloc)]
        if args.platform_replication_coloc:
            cmd += ["--replication-coloc", str(args.platform_replication_coloc)]
        plans = [StagePlan(p.stage, p.name, "READY", cmd,
                           "Cross-platform Olink/SomaScan-style replication summary.")
                 if p.stage == 8 else p for p in plans]
    transcript_args = [args.transcript_protein_mr, args.smr_results,
                       args.transcript_phenotype, args.transcript_output]
    if 9 in stages and all(x is not None for x in transcript_args):
        cmd = [sys.executable, "scripts/run_master_transcriptomics.py",
               "--protein-mr", str(args.transcript_protein_mr),
               "--smr", str(args.smr_results),
               "--phenotype", str(args.transcript_phenotype),
               "--output", str(args.transcript_output)]
        if args.eqtl_coloc:
            cmd += ["--eqtl-coloc", str(args.eqtl_coloc)]
        if args.probe_map:
            cmd += ["--probe-map", str(args.probe_map)]
        plans = [StagePlan(p.stage, p.name, "READY", cmd,
                           "Integrate protein MR with SMR/HEIDI and eQTL colocalization.")
                 if p.stage == 9 else p for p in plans]
    matrix_args = [args.phenotype_manifest, args.phenotype_long_output, args.phenotype_wide_output]
    if any(s in stages for s in (10,12)) and all(x is not None for x in matrix_args):
        cmd = [sys.executable, "scripts/run_master_phenotype_matrix.py",
               "--manifest", str(args.phenotype_manifest),
               "--long-output", str(args.phenotype_long_output),
               "--wide-output", str(args.phenotype_wide_output)]
        plans = [StagePlan(p.stage, p.name, "READY", cmd,
                           "Shared phenotype/subtype MR+coloc matrix engine.")
                 if p.stage in {10,12} else p for p in plans]

    risk_args = [args.risk_disease_mr, args.risk_mr, args.risk_disease_phenotype, args.risk_output]
    if 11 in stages and all(x is not None for x in risk_args):
        cmd = [sys.executable, "scripts/run_master_risk_factor.py",
               "--disease-mr", str(args.risk_disease_mr),
               "--risk-mr", str(args.risk_mr),
               "--disease-phenotype", str(args.risk_disease_phenotype),
               "--output", str(args.risk_output)]
        if args.risk_sign_map:
            cmd += ["--risk-sign-map", str(args.risk_sign_map)]
        plans = [StagePlan(p.stage, p.name, "READY", cmd,
                           "Protein-to-risk-factor mechanism support summary.")
                 if p.stage == 11 else p for p in plans]
    indiv_args = [args.individual_subject_file, args.individual_predictor, args.individual_output]
    if 13 in stages and all(x is not None for x in indiv_args):
        cmd = ["Rscript", "scripts/run_master_individual_validation.R",
               "--subject-file", str(args.individual_subject_file),
               "--predictor", str(args.individual_predictor),
               "--output", str(args.individual_output)]
        if args.baseline_endpoint:
            cmd += ["--baseline-endpoint", str(args.baseline_endpoint)]
        if args.slope_endpoint:
            cmd += ["--slope-endpoint", str(args.slope_endpoint)]
        if args.incident_endpoint:
            cmd += ["--incident-endpoint", str(args.incident_endpoint)]
        if args.individual_covariates:
            cmd += ["--covariates", str(args.individual_covariates)]
        if args.time_to_event and args.event_col:
            cmd += ["--time-to-event", str(args.time_to_event), "--event", str(args.event_col)]
        if args.longitudinal_file:
            cmd += ["--longitudinal-file", str(args.longitudinal_file),
                    "--id-col", str(args.id_col), "--time-col", str(args.time_col),
                    "--egfr-col", str(args.egfr_col)]
        plans = [StagePlan(p.stage, p.name, "READY", cmd,
                           "Individual-level baseline/slope/incident CKD validation; optional Cox/LME.")
                 if p.stage == 13 else p for p in plans]
    loc_args = [args.localization_candidates, args.localization_output]
    if any(s in stages for s in (14,15,16)) and all(x is not None for x in loc_args):
        cmd = [sys.executable, "scripts/run_master_localization.py",
               "--candidates", str(args.localization_candidates),
               "--output", str(args.localization_output)]
        if args.tissue_manifest:
            cmd += ["--tissue-manifest", str(args.tissue_manifest)]
        if args.cell_manifest:
            cmd += ["--cell-manifest", str(args.cell_manifest)]
        if args.spatial_manifest:
            cmd += ["--spatial-manifest", str(args.spatial_manifest)]
        plans = [StagePlan(p.stage, p.name, "READY", cmd,
                           "Shared tissue/single-cell/spatial localization integrator.")
                 if p.stage in {14,15,16} else p for p in plans]
    print_plan(args.disease, plans)
    if not args.execute:
        return 0
    return execute(plans, args.allow_blocked)

if __name__ == "__main__":
    raise SystemExit(main())
