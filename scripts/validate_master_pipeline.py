#!/usr/bin/env python3
import argparse, json
from pathlib import Path

REQUIRED_TOP = {"schema_version","disease","title","ancestry","stages","defaults","outcomes","figures"}
REQUIRED_STAGES = {
    "dataset_audit","exposure","instrument_qc","primary_mr","mr_robustness","coloc","finemap",
    "cross_ancestry","cross_platform","transcriptomics","phenotype_expansion","risk_factor",
    "disease_subtype","individual_validation","tissue","single_cell","spatial","phewas",
    "druggability","evidence_integration"
}

def validate(path: Path):
    cfg = json.loads(path.read_text())
    missing = REQUIRED_TOP - set(cfg)
    if missing:
        raise ValueError(f"{path}: missing top-level keys: {sorted(missing)}")
    missing_stages = REQUIRED_STAGES - set(cfg["stages"])
    if missing_stages:
        raise ValueError(f"{path}: missing stage flags: {sorted(missing_stages)}")
    if cfg["schema_version"] != 1:
        raise ValueError(f"{path}: unsupported schema_version")
    if cfg["disease"] not in {"ckd","is"}:
        raise ValueError(f"{path}: unsupported disease")
    if cfg["defaults"]["f_stat_min"] < 10:
        raise ValueError(f"{path}: f_stat_min must be >=10")
    coloc = cfg["defaults"]["coloc"]
    for k in ("p1","p2","p12","strong_pp_h4"):
        if k not in coloc:
            raise ValueError(f"{path}: missing coloc.{k}")
    if not (0 < coloc["strong_pp_h4"] <= 1):
        raise ValueError(f"{path}: strong_pp_h4 must be in (0,1]")
    if not cfg["outcomes"]:
        raise ValueError(f"{path}: outcomes cannot be empty")
    if cfg["disease"] == "is" and not cfg.get("risk_factors"):
        raise ValueError(f"{path}: IS config requires risk_factors")
    if cfg["disease"] == "ckd" and "individual_validation" not in cfg:
        raise ValueError(f"{path}: CKD config requires individual_validation")
    return cfg

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("configs", nargs="+")
    args = ap.parse_args()
    for raw in args.configs:
        cfg = validate(Path(raw))
        print(f"PASS {raw}: disease={cfg['disease']} stages={sum(bool(v) for v in cfg['stages'].values())}")

if __name__ == "__main__":
    main()
