#!/usr/bin/env python3
"""Fail-closed static protocol validation for EUR/EAS ischemic-stroke study design.

This audits the protocol contract and evidence states, not the original GWAS or
eQTL records. Status "DESIGN_VALID" MUST NOT be described as scientific PASS.
"""
import argparse
import json
import os
import pathlib
import sys
from datetime import datetime, timezone

STATUS = {"PASS", "PARTIAL_HOLD", "HOLD", "BLOCKED", "SOURCE_UNAVAILABLE", "STOP"}
GATE_IDS = ["G0", "G1", "G2", "G3", "G4", "G5"]
REQUIRED_SOURCE_IDS = {
    "GWAS_EUR_IS", "GWAS_EAS_IS", "EQTL_EUR_BLOOD", "EQTL_EAS_BLOOD",
    "LD_EUR", "LD_EAS",
}
REQUIRED_STEP_IDS = ["A0_phenotype_audit", "A1_gwas_qc", "A2_eQTL_ld_qc",
                     "A3_locus_coloc", "A4_compare", "A5_novelty"]


def validate(design, phenotypes):
    errors = []
    def check(condition, message):
        if not condition:
            errors.append(message)

    check(design.get("schema_version") == "1.0", "Schema version must be 1.0")
    check(design.get("protocol_id") == "IS_EUR_EAS_REGULATORY_V1", "Unrecognized protocol")
    check(design.get("locked") is False, "Cannot claim pre-registration without a frozen design")
    check("clinical" in design.get("primary_outcome", {}).get("clinical_basis", "").lower()
          or "cohort" in design.get("primary_outcome", {}).get("clinical_basis", "").lower()
          or "TOAST" in design.get("primary_outcome", {}).get("clinical_basis", ""),
          "Primary outcome needs explicit clinical/ascertainment caveat")
    outcome = design.get("primary_outcome", {})
    check(outcome.get("type") == "binary_case_control_GWAS", "Primary endpoint must not be longitudinal outcome")
    check(outcome.get("name") == "ischemic_stroke_susceptibility", "Primary outcome cannot be AS/CE/PD")
    check(outcome.get("eas_gwas_accession") == "GCST90104545", "EAS IS accession mismatch")
    check("not_independent_replication" in json.dumps(design.get("data_contracts", [])),
          "EAS evaluation must not be called independent replication")
    check(design.get("comparison_design", {}).get("external_validation_claim") is False,
          "Do not claim independent validation using within-consortium strata")
    check(design.get("comparison_design", {}).get("discovery") == "EUR_IS_GWAS_ONLY",
          "Loci must be selected from discovery, not EAS evaluation")
    check(design.get("primary_endpoint", {}).get("unit") == "independent_GWAS_locus",
          "Primary inference unit must be a locus, not a person")
    check("NON_EVALUABLE" in design.get("primary_endpoint", {}).get("missingness_rule", ""),
          "Missing EAS coverage must not imply biological disagreement")
    sources = design.get("data_contracts", [])
    source_ids = [x.get("id") for x in sources]
    check(len(source_ids) == len(set(source_ids)) and set(source_ids) == REQUIRED_SOURCE_IDS,
          "Missing or duplicated source contract")
    for x in sources:
        check(bool(x.get("requires")) and bool(x.get("ancestry")) and bool(x.get("status")),
              "Source contract missing critical fields")
        if x.get("kind") in {"gwas", "eqtl"}:
            for key in ["genome_build", "variant_id", "effect_allele", "other_allele"]:
                check(key in x.get("requires", []), x.get("id", "?") + " missing " + key)
    gates = design.get("decision_gates", [])
    gate_ids = [x.get("id") for x in gates]
    check(gate_ids == GATE_IDS, "Gate order must be G0 through G5")
    check(design.get("run_order") == GATE_IDS, "Run order must match gate order")
    for g in gates:
        check(g.get("status") in STATUS, "Unknown gate status: " + str(g.get("status")))
        check(bool(g.get("required_to_pass")), "Gate has no executable acceptance criteria")
        if g.get("status") == "PASS":
            check(bool(g.get("evidence_url")) and bool(g.get("verified_at_utc")),
                  "PASS requires evidence_url and verified_at_utc: " + str(g.get("id")))
    check(gates[0].get("status") != "PASS", "Cohort diagnostic basis is not yet verified")
    check(gates[1].get("status") != "PASS", "GWAS bytes have not been verified")
    check(gates[2].get("status") != "PASS", "EAS full eQTL not yet audited")
    steps = design.get("analysis_methods", [])
    check([x.get("step") for x in steps] == REQUIRED_STEP_IDS, "Analysis methods order mismatch")
    for step in steps:
        for key in ("question", "method", "deliverable", "stop_rule"):
            check(bool(step.get(key)), step.get("step", "?") + " missing " + key)
    check(phenotypes.get("research_primary_outcome") == "ischemic_stroke",
          "Phenotype register must mark IS only as primary")
    ids = {x.get("accession"): x for x in phenotypes.get("eas_studies", [])}
    check(set(ids) == {"GCST90104544", "GCST90104545", "GCST90104546", "GCST90104547", "GCST90104548"},
          "Incorrect EAS phenotype mapping")
    for key, name in (("GCST90104544", "any_stroke"),
                      ("GCST90104545", "ischemic_stroke"),
                      ("GCST90104546", "cardioembolic_stroke"),
                      ("GCST90104547", "large_artery_stroke"),
                      ("GCST90104548", "small_vessel_stroke")):
        check(ids.get(key, {}).get("phenotype") == name, key + " phenotype mismatch")
    check(ids.get("GCST90104545", {}).get("primary") is True,
          "EAS IS accession must be primary")
    check(all(x.get("case_ascertainment") == "NOT_VERIFIED" for x in ids.values()),
          "No unverified cohort case ascertainment can be silently marked verified")
    return errors


def run(design_path, register_path, output_dir):
    design = json.loads(pathlib.Path(design_path).read_text(encoding="utf-8"))
    register = json.loads(pathlib.Path(register_path).read_text(encoding="utf-8"))
    errors = validate(design, register)
    gates = {x["id"]: x["status"] for x in design["decision_gates"]}
    report = {
        "time_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_id": design["protocol_id"],
        "result": "DESIGN_VALID" if not errors else "DESIGN_INVALID",
        "scientific_readiness": "HOLD",
        "evidence_scope": "STATIC_CONTRACT_ONLY_NOT_SOURCE_DATA_QC",
        "design_errors": errors,
        "gates": gates,
        "tasks_remaining": [k for k, v in gates.items() if v != "PASS"],
        "note": "Do not compute a percentage completed from gate count; no biological inference has run.",
    }
    output = pathlib.Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    (output / "study_design_gate.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    lines = [
        "# Ischemic stroke EUR/EAS research design status",
        "",
        "**Contract:** " + report["result"] + " / **Research:** HOLD",
        "",
        "| Gate | Status |",
        "|---|---|",
    ] + [f"| {k} | {v} |" for k, v in gates.items()]
    if errors:
        lines += ["", "## Contract errors"] + ["- " + x for x in errors]
    lines += ["", "This is a **design/provenance audit only**, not GWAS, eQTL, colocalization, or an independent replication."]
    summary = "\n".join(lines) + "\n"
    (output / "study_design_gate.md").write_text(summary, encoding="utf-8")
    print(summary)
    path = os.environ.get("GITHUB_STEP_SUMMARY")
    if path:
        with open(path, "a", encoding="utf-8") as handle:
            handle.write(summary)
    return 2 if errors else 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--design", default="config/research/is_eur_eas_regulatory_v1.json")
    p.add_argument("--phenotypes", default="data/metadata/stroke_eas_phenotype_register_2026-11-11.json")
    p.add_argument("--output-dir", default="artifacts/is_eas_design_gate")
    args = p.parse_args()
    return run(args.design, args.phenotypes, args.output_dir)


if __name__ == "__main__":
    sys.exit(main())
