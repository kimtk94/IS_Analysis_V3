#!/usr/bin/env python3
"""Validate an optional IS tissue-enrichment design extension.

This checks method/source contracts and emits HOLD statuses. It does NOT run
MAGMA, stratified LDSC, GWAS, or expression enrichment.
"""
import argparse
import json
import os
from pathlib import Path
from datetime import datetime, timezone

EXPECTED_METHODS = ["MAGMA_GENE", "MAGMA_GENE_PROPERTY", "LDSC_SEG"]
EXPECTED_SOURCES = {
    "GWAS_EUR_IS", "GWAS_EAS_IS", "LD_EUR", "LD_EAS", "GENE_ANNOTATION",
    "EXPRESSION_SPECIFICITY", "LDSC_SEG_EUR_ANNOTATION_LDSCORES",
    "LDSC_SEG_EAS_ANNOTATION_LDSCORES", "LDSC_EUR_WEIGHTS", "LDSC_EAS_WEIGHTS"
}
ALLOWED_DECISIONS = {"ADD_AS_SECONDARY", "OPTIONAL_CONTEXT", "CONDITIONAL_SENSITIVITY_ONLY"}


def validate(ext, primary):
    errors = []

    def require(condition, error):
        if not condition:
            errors.append(error)

    require(ext.get("schema_version") == "1.0", "Unsupported extension schema")
    require(ext.get("extends_protocol_id") == primary.get("protocol_id"), "Primary protocol mismatch")
    require(ext.get("primary_endpoint_unchanged") is True, "Primary endpoint cannot change")
    require(ext.get("new_primary_hypothesis") is False, "Optional methods cannot introduce new primary hypothesis")
    require(primary.get("primary_outcome", {}).get("name") == "ischemic_stroke_susceptibility",
            "Primary outcome must remain ischemic stroke susceptibility")
    require(primary.get("primary_endpoint", {}).get("name") == "cross_ancestry_locus_to_gene_concordance",
            "Primary endpoint changed unexpectedly")
    require(primary.get("primary_outcome", {}).get("eur_gwas_accession") == "GCST90104540",
            "Incorrect EUR IS accession in core protocol")
    require(primary.get("primary_outcome", {}).get("eas_gwas_accession") == "GCST90104545",
            "Incorrect EAS IS accession in core protocol")

    sources = ext.get("sources", {})
    require(set(sources) == EXPECTED_SOURCES, "Missing or extra data sources")
    for sid, src in sources.items():
        require(isinstance(src.get("source_qc_verified"), bool), "source_qc_verified must be boolean: " + sid)
        if src.get("source_qc_verified") is True:
            require(bool(src.get("evidence_url")) and bool(src.get("verified_at_utc")),
                    "Cannot verify a source without evidence URL/time: " + sid)
    for sid, ancestry in [("GWAS_EUR_IS", "EUR"), ("GWAS_EAS_IS", "EAS"),
                           ("LD_EUR", "EUR"), ("LD_EAS", "EAS"),
                           ("LDSC_SEG_EUR_ANNOTATION_LDSCORES", "EUR"),
                           ("LDSC_SEG_EAS_ANNOTATION_LDSCORES", "EAS"),
                           ("LDSC_EUR_WEIGHTS", "EUR"), ("LDSC_EAS_WEIGHTS", "EAS")]:
        require(sources.get(sid, {}).get("ancestry") == ancestry, "Reference ancestry mismatch: " + sid)
    require(sources.get("GWAS_EUR_IS", {}).get("accession") == "GCST90104540",
            "EUR GWAS accession mismatch")
    require(sources.get("GWAS_EAS_IS", {}).get("accession") == "GCST90104545",
            "EAS GWAS accession mismatch")
    require(all(sources.get(sid, {}).get("phenotype") == "ischemic_stroke"
                for sid in ("GWAS_EUR_IS", "GWAS_EAS_IS")), "Wrong GWAS phenotype")

    methods = ext.get("analyses", [])
    require([x.get("id") for x in methods] == EXPECTED_METHODS, "Method order/name mismatch")
    required_roles = {
        "MAGMA_GENE": "ADD_AS_SECONDARY",
        "MAGMA_GENE_PROPERTY": "OPTIONAL_CONTEXT",
        "LDSC_SEG": "CONDITIONAL_SENSITIVITY_ONLY",
    }
    for x in methods:
        mid = x.get("id")
        require(x.get("decision") == required_roles.get(mid) and x.get("decision") in ALLOWED_DECISIONS,
                "Unsupported method role: " + str(mid))
        for key in ("question", "questions", "requirements", "outputs", "non_interpretations"):
            if key == "question":
                continue
            require(bool(x.get(key)), "Missing method field " + key + ": " + str(mid))
        require(x.get("eur_status") != "PASS" and x.get("eas_status") != "PASS",
                "Do not mark computational analysis PASS before running real sources: " + str(mid))
    if len(methods) == 3:
        magma, property_test, ldsc = methods
        require("LD_EUR" in magma.get("depends_on", []) and "GWAS_EUR_IS" in magma.get("depends_on", []),
                "MAGMA gene must require EUR GWAS and EUR LD")
        require("MAGMA_GENE" in property_test.get("depends_on", []) and
                "EXPRESSION_SPECIFICITY" in property_test.get("depends_on", []),
                "MAGMA gene-property needs full gene-level results and atlas")
        require("LDSC_SEG_EUR_ANNOTATION_LDSCORES" in ldsc.get("depends_on", []),
                "LDSC SEG requires matching EUR annotation LD scores")
        conditions = " ".join(ldsc.get("requirements", [])).lower()
        require("eas" in conditions and "eur" in conditions and "baseline" in conditions,
                "EAS LDSC SEG must explicitly prohibit EUR LD substitution")
        heuristic = ldsc.get("readiness_heuristic", {})
        require(heuristic.get("min_h2_z_for_consideration", 0) >= 3.0
                and "heuristic" in heuristic.get("interpretation", "").lower(),
                "h2 Z cutoff must be recorded as a heuristic, not biological significance")

    panel = ext.get("prespecified_tissue_panel", {})
    tissues = panel.get("proposed", [])
    require(bool(tissues) and len(tissues) == len(set(tissues)), "Missing or duplicate tissue atlas candidate")
    require(panel.get("status") == "CANDIDATE_PENDING_ATLAS_MAPPING",
            "Tissue panel is not yet matched to an atlas")
    require("primary" in " ".join(ext.get("novelty_constraints", [])).lower(),
            "Novelties must explicitly preserve primary question")
    return errors


def summarize(ext, core):
    sources = ext["sources"]
    core_gates = {x["id"]: x["status"] for x in core["decision_gates"]}
    def ready(*ids):
        return all(sources[s]["source_qc_verified"] for s in ids)

    m_eur = core_gates["G1"] == "PASS" and ready("GWAS_EUR_IS", "LD_EUR", "GENE_ANNOTATION")
    m_eas = core_gates["G1"] == "PASS" and ready("GWAS_EAS_IS", "LD_EAS", "GENE_ANNOTATION")
    gp_eur = m_eur and ready("EXPRESSION_SPECIFICITY")
    gp_eas = m_eas and ready("EXPRESSION_SPECIFICITY")
    # S-LDSC h2 reliability/baseline model validity are not in this manifest yet.
    # Even if input files later pass, the method remains REVIEW pending diagnostics.
    ldsc_eur_inputs = core_gates["G1"] == "PASS" and ready(
        "GWAS_EUR_IS", "LDSC_SEG_EUR_ANNOTATION_LDSCORES", "LDSC_EUR_WEIGHTS"
    )
    ldsc_eas_inputs = core_gates["G1"] == "PASS" and ready(
        "GWAS_EAS_IS", "LDSC_SEG_EAS_ANNOTATION_LDSCORES", "LDSC_EAS_WEIGHTS"
    )
    return {
        "MAGMA_GENE": {"EUR": "ELIGIBLE_FOR_RUNTIME_QC" if m_eur else "HOLD_SOURCE",
                       "EAS": "ELIGIBLE_FOR_RUNTIME_QC" if m_eas else "HOLD_SOURCE"},
        "MAGMA_GENE_PROPERTY": {"EUR": "ELIGIBLE_FOR_RUNTIME_QC" if gp_eur else "HOLD_SOURCE",
                                "EAS": "ELIGIBLE_FOR_RUNTIME_QC" if gp_eas else "HOLD_SOURCE"},
        "LDSC_SEG": {"EUR": "REVIEW_H2_BASELINE" if ldsc_eur_inputs else "HOLD_SOURCE",
                     "EAS": "REVIEW_H2_BASELINE" if ldsc_eas_inputs else "HOLD_SOURCE"}
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--extension", default="config/research/is_optional_magma_ldscseg_v1.json")
    parser.add_argument("--core", default="config/research/is_eur_eas_regulatory_v1.json")
    parser.add_argument("--output-dir", default="artifacts/is_eas_optional_enrichment")
    args = parser.parse_args()
    ext = json.loads(Path(args.extension).read_text(encoding="utf-8"))
    core = json.loads(Path(args.core).read_text(encoding="utf-8"))
    errors = validate(ext, core)
    readiness = summarize(ext, core) if not errors else {}
    report = {
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "extension_id": ext.get("extension_id"),
        "contract_status": "DESIGN_VALID" if not errors else "DESIGN_INVALID",
        "scientific_result": "NOT_RUN",
        "primary_endpoint_unchanged": ext.get("primary_endpoint_unchanged"),
        "readiness": readiness,
        "errors": errors,
        "note": "Source readiness only. This workflow never runs MAGMA or LDSC SEG."
    }
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "optional_enrichment_gate.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    lines = ["# MAGMA / LDSC-SEG optional research methods", "",
             "**Contract:** " + report["contract_status"] + " / **Biological analysis:** NOT_RUN", "",
             "| Module | EUR | EAS |", "|---|---|---|"]
    for method, status in readiness.items():
        lines.append(f"| {method} | {status['EUR']} | {status['EAS']} |")
    if errors:
        lines.extend(["", "## Errors"] + ["- " + x for x in errors])
    lines += ["", "CI success means only the methods and source evidence gates are structured correctly."]
    summary = "\n".join(lines) + "\n"
    (out / "optional_enrichment_gate.md").write_text(summary, encoding="utf-8")
    print(summary)
    if os.getenv("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as stream:
            stream.write(summary)
    return 2 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
