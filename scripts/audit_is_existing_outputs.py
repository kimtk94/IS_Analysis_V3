#!/usr/bin/env python3
"""Read-only audit of frozen IS outputs; never infer missing build, SE or LD provenance."""
from __future__ import annotations
import argparse, csv, json
from collections import Counter, defaultdict
from pathlib import Path

FILES = {
 "cross_ancestry": "stage4_cross_eas/GIGASTROKE_BBJ_HARMONIZED_VARIANTS.tsv",
 "signal_coloc": "stage5_functional/phase9c_convergence/COLOC_SUSIE_SIGNAL_PAIRS.tsv",
 "abf_coloc": "stage5_functional/phase9c_convergence/COLOC_ABF_MASTER_ANNOTATED_V2.tsv",
}
def audit(root):
    report = {"schema_version": 1, "read_only": True, "sources": {}, "loci": {}, "limitations": [
        "Genome build is not explicitly declared in cross-ancestry comparison input.",
        "Per-variant standard errors are absent from cross-ancestry comparison input.",
        "Reference LD ancestry/cohort and LD matrix validation are not declared in SuSiE signal-pair output.",
        "Colocalization support is provisional until molecular QTL LD, harmonization and convergence are independently verified.",
        "Existing results are not recomputed or reclassified as causal."
    ]}
    for key, rel in FILES.items():
        path = root / rel
        if not path.is_file():
            report["sources"][key] = {"status": "MISSING", "path": str(path)}
            continue
        with path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            fields = reader.fieldnames or []
            n = 0
            by_locus = Counter()
            classes = Counter()
            best = defaultdict(float)
            for row in reader:
                n += 1
                locus = row.get("locus", "UNSPECIFIED")
                by_locus[locus] += 1
                if key == "cross_ancestry":
                    classes[row.get("harmonization", "UNKNOWN")] += 1
                if key in ("signal_coloc", "abf_coloc"):
                    raw = row.get("PP.H4.abf", "") if key == "signal_coloc" else row.get("PP.H4", "")
                    try: best[locus] = max(best[locus], float(raw))
                    except (ValueError, TypeError): classes["INVALID_H4"] += 1
            report["sources"][key] = {"status": "READ", "path": str(path), "rows": n,
                "columns": fields, "by_locus": dict(sorted(by_locus.items())),
                "harmonization_or_error_counts": dict(classes)}
            for loc, count in by_locus.items():
                target = report["loci"].setdefault(loc, {})
                target[key + "_rows"] = count
                if key in ("signal_coloc", "abf_coloc"):
                    target[key + "_max_h4_observed"] = best[loc]
                    target[key + "_evidence_status"] = "PROVISIONAL_LD_UNVERIFIED"
    report["gates"] = {
        "variant_canonical_qc": "BLOCK_MISSING_BUILD_AND_SE",
        "ancestry_ld": "UNVERIFIED_REFERENCE_METADATA",
        "signal_coloc": "PROVISIONAL_NOT_RECOMPUTED",
        "functional_annotation": "AWAITING_EXPLICIT_CELL_DATASET_AUDIT",
        "integrated_ranking": "NOT_RECOMPUTED"
    }
    return report

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path("/srv/is-analysis/results/is"))
    p.add_argument("--output", type=Path, required=True)
    a=p.parse_args()
    report=audit(a.root)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"sources": {k: {"rows": v.get("rows", 0), "status": v["status"]} for k,v in report["sources"].items()}, "gates": report["gates"], "output": str(a.output)}, ensure_ascii=False))
if __name__ == "__main__": main()
