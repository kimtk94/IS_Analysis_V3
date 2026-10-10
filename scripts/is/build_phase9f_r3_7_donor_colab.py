#!/usr/bin/env python3
"""Build a self-contained, one-cell R3_7 Colab notebook from verified R3_6.

No source RDS is loaded during generation; embedded R syntax is parsed with
the installed Rscript. The original notebook is never overwritten.
"""
import argparse
import ast
import json
from pathlib import Path
import shutil
import subprocess
import tempfile


def replace_once(s, a, b):
    if s.count(a) != 1:
        raise ValueError(f"Notebook signature count mismatch: {s.count(a)} for {a[:75]!r}")
    return s.replace(a, b, 1)


def extract_embedded_r(py_src):
    tree = ast.parse(py_src)
    for node in tree.body:
        if isinstance(node, ast.Assign):
            if any(isinstance(t, ast.Name) and t.id == "r_code" for t in node.targets):
                r = ast.literal_eval(node.value)
                if not isinstance(r, str):
                    raise TypeError("r_code must be a literal string")
                return r
    raise ValueError("Missing Python literal r_code")


def build_notebook(base, donor_helper, output, rscript="Rscript"):
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite notebook: {output}")
    nb = json.loads(base.read_text(encoding="utf8"))
    cells = nb.get("cells", [])
    if len(cells) != 1 or cells[0]["cell_type"] != "code":
        raise ValueError("Expected exact one-cell Colab notebook")
    code = "".join(cells[0]["source"])
    donor_code = donor_helper.read_text(encoding="utf8")
    if "is_phase9f_donor_pseudobulk <- function" not in donor_code:
        raise ValueError("Wrong donor R helper")

    for a, b in (
        ("R3_6_FGF5_FEATURE_QC_20261010", "R3_7_DONOR_PSEUDOBULK_20261010"),
        ("R3_6_FGF5_FEATURE_QC_PARSEVERIFIED", "R3_7_DONOR_PSEUDOBULK_PARSEVERIFIED"),
    ):
        code = replace_once(code, a, b)
    # R log path occurs in both local and Drive definitions.
    log_src = "R3_6_FGF5_ANALYSIS.log"
    if code.count(log_src) != 2:
        raise ValueError("Unexpected R3_6 analysis log reference count")
    code = code.replace(log_src, "R3_7_DONOR_PSEUDOBULK_ANALYSIS.log")

    # Add the helper source as literal R code at the beginning of the R program,
    # which makes Colab self-contained even when GitHub is inaccessible.
    code = replace_once(
        code,
        "r_code = r'''\nsuppressPackageStartupMessages({",
        "r_code = r'''\n" + donor_code + "\nsuppressPackageStartupMessages({"
    )

    # One count layer is expected and verified in GSE256493. Explicitly fail
    # for a different layer layout instead of double-counting donor/cell combos.
    count_anchor = '''if (length(count_layers)==0) {
  stop("No raw count layer found; normalized data will not be substituted for pseudobulk.")
}'''
    code = replace_once(
        code, count_anchor,
        count_anchor + '''
if (length(count_layers)!=1L) {
  stop("R3_7 donor aggregation is validated for exactly one RNA count layer; separate layer harmonization required")
}'''
    )
    call = '''  donor_result <- is_phase9f_donor_pseudobulk(
    M=M,annotation=layer_annotation,patient=as.character(md[idx_md,"Patient"]),
    library_size=as.numeric(lib_size),genes=targets,
    resolve_feature=resolve_target_hits,min_cells=20L,min_paired_donors=4L
  )
  donor_comparisons <- is_phase9f_paired_comparisons(
    donor_result,is_phase9f_predeclared_comparisons(),min_paired_donors=4L
  )
  is_phase9f_write_donor_outputs(donor_result,donor_comparisons,OUT)
  cat("DONOR_PSEUDOBULK_TABLES_WRITTEN\\\\n")
'''
    code = replace_once(
        code,
        "  rm(M)\n  gc(verbose=FALSE)\n}\n\npresence <-",
        call + "\n  rm(M)\n  gc(verbose=FALSE)\n}\n\npresence <-"
    )
    code = replace_once(
        code,
        'required = [\n    "FGF5_FEATURE_DIAGNOSTIC.tsv",',
        '''required = [
    "HUMAN_DONOR_PSEUDOBULK.tsv",
    "HUMAN_DONOR_CELLTYPE_COVERAGE.tsv",
    "HUMAN_DONOR_FEATURE_STATUS.tsv",
    "HUMAN_DONOR_PAIRED_COMPARISONS.tsv",
    "HUMAN_DONOR_PSEUDOBULK_README.md",
    "FGF5_FEATURE_DIAGNOSTIC.tsv",'''
    )
    code = replace_once(
        code,
        '    "required_r_outputs": required,',
        '''    "required_r_outputs": required,
    "donor_analysis_unit": "Patient",
    "min_cells_per_donor_celltype": 20,
    "min_paired_donors": 4,
    "donor_inference": "DESCRIPTIVE_PSEUDOBULK_NO_DE_NO_STROKE_DISEASE_CONTRAST",'''
    )
    code = replace_once(
        code,
        'print("PHASE9F_E_HUMAN=COMPLETE")',
        '''paired = pd.read_csv(OUTPUT_DIR/"HUMAN_DONOR_PAIRED_COMPARISONS.tsv", sep="\\t")
coverage = pd.read_csv(OUTPUT_DIR/"HUMAN_DONOR_CELLTYPE_COVERAGE.tsv", sep="\\t")
display(paired)
display(coverage)
print("R3_7_DONOR_PSEUDOBULK=DESCRIPTIVE_COMPLETE")
print("PHASE9F_E_HUMAN=COMPLETE")'''
    )

    # Fail if the Python one-cell source or COMPLETE embedded R script is invalid.
    ast.parse(code)
    full_r = extract_embedded_r(code)
    if rscript:
        exe = shutil.which(rscript)
        if exe is None:
            raise FileNotFoundError(f"Rscript executable not found: {rscript}")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "R3_7_embedded.R"
            path.write_text(full_r, encoding="utf8")
            probe = subprocess.run(
                [exe, "-e",
                 "invisible(parse(file=commandArgs(trailingOnly=TRUE)[[1]]));cat('R3_7_R_PARSE_OK\\n')",
                 str(path)],capture_output=True,text=True,timeout=20,
            )
            if probe.returncode != 0:
                raise SyntaxError(f"Embedded R parse failure:\n{probe.stderr}")
            print(probe.stdout.strip())

    cells[0]["source"] = code.splitlines(keepends=True)
    cells[0]["execution_count"] = None
    cells[0]["outputs"] = []
    nb.setdefault("metadata", {}).setdefault("colab", {})["name"] = output.name
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(nb, ensure_ascii=False, indent=1) + "\n",
                      encoding="utf8")
    # Readback test guards against lost newlines in the single cell.
    if "".join(json.loads(output.read_text(encoding="utf8"))["cells"][0]["source"]) != code:
        raise ValueError("Notebook readback mismatch")
    return {"notebook": str(output), "bytes": output.stat().st_size,
            "cells": 1, "embedded_r_lines": len(full_r.splitlines()),
            "donor_rules_embedded": True, "preload_r_parse": "PASS"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base",required=True,type=Path)
    ap.add_argument("--donor-r",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    ap.add_argument("--rscript",default="Rscript")
    a = ap.parse_args()
    print(json.dumps(build_notebook(a.base,a.donor_r,a.output,a.rscript),indent=2))


if __name__=="__main__":
    main()
