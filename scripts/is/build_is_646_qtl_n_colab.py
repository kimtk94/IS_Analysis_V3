#!/usr/bin/env python3
"""Build new one-cell full-646 QTL scalar N stress-test Colab from audited source.

Changes only original Drive-mounted 646 notebook clone and embedded R source;
never overwrites the full original coloc replay notebook or its result folder.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile

def build(template:Path,source_r:Path,dest:Path):
    if dest.exists():raise FileExistsError("DO_NOT_OVERWRITE_EXISTING_NOTEBOOK")
    base=json.loads(template.read_text(encoding="utf8"))
    if len(base["cells"])!=1 or base["cells"][0]["cell_type"]!="code":
        raise ValueError("TEMPLATE_NOT_ONE_CELL")
    py="".join(base["cells"][0]["source"])
    r_code=source_r.read_text(encoding="utf8")
    if r_code.count("IS_QTL_N_SCALAR_SENSITIVITY")!=1:
        raise ValueError("UNEXPECTED_R_SOURCE_SIGNATURE")
    with tempfile.TemporaryDirectory() as tmp:
        r_file=Path(tmp)/"test.R"
        r_file.write_text(r_code,encoding="utf8")
        q=subprocess.run(["Rscript","--vanilla","-e",
           "invisible(parse(file=commandArgs(trailingOnly=TRUE)[[1]]));cat('QTL_N_R_PARSE_OK\\n')",
           str(r_file)],capture_output=True,text=True,timeout=25)
        if q.returncode:raise RuntimeError(q.stderr)
    # Exactly one embedded R string assignment, on one physical Python source line.
    lines=py.splitlines(keepends=True)
    matches=[i for i,l in enumerate(lines) if l.startswith("R_CODE = ")]
    if len(matches)!=1:raise ValueError("TEMPLATE_EMBED_SOURCE_COUNT_CHANGED")
    lines[matches[0]]="R_CODE = "+repr(r_code)+"\n"
    py="".join(lines)
    replacements={
        'WORK = Path("/content/is_phase10_11_646_replay")':
        'WORK = Path("/content/is_phase10_11_646_qtl_n_stress")',
        'OUTPUT = GDRIVE/"results/IS_PHASE10_11_LEGACY_646_SNP_REPLAY_V1"':
        'OUTPUT = GDRIVE/"results/IS_PHASE10_11_646_QTL_N_SCALAR_SENSITIVITY_V1"',
        'RUNLOG = WORK/"IS_PHASE10_11_646_R_LOG.txt"':
        'RUNLOG = WORK/"IS_PHASE10_11_646_QTL_N_R_LOG.txt"',
        'script=WORK/"IS_PHASE10_11_646_batch.R"':
        'script=WORK/"IS_PHASE10_11_646_qtl_n.R"',
        "ORIGINAL_646_R_PARSE_OK":"QTL_N_SENSITIVITY_R_PARSE_OK",
        'status=OUTPUT/"IS_LEGACY_646_SNP_REPLAY_STATUS.tsv"':
        'status=OUTPUT/"IS_646_QTL_N_SCALAR_SENSITIVITY.tsv"',
        '"IS_646_SOURCE_FILE_SHA256.tsv"':'"IS_646_QTL_N_SOURCE_SHA256.tsv"',
        '"schema":"IS_PHASE10_11_646_ORIGINAL_SNP_BATCH_COLAB_V1"':
        '"schema":"IS_PHASE10_11_646_QTL_N_SCALAR_COLAB_V1"',
        '"status":"FULL_646_NUMERIC_REPLAY" if counts["PASS"]==646 else "PARTIAL_OR_FAILED"':
        '"status":"FULL_646_SCALAR_N_STRESS_COMPLETE" if counts["PASS"]==646 else "PARTIAL_OR_FAILED"',
        '"evidence_boundary":"NUMERIC_REPRODUCIBILITY_NOT_CAUSAL_GENE_OR_ANCESTRY_LD"':
        '"evidence_boundary":"N_MIN_MEDIAN_FIRST_SCALAR_STRESS_NOT_PER_SNP_N_MODEL_OR_CAUSAL_GENE"',
        '"IS_646_COLAB_RUN_MANIFEST.json"':'"IS_646_N_SCALAR_COLAB_MANIFEST.json"',
        '"IS_646_FULL_RUN_LOG.txt"':'"IS_646_N_SCALAR_FULL_RUN_LOG.txt"',
        '"IS_PHASE10_11_646_BATCH=RUN_ENDED_AND_STATUS_ACCOUNTED"':
        '"IS_PHASE10_11_646_QTL_N_STRESS=RUN_ENDED_AND_STATUS_ACCOUNTED"',
    }
    for old,new in replacements.items():
        count=py.count(old)
        if count!=1:raise ValueError(f"ORIGINAL_TEMPLATE_SOURCE_SIGNATURE_COUNT_{count}:{old}")
        py=py.replace(old,new)
    # Immutable metadata contract: prevent using altered historical source
    needle='print("MASTER_SHA256=",sha(MASTER))'
    check=('assert sha(INDEX)=="1fe044f1ffd8e9aa6b678718b6772063bd11edd462d7bede4104843cdb0c7414", "INDEX_SOURCE_SHA_CHANGED"\n'
           'assert sha(MASTER)=="0160a01fc6433f9bec839113ac3a2390927fc07db9252d014b0cc4b6b11716ea", "MASTER_SOURCE_SHA_CHANGED"')
    if py.count(needle)!=1:raise ValueError("MISSING_FROZEN_SOURCE_HASH_ANCHOR")
    py=py.replace(needle,needle+"\n"+check)
    py=py.replace("One-cell Phase10/11 original 646-pair SNP-level coloc replay.",
                  "One-cell Phase10/11 original 646-assay QTL scalar N sensitivity (min/median/first).")
    py=py.replace("One-cell", "One-cell",1)
    # Avoid declaring full N sensitivity if original index has fewer files
    needle='counts={k:sum(row["status"]==k for row in results) for k in ("PASS","FAIL","MISSING_INPUT")}'
    if py.count(needle)!=1:raise ValueError("STATUS_CHECK_ANCHOR_MISSING")
    py=py.replace(needle,needle+
       '\nprint("N_SCALAR_646_SUMMARY_FILE=",OUTPUT/"IS_646_QTL_N_SCALAR_SUMMARY.tsv")')
    # Bind the source-only G1 audit into the one-cell workflow. It never
    # independently establishes GWAS effect-allele orientation.
    g1_source=(Path(__file__).parent/"audit_is_646_g1_source_alleles.py").read_text(encoding="utf8")
    if g1_source.count("IS_PHASE10_11_G1_SOURCE_ALLELE_STRUCTURAL_V1")!=1:
        raise ValueError("MISSING_G1_SOURCE_SIGNATURE")
    needle='manifest={'
    if py.count(needle)!=1:raise ValueError("MANIFEST_ANCHOR_CHANGED")
    g1_snippet="G1_CODE = "+repr(g1_source)+"\n"+"""G1_SCRIPT=WORK/"audit_is_g1_source_alleles.py"
G1_SCRIPT.write_text(G1_CODE,encoding="utf8")
G1_OUTPUT=OUTPUT/"G1_ORIGINAL_ALLELE_STRUCTURAL_V1"
g1_proc=subprocess.run(["python3",str(G1_SCRIPT),
     "--index",str(INDEX),"--inputs",str(INPUT),"--out",str(G1_OUTPUT)],
     capture_output=True,text=True,check=False)
print("G1_STRUCTURAL_AUDIT_PROCESS_EXIT=",g1_proc.returncode)
print(g1_proc.stdout[-4000:])
if g1_proc.returncode: raise RuntimeError(g1_proc.stderr[-5000:])
G1_SUMMARY=json.loads((G1_OUTPUT/"IS_G1_SOURCE_STRUCTURAL_SUMMARY.json").read_text(encoding="utf8"))
G1_COUNTS=G1_SUMMARY["status_counts"]
G1_STATUS=("G1_STRUCTURAL_646_PASS_EXTERNAL_EFFECT_ALLELE_UNVERIFIED" if G1_COUNTS=={"PASS_STRUCTURAL_ONLY":646} else "G1_STRUCTURAL_INCOMPLETE_OR_FAIL")
print("G1_SOURCE_QC_STATUS=",G1_STATUS)
"""
    py=py.replace(needle,g1_snippet+needle)
    # Case-specific manifest extension, not original prior/posterior replication manifest.
    needle='    "R_package_coloc":"5.2.3",'
    if py.count(needle)!=1:raise ValueError("METADATA_SIGNATURE_CHANGED")
    py=py.replace(needle,needle+
      '\n    "sensitivity_scenarios":["ORIGINAL_FIRST_EQ_MAX","MEDIAN_SNP_N","MIN_SNP_N"],\n'
      '    "n_assumption_status":"THREE_SCALAR_N_STRESS_NOT_SNP_SPECIFIC_MODEL",')
    needle='    "n_assumption_status":"THREE_SCALAR_N_STRESS_NOT_SNP_SPECIFIC_MODEL",'
    if py.count(needle)!=1:raise ValueError("N_SCENARIO_MANIFEST_MISSING")
    py=py.replace(needle,needle+
       '\n    "g1_source_audit_status":G1_STATUS,'
       '\n    "g1_source_audit_counts":G1_COUNTS,'
       '\n    "g1_script_sha256":sha(G1_SCRIPT),')
    ast.parse(py)
    book={
        "cells":[{"cell_type":"code","execution_count":None,"metadata":{},
                  "outputs":[],"source":py.splitlines(keepends=True)}],
        "metadata":{"colab":{"provenance":[],"name":dest.name},
                    "kernelspec":{"display_name":"Python 3","name":"python3"},
                    "language_info":{"name":"python"}},
        "nbformat":4,"nbformat_minor":0
    }
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(book,ensure_ascii=False,indent=1)+"\n",encoding="utf8")
    if "".join(json.loads(dest.read_text())["cells"][0]["source"])!=py:
        raise ValueError("NOTEBOOK_WRITE_READBACK_MISMATCH")
    print(json.dumps({"notebook":str(dest),"bytes":dest.stat().st_size,
                      "cells":1,
                      "R_SHA256":hashlib.sha256(r_code.encode()).hexdigest(),
                      "notebook_SHA256":hashlib.sha256(dest.read_bytes()).hexdigest()},indent=2))
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--template",type=Path,required=True)
    p.add_argument("--source-r",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    build(a.template,a.source_r,a.out)
if __name__=="__main__":main()
