#!/usr/bin/env python3
"""Build ONE Python-cell Colab notebook embedding the frozen 646-input R script.

The resulting notebook reads original TSV inputs on the user's mounted Drive,
thus avoiding per-file Google Drive API requests and the server's 403 limits.
"""
import argparse
import ast
import json
from pathlib import Path
import subprocess
import tempfile


def build(source_r: Path, dest: Path):
    if dest.exists():
        raise FileExistsError("Do not overwrite existing Colab notebook")
    r_code=source_r.read_text(encoding="utf8")
    if "IS_LEGACY_646_BATCH_REPLAY=COMPLETE_INPUT_ACCOUNTING" not in r_code:
        raise ValueError("Unexpected embedded R batch replay source")
    with tempfile.TemporaryDirectory() as directory:
        p=Path(directory)/"parse_check.R"
        p.write_text(r_code,encoding="utf8")
        proc=subprocess.run(["Rscript","--vanilla","-e",
            "invisible(parse(file=commandArgs(trailingOnly=TRUE)[1]));cat('R_BATCH_PARSE_OK\\n')",
            str(p)],capture_output=True,text=True,timeout=25)
        if proc.returncode:raise RuntimeError(proc.stderr)
        print(proc.stdout.strip())
    py=r'''# One-cell Phase10/11 original 646-pair SNP-level coloc replay.
# Input: original archived BBJ x GTEx V8 four-locus screen; NOT 80 new windows.
from pathlib import Path
import csv
import hashlib
import json
import shutil
import subprocess
import time
from google.colab import drive
drive.mount("/content/drive")
WORK = Path("/content/is_phase10_11_646_replay")
WORK.mkdir(parents=True,exist_ok=True)
GDRIVE = Path("/content/drive/MyDrive/MASTER_DEGREE/IS_COLAB")
INDEX = GDRIVE/"coloc_ready_v2/COLOC_ABF_INPUT_INDEX_V2.tsv"
MASTER = GDRIVE/"results_v2/COLOC_ABF_MASTER_V2.tsv"
INPUT = GDRIVE/"coloc_ready_v2/abf"
OUTPUT = GDRIVE/"results/IS_PHASE10_11_LEGACY_646_SNP_REPLAY_V1"
RUNLOG = WORK/"IS_PHASE10_11_646_R_LOG.txt"
for p in [INDEX, MASTER, INPUT]:
    if not p.exists():
        raise FileNotFoundError(f"Missing original Drive source: {p}")
if OUTPUT.exists():
    raise FileExistsError(f"Refusing existing result folder: {OUTPUT}. Use new version suffix after reviewing past results.")
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(8*1024*1024),b""):
            h.update(block)
    return h.hexdigest()
with INDEX.open(encoding="utf8",newline="") as handle:
    source_rows=list(csv.DictReader(handle,delimiter="\t"))
if len(source_rows)!=646:
    raise ValueError(f"Expected original 646 index rows, got {len(source_rows)}")
input_files=[INPUT/Path(x["file"]).name for x in source_rows]
n_available=sum(p.exists() and p.stat().st_size>0 for p in input_files)
print("INPUTS_INDEXED=",len(input_files))
print("INPUTS_AVAILABLE_ON_DRIVE=",n_available)
print("INPUTS_MISSING_ON_DRIVE=",646-n_available)
print("INDEX_SHA256=",sha(INDEX))
print("MASTER_SHA256=",sha(MASTER))
print("MEMORY=",Path("/proc/meminfo").read_text().splitlines()[0])
if shutil.which("Rscript") is None:
    subprocess.run(["apt-get","update","-qq"],check=True)
    subprocess.run(["apt-get","install","-y","-qq",
                    "r-base","r-base-dev","build-essential",
                    "libcurl4-openssl-dev","libssl-dev","libxml2-dev",
                    "r-cran-data.table"],check=True)
R_INSTALL = r"""
repos <- "https://cloud.r-project.org"
if (!requireNamespace("data.table",quietly=TRUE)) {
  install.packages("data.table",repos=repos,Ncpus=2)
}
ok <- requireNamespace("coloc",quietly=TRUE) &&
      as.character(packageVersion("coloc"))=="5.2.3"
if (!ok) {
  deps <- c("ggplot2","susieR","mixsqp","Rcpp","viridis")
  need <- deps[!vapply(deps,requireNamespace,logical(1),quietly=TRUE)]
  if (length(need)>0) install.packages(need,repos=repos,Ncpus=2)
  # Prefer CRAN current exact source; fall back to its archived release.
  candidates <- c(
    "https://cran.r-project.org/src/contrib/coloc_5.2.3.tar.gz",
    "https://cran.r-project.org/src/contrib/Archive/coloc/coloc_5.2.3.tar.gz")
  installed <- FALSE
  for (url in candidates) {
    tryCatch({
      suppressWarnings(install.packages(url,repos=NULL,type="source"))
      installed <- requireNamespace("coloc",quietly=TRUE) &&
        as.character(packageVersion("coloc"))=="5.2.3"
    },error=function(e) cat("PACKAGE_SOURCE_FAILED:",conditionMessage(e),"\n"))
    if (installed) break
  }
  if (!installed) stop("Unable to pin coloc 5.2.3 from CRAN current or archive")
}
stopifnot(requireNamespace("coloc",quietly=TRUE))
stopifnot(requireNamespace("data.table",quietly=TRUE))
stopifnot(as.character(packageVersion("coloc"))=="5.2.3")
cat("COL0C_PINNED_5.2.3_READY\n")
"""
install=WORK/"install_is_646_deps.R"
install.write_text(R_INSTALL,encoding="utf8")
subprocess.run(["Rscript","--vanilla",str(install)],check=True)
R_CODE = __EMBED_R_CODE__
script=WORK/"IS_PHASE10_11_646_batch.R"
script.write_text(R_CODE,encoding="utf8")
with script.open(encoding="utf8") as f:
    print("R_SCRIPT_SHA256=",sha(script))
subprocess.run(["Rscript","--vanilla","-e",
  "invisible(parse(file=commandArgs(trailingOnly=TRUE)[1]));cat('ORIGINAL_646_R_PARSE_OK\\n')",
  str(script)],check=True)
start=time.time()
with RUNLOG.open("w",encoding="utf8") as log:
    proc=subprocess.run(["Rscript","--vanilla",str(script),
          str(INDEX),str(MASTER),str(INPUT),str(OUTPUT)],
          stdout=log,stderr=subprocess.STDOUT,check=False)
print("COLOC_PROCESS_EXIT=",proc.returncode)
print("RUNTIME_SECONDS=",round(time.time()-start,2))
lines=RUNLOG.read_text(encoding="utf8",errors="replace").splitlines()
print("\n".join(lines[-35:]))
if proc.returncode!=0:
    raise RuntimeError(f"R replay failed (exit={proc.returncode}); inspect {RUNLOG}. Original files unchanged.")
status=OUTPUT/"IS_LEGACY_646_SNP_REPLAY_STATUS.tsv"
if not status.exists():
    raise FileNotFoundError("Runner emitted no audit status")
with status.open(encoding="utf8",newline="") as f:
    results=list(csv.DictReader(f,delimiter="\t"))
if len(results)!=646:
    raise ValueError("Final status must include exactly 646 original rows")
counts={k:sum(row["status"]==k for row in results) for k in ("PASS","FAIL","MISSING_INPUT")}
if sum(counts.values())!=646:
    raise ValueError("Status counts do not sum to 646")
# Hash each accessible SNP input after the audit. The output remains a fresh folder.
with (OUTPUT/"IS_646_SOURCE_FILE_SHA256.tsv").open("w",encoding="utf8",newline="") as audit_file:
    writer=csv.writer(audit_file,delimiter="\t")
    writer.writerow(["locus","dataset_key","gene_base","filename","status","source_sha256"])
    for entry in results:
        source_input=INPUT/entry["filename"]
        file_hash=sha(source_input) if source_input.is_file() else "MISSING"
        writer.writerow([entry["locus"],entry["dataset_key"],entry["gene_base"],
                        entry["filename"],entry["status"],file_hash])
manifest={
    "schema":"IS_PHASE10_11_646_ORIGINAL_SNP_BATCH_COLAB_V1",
    "source":"ARCHIVED_FOUR_BBJ_LOCI_GTEX_NOT_EXPANDED_80",
    "original_code_baseline_priors":{"p1":1e-4,"p2":1e-4,"p12":1e-5},
    "R_package_coloc":"5.2.3",
    "inputs_available_before_run":n_available,
    "counts":counts,
    "status":"FULL_646_NUMERIC_REPLAY" if counts["PASS"]==646 else "PARTIAL_OR_FAILED",
    "evidence_boundary":"NUMERIC_REPRODUCIBILITY_NOT_CAUSAL_GENE_OR_ANCESTRY_LD",
    "source_sha256":{"original_index":sha(INDEX),"original_master":sha(MASTER),
                     "embedded_R_script":sha(script),"output_status":sha(status)},
    "runtime_seconds":round(time.time()-start,2)}
(OUTPUT/"IS_646_COLAB_RUN_MANIFEST.json").write_text(
  json.dumps(manifest,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
shutil.copyfile(RUNLOG,OUTPUT/"IS_646_FULL_RUN_LOG.txt")
print(json.dumps(manifest,indent=2,ensure_ascii=False))
print("IS_PHASE10_11_646_BATCH=RUN_ENDED_AND_STATUS_ACCOUNTED")
'''
    py=py.replace("__EMBED_R_CODE__",repr(r_code))
    ast.parse(py)
    notebook={
      "cells":[{"cell_type":"code","execution_count":None,"metadata":{},
                "outputs":[],"source":py.splitlines(keepends=True)}],
      "metadata":{"colab":{"provenance":[],"name":dest.name},
                  "kernelspec":{"display_name":"Python 3","name":"python3"},
                  "language_info":{"name":"python"}},
      "nbformat":4,"nbformat_minor":0
    }
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(notebook,ensure_ascii=False,indent=1)+"\n",encoding="utf8")
    check=json.loads(dest.read_text(encoding="utf8"))
    if len(check["cells"])!=1 or "".join(check["cells"][0]["source"])!=py:
        raise ValueError("Serialized one-cell source does not match")
    print(json.dumps({"notebook":str(dest),"bytes":dest.stat().st_size,
                      "cells":1,"embedded_R_SHA256":__import__("hashlib").sha256(r_code.encode("utf8")).hexdigest()},
                    ensure_ascii=False,indent=2))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-r",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();build(a.source_r,a.out)
if __name__=="__main__":
    main()
