#!/usr/bin/env python3
import argparse,csv,hashlib,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is")
NAMES=["IS_2225_GENE_EVIDENCE_MATRIX.tsv","IS_2425_GENE_REGION_EVIDENCE.tsv","IS_GENE_EVIDENCE_MATRIX_MANIFEST.json"]
def digest(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for chunk in iter(lambda:f.read(1048576),b""):h.update(chunk)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser()
 for k in ["mode","run_id","attempt","source_commit","audit_script"]:p.add_argument("--"+k.replace("_","-"))
 a=p.parse_args()
 if a.mode not in ("probe","run"):raise ValueError("invalid mode")
 hits=[]
 for root,dirs,files in os.walk(ROOT):
  dirs[:]=sorted(x for x in dirs if not x.startswith(".") and not Path(root,x).is_symlink())
  if all(n in files for n in NAMES):hits.append(Path(root))
  if len(hits)>1:break
 if len(hits)!=1:raise ValueError("P1_SOURCE_ROOT_UNRESOLVED count="+str(len(hits)))
 src=hits[0]
 files={n:src/n for n in NAMES}
 m=json.loads(files[NAMES[2]].read_text())
 if m.get("stable_gene_ids")!=2225 or m.get("gene_region_pairs")!=2425:raise ValueError("P1_CARDINALITY_FAIL")
 sha={n:digest(p) for n,p in files.items()}
 print("SSH_AUTH=PASS",flush=True)
 print("P1_SOURCE_ROOT="+str(src),flush=True)
 print("P1_SHA256="+json.dumps(sha,sort_keys=True),flush=True)
 print("P1_SOURCE_INPUTS=PASS",flush=True)
 if a.mode=="probe":
  print("P1_PROBE_STATUS=PASS",flush=True);return
 if not (a.run_id and a.run_id.isdecimal() and a.attempt and a.attempt.isdecimal() and a.source_commit):raise ValueError("PROVENANCE_MISSING")
 out=ROOT/"audits"/("is_p1_full_universe_gh_"+a.run_id+"_attempt_"+a.attempt)
 if out.exists():raise FileExistsError("IMMUTABLE_OUTPUT_EXISTS")
 start=time.time()
 cmd=[sys.executable,a.audit_script,"--root",str(src),"--out",str(out)]
 x=subprocess.run(cmd,capture_output=True,text=True)
 if x.returncode:
  print(x.stdout[-3000:]);print(x.stderr[-3000:],file=sys.stderr)
  raise RuntimeError("AUDIT_FAILED exit="+str(x.returncode))
 outputs=["IS_2225_NONCAUSAL_WORKFLOW_PRIORITY.tsv","IS_TOP40_NON_G0022_EXISTING_DATA_QUEUE.tsv","IS_2225_WORKFLOW_PRIORITY_MANIFEST.json"]
 if not all((out/n).is_file() for n in outputs):raise ValueError("P1_OUTPUT_MISSING")
 report=json.loads((out/outputs[2]).read_text())
 if report.get("total_source_unique_genes")!=2225 or report.get("total_source_gene_region_pairs")!=2425:raise ValueError("P1_OUTPUT_COUNT_FAIL")
 for key,name in [("manifest",NAMES[2]),("gene_matrix",NAMES[0]),("region_matrix",NAMES[1])]:
  if report.get("source_"+key+"_sha256")!=sha[name]:raise ValueError("HASH_MISMATCH "+key)
 manifest=dict(status="P1_COMPLETE_NONCAUSAL",start_epoch=start,end_epoch=time.time(),source_commit=a.source_commit,command=cmd,run_id=a.run_id,attempt=a.attempt,inputs=sha,output_hashes={n:digest(out/n) for n in outputs},causal_status="NOT_ESTABLISHED")
 (out/"P1_GITHUB_RUN_PROVENANCE.json").write_text(json.dumps(manifest,indent=2)+"\n")
 print("P1_AUDIT_STATUS=COMPLETE",flush=True)
 print("P1_OUTPUT_DIR="+str(out),flush=True)
if __name__=="__main__":
 try:main()
 except Exception as e:
  print("P1_HARNESS_ERROR="+type(e).__name__+": "+str(e),file=sys.stderr,flush=True)
  sys.exit(1)
