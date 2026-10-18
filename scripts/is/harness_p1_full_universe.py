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
 if not hits:raise ValueError("P1_SOURCE_ROOT_UNRESOLVED count=0")
 candidates=[]
 for directory in hits:
  paths={n:directory/n for n in NAMES}
  try:
   m=json.loads(paths[NAMES[2]].read_text())
   hashes={n:digest(p) for n,p in paths.items()}
   valid=(m.get("stable_gene_ids")==2225 and m.get("gene_region_pairs")==2425)
   candidates.append((directory,hashes,valid))
  except (OSError,ValueError) as exc:
   print("P1_CANDIDATE_INVALID="+str(directory)+" "+str(type(exc).__name__),flush=True)
   candidates.append((directory,None,False))
 for directory,hashes,valid in candidates:
  print("P1_CANDIDATE="+json.dumps(dict(path=str(directory),sha256=hashes,valid_manifest=valid),sort_keys=True),flush=True)
 valid=[item for item in candidates if item[2]]
 if len(valid)==2 and len(valid)==len(candidates):
  from collections import Counter
  a,b=valid[0][0],valid[1][0]
  def tsv_profile(path,key):
   with path.open(newline="") as stream:
    rd=csv.DictReader(stream,delimiter="\t")
    columns=rd.fieldnames or []
    if key not in columns:raise ValueError("MISSING_ID_COLUMN "+key)
    counts=Counter()
    rows={}
    for row in rd:
     ident=row[key]
     counts[ident]+=1
     rows.setdefault(ident,[]).append(row)
    return columns,counts,rows
  try:
   for name,key in [(NAMES[0],"gene_id_stable"),(NAMES[1],"gene_id_stable")]:
    ca,na,ra=tsv_profile(a/name,key)
    cb,nb,rb=tsv_profile(b/name,key)
    shared=ra.keys() & rb.keys()
    changed=sum(ra[k]!=rb[k] for k in shared)
    summary=dict(file=name,left=str(a),right=str(b),
      left_rows=sum(na.values()),right_rows=sum(nb.values()),
      left_unique_ids=len(ra),right_unique_ids=len(rb),
      only_left_ids=len(ra.keys()-rb.keys()),only_right_ids=len(rb.keys()-ra.keys()),
      changed_shared_ids=changed,left_duplicate_ids=sum(v>1 for v in na.values()),
      right_duplicate_ids=sum(v>1 for v in nb.values()),
      column_schema_equal=(ca==cb),
      common_changed_columns=sorted({col for k in shared for r1,r2 in zip(ra[k],rb[k])
        for col in set(r1)&set(r2) if r1[col]!=r2[col]})[:60])
    print("P1_SOURCE_DIFF="+json.dumps(summary,sort_keys=True),flush=True)
   ma=json.loads((a/NAMES[2]).read_text())
   mb=json.loads((b/NAMES[2]).read_text())
   print("P1_MANIFEST_DIFF="+json.dumps(dict(
    only_left_keys=sorted(ma.keys()-mb.keys()),
    only_right_keys=sorted(mb.keys()-ma.keys()),
    changed_shared_keys=sorted(k for k in ma.keys()&mb.keys() if ma[k]!=mb[k])
   ),sort_keys=True),flush=True)
  except Exception as exc:
   print("P1_SOURCE_COMPARISON_ERROR="+type(exc).__name__+": "+str(exc),flush=True)
 signatures={tuple(item[1][n] for n in NAMES) for item in valid}
 if len(valid)!=len(candidates) or len(signatures)!=1:
  raise ValueError("P1_SOURCE_ROOT_CONFLICT candidates="+str(len(candidates))+" valid="+str(len(valid))+" distinct_hash_sets="+str(len(signatures)))
 src=sorted((item[0] for item in valid),key=lambda p:(len(p.parts),str(p)))[0]
 files={n:src/n for n in NAMES}
 sha={n:digest(p) for n,p in files.items()}
 print("P1_EQUIVALENT_SOURCE_COPIES="+str(len(valid)),flush=True)
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
