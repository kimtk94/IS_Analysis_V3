#!/usr/bin/env python3
"""CI-only Harness status. No server dependency, no causal inference."""
import datetime as dt,json,os,pathlib,subprocess,sys
def main():
 out=pathlib.Path(os.environ.get("GITHUB_WORKSPACE","."))/"ci_harness_status.json"
 jobs=[
  ("harness_regression", [sys.executable,"-m","unittest","discover","-s","tests","-p","test_is_harness*.py","-v"]),
  ("auto_next_regression", [sys.executable,"-m","unittest","discover","-s","tests","-p","test_is_auto_next*.py","-v"]),
  ("official_gwas_source_regression", [sys.executable,"-m","unittest","discover","-s","tests","-p","test_is_compare_official_gwas_source.py","-v"])
 ]
 items=[]
 for name,cmd in jobs:
  p=subprocess.run(cmd,capture_output=True,text=True,timeout=300)
  combined=(p.stdout+"\n"+p.stderr)
  items.append({"job":name,"status":"PASS" if p.returncode==0 else "FAIL","returncode":p.returncode,
               "summary_tail":combined[-1600:]})
 state={"schema":"is-harness-ci-only-v1","time_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "mode":"CI_ONLY_NO_SSH","source_commit":os.environ.get("GITHUB_SHA",""),
        "ci_jobs":items,"remote_scientific_analysis":"BLOCKED_REMOTE_INPUTS_NOT_AVAILABLE",
        "p1_to_p5_last_remote_status":"HISTORICAL_PASS_NOT_REEXECUTED",
        "independent_replication":"NOT_ESTABLISHED",
        "no_live_operations":True,
        "completed":all(x["status"]=="PASS" for x in items)}
 out.write_text(json.dumps(state,ensure_ascii=False,indent=2)+"\n")
 print("CI_HARNESS_SUMMARY="+json.dumps({k:v for k,v in state.items() if k!="ci_jobs"},sort_keys=True))
 for j in items:print("CI_HARNESS_JOB="+j["job"]+" "+j["status"],flush=True)
 if not state["completed"]:raise SystemExit(1)
if __name__=="__main__":main()
