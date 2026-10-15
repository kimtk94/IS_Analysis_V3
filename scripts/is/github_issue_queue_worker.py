#!/usr/bin/env python3
"""Restricted GitHub issue queue runner for IS research. No arbitrary shell execution."""
import argparse, datetime as dt, json, os, pathlib, subprocess, time, urllib.request, urllib.error
ROOT=pathlib.Path("/srv/is-analysis")
REPO=ROOT/"IS_Analysis_V3"
RESULTS=ROOT/"results/is"
API="https://api.github.com/repos/kimtk94/IS_Analysis_V3"
ALLOWED={"status","inventory","tests","mkdir-output"}
def api(method,path,token,data=None):
    body=None if data is None else json.dumps(data).encode()
    req=urllib.request.Request(API+path,data=body,method=method,headers={"Authorization":"Bearer "+token,"Accept":"application/vnd.github+json","X-GitHub-Api-Version":"2022-11-28","User-Agent":"is-research-queue"})
    with urllib.request.urlopen(req,timeout=20) as response:
        payload=response.read()
        return json.loads(payload) if payload else None
def run_job(op):
    if op=="status":
        commands=[["git","-C",str(REPO),"status","--short"],["git","-C",str(REPO),"branch","--show-current"]]
    elif op=="inventory":
        commands=[["find",str(RESULTS/"audits"),"-maxdepth","3","-type","f","-name","IS_EAS5*"]]
    elif op=="tests":
        commands=[["python3","-m","pytest","-q","tests/test_is_public_nominal_qtl_readiness.py","tests/test_is_public_nominal_qtl_region_queue.py","tests/test_is_gwas_qtl_snp_harmonization.py","tests/test_is_inventory_snp_inputs.py"]]
    elif op=="mkdir-output":
        (RESULTS/"github_actions_runs").mkdir(parents=True,exist_ok=True)
        return "WRITE_ACCESS_CONFIRMED"
    else: raise ValueError("Operation not allowed")
    output=[]
    for command in commands:
        p=subprocess.run(command,cwd=REPO,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240,check=False)
        output.append("$ "+" ".join(command)+"\n"+p.stdout[-10000:])
        if p.returncode: raise RuntimeError("Command failed (exit "+str(p.returncode)+")\n"+"\n".join(output))
    return "\n".join(output)
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--once",action="store_true")
    p.add_argument("--interval",type=int,default=60)
    args=p.parse_args()
    if args.interval<30: p.error("interval must be >= 30")
    token=os.environ.get("IS_QUEUE_GITHUB_TOKEN","")
    if not token: p.error("IS_QUEUE_GITHUB_TOKEN missing")
    while True:
        try:
            issues=api("GET","/issues?state=open&labels=is-harness-queued&per_page=30",token)
            print(f"queue poll: {len(issues)} queued issue(s)", flush=True)
            for issue in issues:
                if "pull_request" in issue: continue
                num=issue["number"]
                # Require trusted author; avoid letting arbitrary issue authors execute server operations.
                author=issue.get("user",{}).get("login","")
                if author.lower()!="kimtk94":
                    print(f"issue #{num}: skipped (untrusted author)",flush=True)
                    continue
                body=issue.get("body") or ""
                try:
                    spec=json.loads(body)
                    op=spec["operation"]
                    if op not in ALLOWED or set(spec)!={"operation"}:
                        print(f"issue #{num}: skipped (unsupported request)",flush=True)
                        continue
                except (ValueError,TypeError,KeyError):
                    print(f"issue #{num}: skipped (invalid JSON)",flush=True)
                    continue
                # Claim by removing queued label BEFORE executing. If claim fails, do not execute.
                try: api("DELETE",f"/issues/{num}/labels/is-harness-queued",token)
                except urllib.error.HTTPError as e:
                    print(f"issue #{num}: claim failed HTTP {e.code}",flush=True)
                    continue
                api("POST",f"/issues/{num}/labels",token,{"labels":["is-harness-running"]})
                try:
                    result=run_job(op)
                    status="PASS"
                except Exception as e:
                    result=str(e)
                    status="FAIL"
                comment=f"IS HARNESS {status} | operation={op} | UTC={dt.datetime.now(dt.timezone.utc).isoformat()}\n\n```text\n{result[:12000]}\n```"
                api("POST",f"/issues/{num}/comments",token,{"body":comment})
                try: api("DELETE",f"/issues/{num}/labels/is-harness-running",token)
                except urllib.error.HTTPError: pass
                api("POST",f"/issues/{num}/labels",token,{"labels":["is-harness-done" if status=="PASS" else "is-harness-failed"]})
                print(f"issue #{num}: {status}",flush=True)
        except Exception as e:
            print("queue poll error:",type(e).__name__,str(e)[:300],flush=True)
            if args.once: raise SystemExit(1) from e
        if args.once: break
        time.sleep(args.interval)
if __name__=="__main__": main()
