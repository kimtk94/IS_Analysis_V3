#!/usr/bin/env python3
"""Restricted GitHub issue queue runner for IS research. No arbitrary shell execution."""
import argparse, datetime as dt, json, os, pathlib, subprocess, time, urllib.request, urllib.error, gzip, sys, re
ROOT=pathlib.Path("/srv/is-analysis")
REPO=ROOT/"IS_Analysis_V3"
RESULTS=ROOT/"results/is"
API="https://api.github.com/repos/kimtk94/IS_Analysis_V3"
WORKER_VERSION="20261015-30minute-scheduler-v1"
ALLOWED={"status","inventory","tests","mkdir-output","headers","eas5-audit","candidate-2225-audit","figure-manifest","research-refresh"}
def is_trusted_issue(author, operation, title):
    """Bot may only enqueue restricted scheduled research refresh jobs."""
    if author.lower() == "kimtk94":
        return True
    if author.lower() != "github-actions[bot]" or operation != "research-refresh":
        return False
    title = title or ""
    legacy_daily = r"IS Harness daily research refresh \d{4}-\d{2}-\d{2} KST"
    interval = (
        r"IS Harness scheduled research refresh "
        r"\d{4}-\d{2}-\d{2} (?:[01]\d|2[0-3]):(?:17|47) KST"
    )
    return bool(re.fullmatch(legacy_daily, title) or re.fullmatch(interval, title))

def api(method,path,token,data=None):
    body=None if data is None else json.dumps(data).encode()
    req=urllib.request.Request(API+path,data=body,method=method,headers={"Authorization":"Bearer "+token,"Accept":"application/vnd.github+json","X-GitHub-Api-Version":"2022-11-28","User-Agent":"is-research-queue"})
    with urllib.request.urlopen(req,timeout=20) as response:
        payload=response.read()
        return json.loads(payload) if payload else None
def run_job(op):
    if op in {"eas5-audit", "candidate-2225-audit", "figure-manifest", "research-refresh"}:
        script=pathlib.Path(__file__).resolve().with_name("harness_research_jobs.py")
        if not script.is_file():
            raise RuntimeError("Harness research helper missing; install companion file first")
        proc=subprocess.run([sys.executable, str(script), "--operation", op],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, timeout=260, check=False, cwd=script.parent)
        if proc.returncode:
            raise RuntimeError(f"Research audit failed (exit {proc.returncode}): "+proc.stdout[-3500:])
        return proc.stdout[-9000:]
    if op=="status":
        commands=[["git","-C",str(REPO),"status","--short"],["git","-C",str(REPO),"branch","--show-current"]]
    elif op=="inventory":
        commands=[["find",str(RESULTS/"audits"),"-maxdepth","3","-type","f","-name","IS_EAS5*"]]
    elif op=="tests":
        commands=[["python3","-m","pytest","-q","tests/test_is_public_nominal_qtl_readiness.py","tests/test_is_public_nominal_qtl_region_queue.py","tests/test_is_gwas_qtl_snp_harmonization.py","tests/test_is_inventory_snp_inputs.py"]]
    elif op=="headers":
        base=RESULTS/"audits"
        targets=[
            ("is_eas5_qtl_gwas_official_alt_harmonization_20261011_v1","IS_EAS5_OFFICIAL_QTL_ALT_GWAS_HARMONIZED_VARIANTS.tsv"),
            ("is_eas5_qtl_gwas_official_alt_harmonization_20261011_v1","IS_EAS5_GENE_TISSUE_ALLELE_HARMONIZED_READINESS.tsv"),
            ("is_eas5_qtl_gwas_snp_ref37_20261011_v2","IS_EAS5_GENE_TISSUE_GWAS_SNP_OVERLAP_QC.tsv"),
            ("is_eas5_98_egene_fullcis_gwas_overlap_20261011_v3_fastaanchored","IS_EAS5_98_FDR_EGENE_FULLCIS_BBJ_EAS_MATCHED_VARIANTS.tsv.gz"),
        ]
        lines=[]
        for directory,name in targets:
            p=base/directory/name
            if not p.is_file():
                lines.append(f"MISSING\t{directory}/{name}")
                continue
            opener=gzip.open if p.suffix==".gz" else open
            with opener(p,"rt",encoding="utf-8",errors="replace") as stream:
                header=stream.readline().rstrip("\n\r")
            lines.append(f"FILE\t{directory}/{name}\\tbytes={p.stat().st_size}\\tcolumns={len(header.split(chr(9)))}")
            lines.append(f"HEADER\t{header[:3000]}")
        return "\n".join(lines)
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
    print(f"worker version: {WORKER_VERSION}; allowed: {sorted(ALLOWED)}",flush=True)
    token=os.environ.get("IS_QUEUE_GITHUB_TOKEN","")
    if not token: p.error("IS_QUEUE_GITHUB_TOKEN missing")
    while True:
        try:
            issues=api("GET","/issues?state=open&labels=is-harness-queued&per_page=30",token)
            print(f"queue poll: {len(issues)} queued issue(s)", flush=True)
            for issue in issues:
                if "pull_request" in issue: continue
                num=issue["number"]
                # Only the exact scheduled GitHub Actions bot title can request research-refresh.
                author=issue.get("user",{}).get("login","")
                body=issue.get("body") or ""
                try:
                    spec=json.loads(body)
                    op=spec["operation"]
                    if op not in ALLOWED or set(spec)!={"operation"}:
                        print(f"issue #{num}: skipped (unsupported request, operation={op!r}, fields={sorted(spec)})",flush=True)
                        continue
                except (ValueError,TypeError,KeyError):
                    print(f"issue #{num}: skipped (invalid JSON)",flush=True)
                    continue
                if not is_trusted_issue(author, op, issue.get("title","")):
                    print(f"issue #{num}: skipped (untrusted author for {op})",flush=True)
                    continue
                # Claim by removing queued label BEFORE executing. If claim fails, do not execute.
                try: api("DELETE",f"/issues/{num}/labels/is-harness-queued",token)
                except urllib.error.HTTPError as e:
                    detail=e.read().decode("utf-8",errors="replace")[:700]
                    print(f"issue #{num}: claim failed HTTP {e.code}: {detail}",flush=True)
                    if args.once: raise SystemExit(1) from e
                    continue
                api("POST",f"/issues/{num}/labels",token,{"labels":["is-harness-running"]})
                started_at=dt.datetime.now(dt.timezone.utc)
                start_iso=started_at.isoformat()
                print(f"issue #{num}: START {start_iso}",flush=True)
                api("POST",f"/issues/{num}/comments",token,{"body":f"IS HARNESS START | operation={op} | started_at={start_iso}"})
                try:
                    result=run_job(op)
                    status="PASS"
                except Exception as e:
                    result=str(e)
                    status="FAIL"
                ended_at=dt.datetime.now(dt.timezone.utc)
                end_iso=ended_at.isoformat()
                duration_s=(ended_at-started_at).total_seconds()
                comment=f"IS HARNESS {status} | operation={op} | started_at={start_iso} | ended_at={end_iso} | duration_seconds={duration_s:.3f} | UTC={end_iso}\n\n```text\n{result[:12000]}\n```"
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
