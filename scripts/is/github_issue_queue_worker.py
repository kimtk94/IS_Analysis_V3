#!/usr/bin/env python3
"""Restricted GitHub issue queue runner for IS research. No arbitrary shell execution."""
import argparse, datetime as dt, json, os, pathlib, subprocess, time, urllib.request, urllib.error, gzip, sys, re, importlib.util, tempfile
ROOT=pathlib.Path("/srv/is-analysis")
REPO=ROOT/"IS_Analysis_V3"
RESULTS=ROOT/"results/is"
API="https://api.github.com/repos/kimtk94/IS_Analysis_V3"
WORKER_VERSION="20261016-science-heartbeat-v2"
ALLOWED={"status","inventory","tests","mkdir-output","headers","eas5-audit","candidate-2225-audit","figure-manifest","research-refresh"}
def is_trusted_issue(author, operation, title):
    """Bot may only enqueue restricted scheduled research refresh jobs."""
    if author.lower() == "kimtk94":
        return True
    if author.lower() != "github-actions[bot]" or operation != "research-refresh":
        return False
    title = title or ""
    legacy_daily = r"IS Harness daily research refresh \d{4}-\d{2}-\d{2} KST"
    interval_legacy = (
        r"IS Harness scheduled research refresh "
        r"\d{4}-\d{2}-\d{2} (?:[01]\d|2[0-3]):(?:17|47) KST"
    )
    interval_5m = (
        r"IS Harness scheduled research refresh "
        r"\d{4}-\d{2}-\d{2} (?:[01]\d|2[0-3]):[0-5][05] KST"
    )
    return bool(
        re.fullmatch(legacy_daily, title)
        or re.fullmatch(interval_legacy, title)
        or re.fullmatch(interval_5m, title)
    )

def api(method,path,token,data=None):
    body=None if data is None else json.dumps(data).encode()
    req=urllib.request.Request(API+path,data=body,method=method,headers={"Authorization":"Bearer "+token,"Accept":"application/vnd.github+json","X-GitHub-Api-Version":"2022-11-28","User-Agent":"is-research-queue"})
    with urllib.request.urlopen(req,timeout=20) as response:
        payload=response.read()
        return json.loads(payload) if payload else None
def scheduled_slot_title(now_utc):
    """Latest 5-min KST bucket; a late systemd invocation may still recover it."""
    if now_utc.tzinfo is None:
        raise ValueError("Timezone-aware datetime required")
    slot=now_utc.replace(minute=(now_utc.minute//5)*5,second=0,microsecond=0)
    slot=slot.astimezone(dt.timezone.utc)+dt.timedelta(hours=9)
    return "IS Harness scheduled research refresh "+slot.strftime("%Y-%m-%d %H:%M KST")


def unchanged_inputs_within_hour(now_utc):
    """Suppress costly repeat audits, but force full input hashing at least hourly.

    This is only a metadata fast path; a same-size/same-mtime rewrite may escape
    the fast check for at most 60 minutes. Full audit records SHA256 content.
    """
    path=RESULTS/"github_actions_runs"/"evidence"/"research_manifest.json"
    try:
        obj=json.loads(path.read_text(encoding="utf-8"))
        if obj.get("schema")!="is-research-evidence-v1" or obj.get("integrity")!="PASS":
            return False
        generated=dt.datetime.fromisoformat(obj["generated_at"].replace("Z","+00:00"))
        age=(now_utc-generated).total_seconds()
        if age<0 or age>=3600:
            return False
        helper=pathlib.Path(__file__).resolve().with_name("harness_research_jobs.py")
        spec=importlib.util.spec_from_file_location("is_harness_signature_check",helper)
        if spec is None or spec.loader is None: return False
        module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        expected=obj.get("input_metadata_signature")
        return bool(expected and expected==module.input_metadata_signature())
    except (OSError,ValueError,KeyError,TypeError,AttributeError):
        return False


def is_pending_research_issue(issue):
    if "pull_request" in issue:
        return False
    labels={label.get("name") for label in issue.get("labels",[])}
    if not labels.intersection({"is-harness-queued","is-harness-running"}):
        return False
    try:
        return json.loads(issue.get("body") or "{}")=={"operation":"research-refresh"}
    except (ValueError,TypeError):
        return False

def ensure_scheduled_research_issue(token,now_utc):
    """Server-side fallback if GitHub cron is delayed; no catch-up or arbitrary jobs.

    Called from the existing locked user systemd unit, at most once per minute.
    Existing GitHub-scheduled Issues use exactly the same title and JSON body.
    """
    title=scheduled_slot_title(now_utc)
    # Five-minute heartbeat checks local source metadata first.
    # No new Issue is needed when all sources match the latest full audit.
    if unchanged_inputs_within_hour(now_utc):
        return "UNCHANGED_SOURCES_FAST_CHECK"
    recent=api("GET","/issues?state=all&sort=created&direction=desc&per_page=100",token)
    if any("pull_request" not in issue and issue.get("title")==title for issue in recent):
        return "SLOT_ALREADY_ENQUEUED"
    # Block backlogs if any previous refresh is still queued or running.
    active=api("GET","/issues?state=open&sort=updated&direction=desc&per_page=100",token)
    if any(is_pending_research_issue(issue) for issue in active):
        return "PENDING_REFRESH_EXISTS"
    created=api("POST","/issues",token,{
        "title":title,
        "body":json.dumps({"operation":"research-refresh"}),
        "labels":["is-harness-queued"],
    })
    return "FALLBACK_ENQUEUED_ISSUE_"+str(created["number"])

def record_heartbeat(now_utc, decision, path=None):
    """One local record per slot/status, separate from research results.

    A heartbeat establishes only that the watcher ran, not that science data
    changed or an MR/coloc was computed. Never publishes arbitrary server data.
    """
    path = pathlib.Path(path) if path is not None else RESULTS/"github_actions_runs"/"evidence"/"harness_heartbeat.json"
    slot = scheduled_slot_title(now_utc)
    event = {
        "schema":"is-harness-heartbeat-v1",
        "observed_at":now_utc.isoformat(),
        "slot_kst":slot,
        "watcher_status":decision,
        "is_new_scientific_result":False,
        "interpretation":"Input-change watcher only; scientific progress is separately verified.",
    }
    try:
        previous=json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
        if previous.get("slot_kst")==slot and previous.get("watcher_status")==decision:
            return "UNCHANGED_HEARTBEAT"
        path.parent.mkdir(parents=True,exist_ok=True)
        with tempfile.NamedTemporaryFile("w",dir=path.parent,prefix=".hb_",suffix=".json",
                                         encoding="utf-8",delete=False) as handle:
            json.dump(event,handle,ensure_ascii=False,indent=2)
            handle.write("\n")
            tmp=pathlib.Path(handle.name)
        try:
            os.replace(tmp,path)
        finally:
            tmp.unlink(missing_ok=True)
        return "HEARTBEAT_WRITTEN"
    except (OSError,ValueError,TypeError):
        return "HEARTBEAT_UNAVAILABLE"


def run_job(op):
    """Serialize approved research jobs across the queue worker and AUTO_NEXT.

    Both controllers must use this version of this module on the same host.
    This does not coordinate with other independently invoked research jobs.
    """
    import fcntl
    lock_path = pathlib.Path(os.environ.get(
        "IS_HARNESS_EXECUTION_LOCK",
        "/srv/is-analysis/results/is/github_actions_runs/is_research_execution.lock"))
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("SHARED_EXECUTION_LOCK_BUSY") from exc
        return _run_job_unlocked(op)


def _run_job_unlocked(op):
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
            try:
                poll_time=dt.datetime.now(dt.timezone.utc)
                fallback=ensure_scheduled_research_issue(token,poll_time)
                heartbeat_state=record_heartbeat(poll_time,fallback)
                if fallback not in {"UNCHANGED_SOURCES_FAST_CHECK"}:
                    print("five-minute scheduling:",fallback,"heartbeat:",heartbeat_state,flush=True)
            except Exception as fallback_error:
                # Scheduling error must never block manual or GitHub-enqueued jobs.
                print("five-minute scheduling unavailable:",
                      type(fallback_error).__name__,str(fallback_error)[:300],flush=True)
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
                    status="RETRY" if "SHARED_EXECUTION_LOCK_BUSY" in result else "FAIL"
                ended_at=dt.datetime.now(dt.timezone.utc)
                end_iso=ended_at.isoformat()
                duration_s=(ended_at-started_at).total_seconds()
                comment=f"IS HARNESS {status} | operation={op} | started_at={start_iso} | ended_at={end_iso} | duration_seconds={duration_s:.3f} | UTC={end_iso}\n\n```text\n{result[:12000]}\n```"
                api("POST",f"/issues/{num}/comments",token,{"body":comment})
                if status=="RETRY":
                    # Requeue only when the existing shared lock is currently busy.
                    # Other failures remain explicit failures for operator review.
                    api("POST",f"/issues/{num}/labels",token,{"labels":["is-harness-queued"]})
                    try: api("DELETE",f"/issues/{num}/labels/is-harness-running",token)
                    except urllib.error.HTTPError: pass
                else:
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
