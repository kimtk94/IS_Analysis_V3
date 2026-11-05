#!/usr/bin/env python3
from __future__ import annotations
import json,os
from pathlib import Path
REPO=Path(os.environ.get('IS_ANALYSIS_REPO',str(Path(__file__).resolve().parents[3])))
MAN=REPO/'docs/metabolic_resilience/DRIVE_SOURCE_MANIFEST.json'
CODE=REPO/'scripts/metabolic_resilience/drive_v2'
if not MAN.exists(): raise SystemExit('Drive manifest missing')
p=json.loads(MAN.read_text())
drive={x['title'] for x in p.get('files',[]) if x['title'].endswith(('.py','.sh','.R'))}
git={x.name for x in CODE.iterdir() if x.is_file() and x.suffix in {'.py','.sh','.R'}}
repo_native_allow={'audit_metabolic_resilience_koges_inputs.py'}
git_mirror=git-repo_native_allow
res={'drive_only':sorted(drive-git_mirror),'git_only':sorted(git_mirror-drive),'repo_native_allow':sorted(repo_native_allow & git),'matched':len(drive&git_mirror),'status':'PASS' if drive==git_mirror else 'DRIFT'}
out=REPO/'docs/metabolic_resilience/DRIVE_GIT_DRIFT_AUDIT.json';out.write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2));raise SystemExit(0 if drive==git_mirror else 2)
