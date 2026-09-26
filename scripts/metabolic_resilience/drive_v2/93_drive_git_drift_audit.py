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
res={'drive_only':sorted(drive-git),'git_only':sorted(git-drive),'matched':len(drive&git),'status':'PASS' if drive==git else 'DRIFT'}
out=REPO/'docs/metabolic_resilience/DRIVE_GIT_DRIFT_AUDIT.json';out.write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2));raise SystemExit(0 if drive==git else 2)
