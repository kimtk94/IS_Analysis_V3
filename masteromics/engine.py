"""DAG execution with content-addressed checkpoints and validated atomic outputs."""
from __future__ import annotations
import concurrent.futures as futures
import csv
import fcntl
import hashlib
import importlib.metadata
import math
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()


def atomic_json(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile('w', dir=path.parent, delete=False) as f:
        json.dump(value, f, indent=2, allow_nan=False); tmp = f.name
    os.replace(tmp, path)


def validate(path, contract):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0: raise ValueError(f'Missing/empty output: {path}')
    if contract['kind'] == 'json':
        value = json.loads(path.read_text())
        if not set(contract.get('keys', [])).issubset(value): raise ValueError(f'JSON keys: {path}')
        return {'sha256': sha(path), 'size_bytes': path.stat().st_size}
    if contract['kind'] != 'tsv': raise ValueError('Only explicit json/tsv contracts supported')
    with path.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not set(contract.get('columns', [])).issubset(reader.fieldnames or []):
            raise ValueError(f'Output schema mismatch: {path}')
        rows = 0
        for row in reader:
            rows += 1
            for col in contract.get('finite_columns', []):
                try: value = float(row[col])
                except (KeyError, TypeError, ValueError): raise ValueError(f'Invalid numeric output {col}: {path}')
                if not math.isfinite(value): raise ValueError(f'Nonfinite output {col}: {path}')
            for col in contract.get('probability_columns', []):
                if not 0 <= float(row[col]) <= 1: raise ValueError(f'Invalid probability {col}: {path}')
    if rows < contract.get('min_rows', 1): raise ValueError(f'Header-only/too few rows: {path}')
    return {'rows': rows, 'sha256': sha(path), 'size_bytes': path.stat().st_size}


def code_digest(root):
    h = hashlib.sha256()
    for p in sorted((root / 'masteromics').rglob('*')):
        if p.suffix in {'.py', '.R'}:
            h.update(str(p.relative_to(root)).encode()); h.update(p.read_bytes())
    return h.hexdigest()


def check_graph(stages):
    ids = [s['id'] for s in stages]
    if len(ids) != len(set(ids)): raise ValueError('Duplicate stage id')
    seen = set()
    for s in stages:
        if not s['id'].replace('_','').isalnum(): raise ValueError('Unsafe stage id')
        if not set(s.get('depends', [])).issubset(seen): raise ValueError('Stages must be topologically ordered')
        if not s.get('outputs'): raise ValueError('Every stage needs output contracts')
        seen.add(s['id'])


def execute(config_path, root, jobs=1, plan=False):
    cfg = json.loads(Path(config_path).read_text()); stages = cfg['stages']; check_graph(stages)
    variables = {**cfg.get('variables', {}), 'code': str(root), 'python': sys.executable}
    project = cfg['project']
    if not project.replace('_','').isalnum(): raise ValueError('Unsafe project')
    outroot = Path(cfg['output_root'].format_map(variables)).resolve()
    variables['results'] = str(outroot)
    def expand(v): return v.format_map(variables)
    rendered = []
    for s in stages:
        rendered.append({**s, 'argv': [expand(x) for x in s['argv']],
            'inputs': [expand(x) for x in s.get('inputs', [])],
            'outputs': [{**o, 'path': expand(o['path'])} for o in s['outputs']]})
    if plan:
        print(json.dumps({'project': project, 'stages': rendered}, indent=2)); return 0
    outroot.mkdir(parents=True, exist_ok=True)
    lock = (outroot / '.run.lock').open('w')
    try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError: raise RuntimeError(f'Project already running: {outroot}')
    cfgsha = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()
    code = code_digest(root); states = {}; failed = False
    runtime = {'python':sys.version, 'packages':{p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy']}}
    if any('Rscript' in s['argv'] or 'regional_gate' in s['argv'] for s in rendered):
        packages = ['TwoSampleMR','coloc','susieR','jsonlite']
        if any(s['id'] in ['cohort','incident'] for s in rendered): packages += ['lme4','survival']
        expression = 'cat(R.version.string); for(p in c(' + ','.join(json.dumps(p) for p in packages) + ')) cat(p,as.character(packageVersion(p)))'
        try: runtime['R'] = subprocess.check_output(['Rscript','-e',expression], text=True, stderr=subprocess.STDOUT)
        except (OSError, subprocess.CalledProcessError) as e:
            atomic_json(outroot/'run_summary.json',{'project':project,'status':'FAILED_ENVIRONMENT','error':str(e)})
            lock.close(); return 1
    atomic_json(outroot/'run_manifest.json', {'project':project, 'config_sha256':cfgsha,
        'code_sha256':code, 'runtime':runtime, 'started_at':time.time(), 'status':'RUNNING'})
    def run(s):
        started = time.time(); checkpoint = outroot/'state'/f"{s['id']}.json"
        record = {'stage':s['id'], 'started_at':started, 'argv':s['argv']}
        tmpdir = None
        try:
            inputs = {p:sha(p) for p in s['inputs']}
            # Executables/scripts outside core are part of identity too.
            external = {p:sha(p) for p in s['argv'] if Path(p).is_file()}
            fingerprint = hashlib.sha256(json.dumps([cfgsha,code,runtime,inputs,external,
                {d:states[d]['fingerprint'] for d in s.get('depends',[])}],sort_keys=True).encode()).hexdigest()
            record.update(fingerprint=fingerprint, inputs=inputs, external_code=external)
            if checkpoint.exists():
                old = json.loads(checkpoint.read_text())
                if old.get('fingerprint') == fingerprint and old.get('status') in ['SUCCESS_NONEMPTY','SKIPPED_EXISTING_VALID']:
                    try: current = {o['path']:validate(o['path'],o) for o in s['outputs']}
                    except (ValueError, OSError, json.JSONDecodeError): current = None
                    if current == old.get('outputs'):
                        return {**old, 'status':'SKIPPED_EXISTING_VALID'}
            # Child processes write into a private stage directory; {out0}, etc. are argv tokens.
            tmpdir = Path(tempfile.mkdtemp(prefix=s['id']+'-', dir=outroot))
            destinations = {f'@out{i}':str(tmpdir/Path(o['path']).name) for i,o in enumerate(s['outputs'])}
            argv = [destinations.get(a,a) for a in s['argv']]
            logs = outroot/'logs'; logs.mkdir(exist_ok=True)
            with (logs/f"{s['id']}.log").open('w') as log:
                subprocess.run(argv, cwd=root, stdout=log, stderr=subprocess.STDOUT, check=True,
                    timeout=s.get('timeout_seconds',86400))
            validated = {o['path']:validate(destinations[f'@out{i}'],o) for i,o in enumerate(s['outputs'])}
            for i,o in enumerate(s['outputs']):
                target = Path(o['path']); target.parent.mkdir(parents=True,exist_ok=True)
                os.replace(destinations[f'@out{i}'],target)
            record.update(status='SUCCESS_NONEMPTY', outputs=validated)
        except Exception as e:
            record.update(status='FAILED_RUNTIME', error=f'{type(e).__name__}: {e}')
        finally:
            if tmpdir: shutil.rmtree(tmpdir,ignore_errors=True)
        record['completed_at']=time.time(); atomic_json(checkpoint,record)
        return record
    pending = {s['id']:s for s in rendered}
    with futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        active = {}
        while pending or active:
            for key,s in list(pending.items()):
                deps = s.get('depends',[])
                if any(d in states and states[d]['status'].startswith(('FAILED','BLOCKED')) for d in deps):
                    states[key]={'status':'BLOCKED_DEPENDENCY'}; del pending[key]; failed=True
                elif all(d in states for d in deps) and len(active)<jobs:
                    active[pool.submit(run,s)] = key; del pending[key]
            if not active:
                if pending: raise RuntimeError('Unresolvable graph')
                break
            done,_ = futures.wait(active,return_when=futures.FIRST_COMPLETED)
            for f in done:
                key=active.pop(f); states[key]=f.result()
                print(f"[{project}/{key}] {states[key]['status']}",flush=True)
                failed |= states[key]['status'].startswith('FAILED')
    atomic_json(outroot/'run_summary.json',{'project':project,'status':'FAILED' if failed else 'SUCCESS', 'stages':states})
    atomic_json(outroot/'run_manifest.json',{'project':project,'status':'FAILED' if failed else 'SUCCESS',
        'config_sha256':cfgsha,'code_sha256':code,'completed_at':time.time(),'runtime':runtime})
    lock.close(); return 1 if failed else 0
