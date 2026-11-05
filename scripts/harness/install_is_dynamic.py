#!/usr/bin/env python3
"""Install the IS-only dynamic runner into existing standalone Harness.

Only changes the IS route of experiment_runner.py; Kalman routes are untouched.
No execution of research experiments during installation.
"""
from pathlib import Path
import argparse
import os
import py_compile
import shutil
import sys
import tempfile

DEFAULT_HARNESS = Path("/home/taehoon/agent-harness")
ANCHOR = " if not 1<=x.steps<=20:a.error('--steps must be 1..20')"
MARKER = "# IS_DYNAMIC_ROUTER_20261013"
ROUTE = ''' # IS_DYNAMIC_ROUTER_20261013
 if x.job == 'is-research':
  import subprocess as _is_sp
  _cmd = [sys.executable, str(BASE/'is_dynamic_loop.py'), '--steps', str(x.steps)]
  _cmd.append('--run' if x.run else '--status')
  _response = _is_sp.run(_cmd, cwd=str(BASE), capture_output=True, text=True)
  if _response.returncode != 0:
   print(json.dumps({'jobs':{'is-research':{'state':'IS_DYNAMIC_RUNNER_ERROR',
       'exit':_response.returncode, 'stderr':_response.stderr[-1500:]}}}))
   raise SystemExit(_response.returncode)
  try:
   _out = json.loads(_response.stdout)
  except ValueError:
   raise SystemExit('IS dynamic runner returned invalid JSON')
  print(json.dumps({'jobs':{'is-research':_out}},ensure_ascii=False))
  return
'''
def install(harness):
    destination = harness / "experiment_runner.py"
    source = Path(__file__).with_name("is_dynamic_loop.py")
    if not source.is_file() or not destination.is_file():
        raise SystemExit("missing script or Harness runner")
    current = destination.read_text()
    if MARKER in current:
        print("ALREADY_INSTALLED")
        return
    if current.count(ANCHOR) != 1:
        raise SystemExit("runner version mismatch: refusing patch")
    patched = current.replace(ANCHOR, ROUTE + ANCHOR, 1)
    backup = destination.with_suffix(".py.pre_is_dynamic_20261013.bak")
    if backup.exists():
        raise SystemExit("backup already exists; refusing unexpected overwrite")
    # First validate both sources without writing to the live runner.
    compile(patched, str(destination), "exec")
    compile(source.read_text(), str(source), "exec")
    shutil.copy2(destination, backup)
    shutil.copy2(source, harness / "is_dynamic_loop.py")
    temp = destination.with_suffix(".is_dynamic_install.tmp")
    try:
        temp.write_text(patched)
        os.replace(temp, destination)
    except Exception:
        if temp.exists():
            temp.unlink()
        shutil.copy2(backup, destination)
        raise
    print("INSTALLED_IS_ONLY", "backup=", backup, "module=", harness / "is_dynamic_loop.py")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--harness", type=Path, default=DEFAULT_HARNESS)
    p.add_argument("--apply", action="store_true")
    a = p.parse_args()
    if not a.apply:
        print("DRY_RUN: pass --apply to patch IS route only; Kalman is untouched")
        return
    install(a.harness)

if __name__ == "__main__":
    main()
