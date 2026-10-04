import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "workflow" / "run_master.py"

def run(*args):
    return subprocess.run(
        [sys.executable, str(RUNNER), *args],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )

def test_ckd_stage_0_6_dry_run():
    p = run("--disease","ckd","--stages","0-6","--dry-run")
    assert p.returncode == 0, p.stderr
    for stage in range(7):
        assert f"[{stage:02d}]" in p.stdout
    assert "primary_mr" in p.stdout
    assert "coloc" in p.stdout
    assert "finemap" in p.stdout

def test_is_stage_0_6_exposes_gaps():
    p = run("--disease","is","--stages","0-6","--dry-run")
    assert p.returncode == 0, p.stderr
    assert "BLOCKED" in p.stdout
    assert "exposure" in p.stdout
    assert "READY" in p.stdout

def test_execute_refuses_blocked_plan():
    p = run("--disease","is","--stages","0-3","--execute")
    assert p.returncode == 2
    assert "non-executable" in p.stderr

def test_invalid_stage_fails():
    p = run("--disease","ckd","--stages","99","--dry-run")
    assert p.returncode != 0


def test_stage4_becomes_ready_with_explicit_io(tmp_path):
    inp=tmp_path/"harm.tsv"
    out=tmp_path/"robust.tsv"
    inp.write_text("protein_id\tgene_symbol\tancestry\tphenotype\tbeta_exposure\tbeta_outcome\tse_outcome\n",encoding="utf-8")
    p=run("--disease","ckd","--stages","4","--dry-run",
          "--harmonized-input",str(inp),"--robust-output",str(out))
    assert p.returncode==0,p.stderr
    assert "[04]" in p.stdout
    assert "READY" in p.stdout
    assert "run_master_mr_robustness.py" in p.stdout


def test_stage6_becomes_ready_with_explicit_io(tmp_path):
    inp=tmp_path/"susie_input"; ld=tmp_path/"ld"; out=tmp_path/"susie_out"
    inp.mkdir(); ld.mkdir()
    abf=tmp_path/"abf.tsv"
    abf.write_text("locus\tPP.H3\tPP.H4\nGENE1\t0.1\t0.9\n",encoding="utf-8")
    p=run(
      "--disease","is","--stages","6","--dry-run",
      "--susie-input-dir",str(inp),
      "--ld-dir",str(ld),
      "--abf-file",str(abf),
      "--susie-output-dir",str(out),
      "--ancestry","EAS",
      "--outcome-type","cc",
    )
    assert p.returncode==0,p.stderr
    assert "[06]" in p.stdout
    assert "READY" in p.stdout
    assert "run_master_susie.R" in p.stdout
    assert " EAS cc" in p.stdout
