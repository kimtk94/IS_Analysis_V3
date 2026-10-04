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


def test_stage7_becomes_ready_with_explicit_io(tmp_path):
    disc=tmp_path/"disc.tsv"; repl=tmp_path/"repl.tsv"; out=tmp_path/"cross.tsv"
    header="protein_id\tgene_symbol\tphenotype\tbeta\tse\tp\n"
    disc.write_text(header+"P1\tG1\tIS\t0.2\t0.05\t0.001\n",encoding="utf-8")
    repl.write_text(header+"P1\tG1\tIS\t0.1\t0.06\t0.05\n",encoding="utf-8")
    p=run(
      "--disease","is","--stages","7","--dry-run",
      "--discovery-mr",str(disc),
      "--replication-mr",str(repl),
      "--cross-ancestry-output",str(out),
      "--discovery-ancestry","EUR",
      "--replication-ancestry","EAS",
    )
    assert p.returncode==0,p.stderr
    assert "[07]" in p.stdout
    assert "READY" in p.stdout
    assert "run_master_cross_ancestry.py" in p.stdout


def test_stage8_becomes_ready_with_explicit_io(tmp_path):
    d=tmp_path/"olink.tsv"; r=tmp_path/"soma.tsv"; out=tmp_path/"platform.tsv"
    header="gene_symbol\tphenotype\tbeta\tse\tp\n"
    d.write_text(header+"F11\tIS\t0.2\t0.05\t0.001\n",encoding="utf-8")
    r.write_text(header+"F11\tIS\t0.1\t0.06\t0.05\n",encoding="utf-8")
    p=run(
      "--disease","is","--stages","8","--dry-run",
      "--platform-discovery-mr",str(d),
      "--platform-replication-mr",str(r),
      "--platform-output",str(out),
      "--discovery-platform","Olink",
      "--replication-platform","SomaScan",
    )
    assert p.returncode==0,p.stderr
    assert "[08]" in p.stdout
    assert "READY" in p.stdout
    assert "run_master_cross_platform.py" in p.stdout
