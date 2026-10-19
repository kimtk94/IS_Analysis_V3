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


def test_stage9_becomes_ready_with_explicit_io(tmp_path):
    pmr=tmp_path/"pmr.tsv"; smr=tmp_path/"smr.tsv"; out=tmp_path/"transcript.tsv"
    pmr.write_text("gene_symbol\tphenotype\tbeta\tse\tp\nG1\tCKD\t0.2\t0.05\t0.001\n",encoding="utf-8")
    smr.write_text("Gene\tProbeID\tb_SMR\tse_SMR\tp_SMR\tp_HEIDI\nG1\tP1\t0.3\t0.08\t0.002\t0.5\n",encoding="utf-8")
    p=run(
      "--disease","ckd","--stages","9","--dry-run",
      "--transcript-protein-mr",str(pmr),
      "--smr-results",str(smr),
      "--transcript-phenotype","CKD",
      "--transcript-output",str(out),
    )
    assert p.returncode==0,p.stderr
    assert "[09]" in p.stdout
    assert "READY" in p.stdout
    assert "run_master_transcriptomics.py" in p.stdout


def test_stage10_and_12_become_ready_with_manifest(tmp_path):
    m=tmp_path/"manifest.tsv"; long=tmp_path/"long.tsv"; wide=tmp_path/"wide.tsv"
    m.write_text("mr_file\tlabel\tgroup\tcoloc_file\n",encoding="utf-8")
    p=run(
      "--disease","is","--stages","10,12","--dry-run",
      "--phenotype-manifest",str(m),
      "--phenotype-long-output",str(long),
      "--phenotype-wide-output",str(wide),
    )
    assert p.returncode==0,p.stderr
    assert "[10]" in p.stdout and "[12]" in p.stdout
    assert p.stdout.count("READY")>=2
    assert "run_master_phenotype_matrix.py" in p.stdout

def test_stage11_becomes_ready_with_explicit_io(tmp_path):
    d=tmp_path/"d.tsv"; r=tmp_path/"r.tsv"; out=tmp_path/"risk.tsv"
    header="gene_symbol\tphenotype\tbeta\tse\tp\n"
    d.write_text(header+"FURIN\tIS\t0.2\t0.05\t0.001\n",encoding="utf-8")
    r.write_text(header+"FURIN\tSBP\t0.3\t0.08\t0.002\n",encoding="utf-8")
    p=run(
      "--disease","is","--stages","11","--dry-run",
      "--risk-disease-mr",str(d),
      "--risk-mr",str(r),
      "--risk-disease-phenotype","IS",
      "--risk-output",str(out),
    )
    assert p.returncode==0,p.stderr
    assert "[11]" in p.stdout
    assert "READY" in p.stdout
    assert "run_master_risk_factor.py" in p.stdout


def test_stage13_becomes_ready_with_explicit_io(tmp_path):
    subj=tmp_path/"subjects.tsv"; out=tmp_path/"individual.tsv"
    subj.write_text(
      "participant_id\tdosage\tbaseline_f0_egfr\tegfr_slope\tincident_ckd\tage\tsex\n"
      "1\t0\t90\t-1.0\t0\t50\t1\n",
      encoding="utf-8"
    )
    p=run(
      "--disease","ckd","--stages","13","--dry-run",
      "--individual-subject-file",str(subj),
      "--individual-predictor","dosage",
      "--individual-output",str(out),
      "--baseline-endpoint","baseline_f0_egfr",
      "--slope-endpoint","egfr_slope",
      "--incident-endpoint","incident_ckd",
      "--individual-covariates","age,sex",
    )
    assert p.returncode==0,p.stderr
    assert "[13]" in p.stdout
    assert "READY" in p.stdout
    assert "run_master_individual_validation.R" in p.stdout


def test_stage14_16_become_ready_with_localization_io(tmp_path):
    cand=tmp_path/"cand.tsv"; out=tmp_path/"loc.tsv"; tm=tmp_path/"tm.tsv"
    cand.write_text("gene_symbol\nUMOD\n",encoding="utf-8")
    tm.write_text("source_file\tsource_name\ttissue_or_region\n",encoding="utf-8")
    p=run(
      "--disease","ckd","--stages","14,15,16","--dry-run",
      "--localization-candidates",str(cand),
      "--tissue-manifest",str(tm),
      "--localization-output",str(out),
    )
    assert p.returncode==0,p.stderr
    assert "[14]" in p.stdout and "[15]" in p.stdout and "[16]" in p.stdout
    assert p.stdout.count("READY")>=3
    assert "run_master_localization.py" in p.stdout


def test_stage17_becomes_ready_with_manifest(tmp_path):
    m=tmp_path/"phewas_manifest.tsv"; long=tmp_path/"phewas_long.tsv"; summary=tmp_path/"phewas_summary.tsv"
    m.write_text("source_file\tsource_name\tphenotype_group\n",encoding="utf-8")
    p=run(
      "--disease","is","--stages","17","--dry-run",
      "--phewas-manifest",str(m),
      "--phewas-long-output",str(long),
      "--phewas-summary-output",str(summary),
    )
    assert p.returncode==0,p.stderr
    assert "[17]" in p.stdout
    assert "READY" in p.stdout
    assert "run_master_phewas.py" in p.stdout

def test_stage18_becomes_ready_with_drug_table(tmp_path):
    d=tmp_path/"drug.tsv"; out=tmp_path/"drug_out.tsv"
    d.write_text("gene_symbol\tdrug_name\nF11\tDrugA\n",encoding="utf-8")
    p=run(
      "--disease","is","--stages","18","--dry-run",
      "--drug-table",str(d),
      "--drug-output",str(out),
    )
    assert p.returncode==0,p.stderr
    assert "[18]" in p.stdout
    assert "READY" in p.stdout
    assert "run_master_druggability.py" in p.stdout

def test_stage19_becomes_ready_with_evidence_manifest(tmp_path):
    cand=tmp_path/"cand.tsv"; m=tmp_path/"ev.tsv"; out=tmp_path/"final.tsv"
    cand.write_text("gene_symbol\nF11\n",encoding="utf-8")
    m.write_text("stage\tfile\tkey_column\n",encoding="utf-8")
    p=run(
      "--disease","is","--stages","19","--dry-run",
      "--evidence-candidates",str(cand),
      "--evidence-manifest",str(m),
      "--evidence-output",str(out),
    )
    assert p.returncode==0,p.stderr
    assert "[19]" in p.stdout
    assert "READY" in p.stdout
    assert "run_master_evidence_integration.py" in p.stdout


def test_stage0_becomes_ready_with_registry_io(tmp_path):
    manifest=tmp_path/"datasets.tsv"; reg=tmp_path/"registry.tsv"; summary=tmp_path/"summary.json"
    manifest.write_text("dataset_id\trole\tpath\tancestry\tgenome_build\tphenotype\tplatform\tsource_name\n",encoding="utf-8")
    p=run(
      "--disease","ckd","--stages","0","--dry-run",
      "--dataset-manifest",str(manifest),
      "--dataset-registry-output",str(reg),
      "--dataset-registry-summary",str(summary),
    )
    assert p.returncode==0,p.stderr
    assert "[00]" in p.stdout and "READY" in p.stdout
    assert "run_master_dataset_registry.py" in p.stdout

def test_stage1_becomes_ready_with_exposure_io(tmp_path):
    src=tmp_path/"src.tsv"; cmap=tmp_path/"map.tsv"; out=tmp_path/"norm.tsv"
    src.write_text("x\n",encoding="utf-8")
    cmap.write_text("canonical_column\tsource_column\n",encoding="utf-8")
    p=run(
      "--disease","is","--stages","1","--dry-run",
      "--exposure-input",str(src),
      "--exposure-column-map",str(cmap),
      "--exposure-output",str(out),
      "--exposure-ancestry","EUR",
      "--exposure-genome-build","GRCh37",
      "--exposure-platform","Olink",
    )
    assert p.returncode==0,p.stderr
    assert "[01]" in p.stdout and "READY" in p.stdout
    assert "run_master_exposure_normalize.py" in p.stdout

def test_stage2_becomes_ready_with_instrument_io(tmp_path):
    src=tmp_path/"norm.tsv"; out=tmp_path/"qc.tsv"
    src.write_text("protein_id\tgene_symbol\trsid\teffect_allele\tother_allele\tbeta\tse\n",encoding="utf-8")
    p=run(
      "--disease","ckd","--stages","2","--dry-run",
      "--instrument-input",str(src),
      "--instrument-output",str(out),
    )
    assert p.returncode==0,p.stderr
    assert "[02]" in p.stdout and "READY" in p.stdout
    assert "run_master_instrument_qc.py" in p.stdout


def test_stage3_becomes_ready_with_generic_mr_io(tmp_path):
    exp=tmp_path/"exp.tsv"; outc=tmp_path/"out.tsv"; harm=tmp_path/"harm.tsv"; mr=tmp_path/"mr.tsv"
    exp.write_text("protein_id\tgene_symbol\trsid\teffect_allele\tother_allele\tbeta\tse\n",encoding="utf-8")
    outc.write_text("rsid\teffect_allele\tother_allele\tbeta\tse\n",encoding="utf-8")
    p=run(
      "--disease","is","--stages","3","--dry-run",
      "--mr-exposure",str(exp),
      "--mr-outcome",str(outc),
      "--mr-harmonized-output",str(harm),
      "--mr-output",str(mr),
      "--mr-ancestry","EUR",
      "--mr-phenotype","IS",
    )
    assert p.returncode==0,p.stderr
    assert "[03]" in p.stdout and "READY" in p.stdout
    assert "run_master_mr.py" in p.stdout

def test_stage5_becomes_ready_with_generic_coloc_io(tmp_path):
    inp=tmp_path/"coloc_in"; out=tmp_path/"coloc_out"
    inp.mkdir()
    p=run(
      "--disease","is","--stages","5","--dry-run",
      "--coloc-input-dir",str(inp),
      "--coloc-output-dir",str(out),
      "--coloc-outcome-type","cc",
    )
    assert p.returncode==0,p.stderr
    assert "[05]" in p.stdout and "READY" in p.stdout
    assert "run_master_coloc.R" in p.stdout
    assert " cc " in p.stdout
