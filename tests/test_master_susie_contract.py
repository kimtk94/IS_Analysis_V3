from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"run_master_susie.R"

def test_master_susie_contract():
    text=SCRIPT.read_text(encoding="utf-8")
    assert "coloc.susie" in text
    assert "runsusie" in text
    assert "case_fraction" in text
    assert "outcome_type" in text
    assert "LD ancestry mismatch" in text
    assert "effect_vs_ld_major_sign" in text
    assert "MASTER_SUSIE_DEFAULT.tsv" in text
    assert "MASTER_SUSIE_FAILURES.tsv" in text
