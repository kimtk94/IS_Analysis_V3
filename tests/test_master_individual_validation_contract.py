from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"run_master_individual_validation.R"

def test_individual_validation_r_contract():
    text=SCRIPT.read_text(encoding="utf-8")
    assert "baseline_egfr" in text
    assert "egfr_slope" in text
    assert "incident_ckd_logistic" in text
    assert "survival::coxph" in text
    assert "lme4::lmer" in text
    assert "time-by-predictor interaction" in text
    assert "MASTER_INDIVIDUAL_VALIDATION_PASS" in text
