import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = [ROOT/"config/master_ckd.json", ROOT/"config/master_is.json"]

def test_master_configs_are_valid_json():
    for path in CONFIGS:
        with path.open() as f:
            cfg = json.load(f)
        assert cfg["schema_version"] == 1
        assert cfg["disease"] in {"ckd","is"}

def test_master_validator_passes():
    cmd = [sys.executable, str(ROOT/"scripts/validate_master_pipeline.py"), *map(str, CONFIGS)]
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + "\n" + p.stderr
    assert p.stdout.count("PASS") == 2

def test_stage_count_is_complete():
    expected = 20
    for path in CONFIGS:
        cfg = json.loads(path.read_text())
        assert len(cfg["stages"]) == expected
