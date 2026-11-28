"""Exercise R cache behavior with a synthetic fitter and base R only."""
from pathlib import Path
import shutil
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[1]


class SusieCacheTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which("Rscript"), "Rscript unavailable; base-R cache regression runs in CI")
    def test_cache_invalidation_and_legacy_migration(self):
        result = subprocess.run(
            [shutil.which("Rscript"), str(ROOT / "tests" / "test_ckd_stage2c_susie_cache.R"), str(ROOT)],
            cwd=ROOT, text=True, capture_output=True, timeout=60,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("CKD_SUSIE_CACHE_TEST_PASS", result.stdout)


if __name__ == "__main__":
    unittest.main()
