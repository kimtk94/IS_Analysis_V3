import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import publication_figures as pf

ROOT = Path(__file__).resolve().parent / "fixtures" / "publication_figures"
MANIFEST = ROOT / "demo.json"


class PublicationFiguresTests(unittest.TestCase):
    def test_fixture(self):
        spec, records = pf.load(MANIFEST, ROOT)
        self.assertEqual(spec["data_status"], "synthetic")
        self.assertEqual({r["type"] for r in records}, pf.TYPES)
        self.assertTrue(all(len(r["sha256"]) == 64 for r in records if "source" in r))

    def test_invalid_probabilities(self):
        for v in (0, -.1, 1.1, "nan", "inf"):
            with self.subTest(v=v), self.assertRaises(ValueError):
                pf.prob(v)

    def test_path_escape(self):
        for v in ("../gigastroke.tsv", str(MANIFEST)):
            with self.assertRaises(ValueError):
                pf.source_path(ROOT, v)

    def test_documented_edge_reference(self):
        spec = json.loads(MANIFEST.read_text())
        spec["panels"][0]["edges"][0]["kind"] = "documented"
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "spec.json"
            path.write_text(json.dumps(spec))
            with self.assertRaisesRegex(ValueError, "reference"):
                pf.load(path, ROOT)

    def test_no_individual_input(self):
        spec = json.loads(MANIFEST.read_text())
        spec["data_status"] = "individual"
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "spec.json"
            path.write_text(json.dumps(spec))
            with self.assertRaises(ValueError):
                pf.load(path, ROOT)

    @unittest.skipUnless(importlib.util.find_spec("matplotlib") and importlib.util.find_spec("PIL"),
                         "Optional publication dependencies not installed")
    def test_render_and_refuse_overwrite(self):
        from PIL import Image
        import xml.etree.ElementTree as ET
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "result"
            pf.render(MANIFEST, ROOT, out)
            report = json.loads((out / "provenance.json").read_text())
            self.assertEqual(report["status"], "DRAFT_REVIEW_REQUIRED")
            self.assertEqual(set(report["exports"]), {"figure.svg", "figure.pdf", "figure.png", "figure.tiff"})
            for name, meta in report["exports"].items():
                self.assertEqual(pf.digest(out / name), meta["sha256"])
            svg = (out / "figure.svg").read_text()
            ET.fromstring(svg)
            self.assertIn("SYNTHETIC DEMO", svg)
            self.assertEqual((out / "figure.pdf").read_bytes()[:5], b"%PDF-")
            for ext in ("png", "tiff"):
                with Image.open(out / ("figure." + ext)) as im:
                    self.assertAlmostEqual(im.info["dpi"][0], 600, delta=1)
            with self.assertRaisesRegex(ValueError, "refusing overwrite"):
                pf.render(MANIFEST, ROOT, out)


if __name__ == "__main__":
    unittest.main()
