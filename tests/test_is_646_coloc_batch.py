"""Fail-closed contracts for original 646-coloc batch and EAS LD coverage."""
import ast
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
R=ROOT/"scripts/is/replay_is_legacy_646_batch.R"
BUILD=ROOT/"scripts/is/build_is_646_snp_replay_colab.py"
LD=ROOT/"scripts/is/audit_is_1000g_eas_ld_priority8.py"
NOTEBOOK=ROOT/"notebooks/is/IS_PHASE10_11_646_SNP_COLOC_REPLAY_ONECELL.ipynb"

def load(path):
    s=importlib.util.spec_from_file_location(path.stem,path)
    m=importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m

class TestIS646Batch(unittest.TestCase):
    def test_onecell_embeds_exact_frozen_R_source(self):
        n=json.loads(NOTEBOOK.read_text(encoding="utf8"))
        self.assertEqual(len(n["cells"]),1)
        code="".join(n["cells"][0]["source"])
        ast.parse(code)
        self.assertIn("IS_PHASE10_11_646_BATCH=RUN_ENDED_AND_STATUS_ACCOUNTED",code)
        self.assertIn("MISSING_INPUT",code)
        self.assertIn("GTEX",code.upper())
        self.assertIn("coloc_5.2.3.tar.gz",code)
        self.assertIn("IS_646_SOURCE_FILE_SHA256.tsv",code)
        self.assertIn("viridis",code)
        embedded=R.read_text(encoding="utf8")
        self.assertIn(repr(embedded),code)
        self.assertIn('n_available=sum',code)
        self.assertIn("INPUTS_MISSING_ON_DRIVE=",code)

    def test_r_full_batch_parse(self):
        r=subprocess.run(["Rscript","--vanilla","-e",
            "invisible(parse(file=commandArgs(trailingOnly=TRUE)[[1]]));cat('R_646_PARSE_OK\\n')",
            str(R)],capture_output=True,text=True,timeout=20)
        self.assertEqual(r.returncode,0,msg=r.stderr)
        self.assertIn("R_646_PARSE_OK",r.stdout)

    def test_builder_no_overwrite(self):
        m=load(BUILD)
        with tempfile.TemporaryDirectory() as t:
            out=Path(t)/"one.ipynb"
            m.build(R,out)
            n=json.loads(out.read_text(encoding="utf8"))
            self.assertEqual(len(n["cells"]),1)
            before=out.read_bytes()
            with self.assertRaises(FileExistsError):
                m.build(R,out)
            self.assertEqual(before,out.read_bytes())

    def test_ld_audit_fails_closed_on_nonempty_output(self):
        m=load(LD)
        with tempfile.TemporaryDirectory() as t:
            out=Path(t)/"out"
            out.mkdir()
            sentinel=out/"historical.txt"
            sentinel.write_text("DO_NOT_OVERWRITE")
            with self.assertRaises(FileExistsError):
                m.run(Path(t)/"missing_input",Path(t)/"missing_report",
                      Path(t)/"missing_ld",out)
            self.assertEqual(sentinel.read_text(),"DO_NOT_OVERWRITE")

    def test_ld_audit_missing_input_cannot_generate_pass(self):
        m=load(LD)
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            report=root/"empty.tsv"
            report.write_text("gene\tlocus\tdataset_key\tnsnps\n")
            out=root/"out"
            with self.assertRaises(ValueError):
                m.run(root,report,root,out)
            self.assertFalse(out.exists())

if __name__=="__main__":
    unittest.main()
