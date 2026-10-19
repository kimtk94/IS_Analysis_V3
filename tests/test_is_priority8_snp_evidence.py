"""Priority-eight original SNP-level ABF audit contracts."""
import importlib.util
import subprocess
import tempfile
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
PY=ROOT/"scripts/is/summarize_is_priority8_snp_evidence.py"
RR=ROOT/"scripts/is/replay_is_phase10_abf_priority8_snp_inputs.R"
MAF=ROOT/"scripts/is/audit_is_legacy_crossancestry_maf_priority8.R"
COMPARATOR=ROOT/"scripts/is/audit_is_priority8_snp_vs_summary_prior_grid.py"


def mod(path):
    spec=importlib.util.spec_from_file_location(path.stem,path)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


class PriorityEightGates(unittest.TestCase):
    def test_prespecified_eight_unique_gene_ids(self):
        m=mod(PY)
        self.assertEqual(len(m.GENES),8)
        self.assertEqual(len(set(m.GENES)),8)
        self.assertEqual(len(m.PRIORS),5)

    def test_never_overwrite_report_folder(self):
        m=mod(PY)
        with tempfile.TemporaryDirectory() as t:
            folder=Path(t)
            output=folder/"prior"
            output.mkdir()
            (output/"old.txt").write_text("DO_NOT_MUTATE")
            with self.assertRaises(FileExistsError):
                m.run(folder/"missing1",folder/"missing2",folder/"missing3",folder/"missing4",output)
            self.assertEqual((output/"old.txt").read_text(),"DO_NOT_MUTATE")

    def test_missing_sources_cause_failure(self):
        m=mod(PY)
        with tempfile.TemporaryDirectory() as t:
            folder=Path(t)
            out=folder/"new"
            with self.assertRaises(FileNotFoundError):
                m.run(folder/"missing1",folder/"missing2",folder/"missing3",folder/"missing4",out)
            self.assertFalse(out.exists())

    def test_r_scripts_parse_clean(self):
        for r in [RR,MAF]:
            proc=subprocess.run(["Rscript","--vanilla","-e",
                "parse(file=commandArgs(trailingOnly=TRUE)[1]);cat('PARSE_OK\\n')",
                str(r)],capture_output=True,text=True,timeout=30)
            if proc.returncode!=0 and "No such file" in proc.stderr:
                self.skipTest("Rscript not installed")
            self.assertEqual(proc.returncode,0,proc.stderr)
            self.assertIn("PARSE_OK",proc.stdout)

    def test_summary_prior_match_script_imports(self):
        m=mod(COMPARATOR)
        self.assertEqual(len(m.GENE_ID),8)
        self.assertEqual(set(m.GENE_ID),set(mod(PY).GENES))


if __name__=="__main__":
    unittest.main()
