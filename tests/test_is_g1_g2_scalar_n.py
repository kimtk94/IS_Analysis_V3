"""G1 structural allele and G2 scalar-QTL-N stress contracts, fail closed."""
import ast
import importlib.util
import json
from pathlib import Path
import tempfile
import subprocess
import unittest

ROOT=Path(__file__).resolve().parents[1]
G1=ROOT/"scripts/is/audit_is_646_g1_source_alleles.py"
R=ROOT/"scripts/is/is_646_qtl_scalar_n_sensitivity.R"
BUILDER=ROOT/"scripts/is/build_is_646_qtl_n_colab.py"
TEMPLATE=ROOT/"notebooks/is/IS_PHASE10_11_646_SNP_COLOC_REPLAY_ONECELL.ipynb"
NOTEBOOK=ROOT/"notebooks/is/IS_PHASE10_11_646_QTL_N_SENSITIVITY_ONECELL.ipynb"

def load(path):
    spec=importlib.util.spec_from_file_location(path.stem,path)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def record():
    return dict(
        variant="chr10_103079395_T_C",
        chromosome="10",position="103079395",ref="T",alt="C",
        match_key="10:103079395:T:C",
        variant_id_hg38="10:103079395:T:C",
        gwas_eaf="0.37717391182341",gwas_maf="0.37717391182341",
        eqtl_maf="0.471429",maf="0.471429",
        eqtl_beta="0.0599577",beta="0.0599577",
        eqtl_se="0.0964537",se="0.0964537",
        harmonization="EXACT_REF_ALT_GRCH38",
        variant_id="10:104839152:T:C")

class TestISG1Structural(unittest.TestCase):
    def test_readable_valid_record(self):
        r=load(G1).parse_record(record())
        self.assertEqual(r["is_palindromic"],0)
        self.assertEqual(r["is_gwas_gtEx_maf_diff_gt_0_1"],0)

    def test_would_catch_grch38_position_or_ref_alt_mismatch(self):
        r=record();r["variant_id_hg38"]="10:103079395:C:T"
        with self.assertRaisesRegex(ValueError,"GTEX_GRCH38"):
            load(G1).parse_record(r)

    def test_would_catch_maf_effect_source_disagreement(self):
        r=record();r["eqtl_beta"]="0.1"
        with self.assertRaisesRegex(ValueError,"EQTL_BETA"):
            load(G1).parse_record(r)

    def test_gwas_eaf_maf_internal_mismatch_fails(self):
        r=record();r["gwas_maf"]="0.1"
        with self.assertRaisesRegex(ValueError,"GWAS_EAF_MAF"):
            load(G1).parse_record(r)

    def test_existing_directory_not_overwritten(self):
        with tempfile.TemporaryDirectory() as t:
            out=Path(t)/"out";out.mkdir()
            f=out/"sentinel.txt";f.write_text("ORIGINAL")
            with self.assertRaises(FileExistsError):
                load(G1).run(Path(t)/"absent",Path(t)/"absent",out)
            self.assertEqual(f.read_text(),"ORIGINAL")

class TestISNScalarColab(unittest.TestCase):
    def test_embedded_r_exact_and_singlecell_ast(self):
        note=json.loads(NOTEBOOK.read_text(encoding="utf8"))
        self.assertEqual(len(note["cells"]),1)
        py="".join(note["cells"][0]["source"])
        ast.parse(py)
        self.assertIn(repr(R.read_text(encoding="utf8")),py)
        self.assertIn("IS_PHASE10_11_646_QTL_N_SCALAR_SENSITIVITY_V1",py)
        self.assertIn("THREE_SCALAR_N_STRESS_NOT_SNP_SPECIFIC_MODEL",py)
        self.assertIn("SOURCE_SHA_CHANGED",py)
        self.assertIn("IS_646_QTL_N_SCALAR_SENSITIVITY.tsv",py)
        self.assertIn("IS_646_N_SCALAR_COLAB_MANIFEST.json",py)
        self.assertIn("G1_ORIGINAL_ALLELE_STRUCTURAL_V1",py)
        self.assertIn(repr(G1.read_text(encoding="utf8")),py)
        self.assertIn("g1_source_audit_status",py)
        self.assertIn("G1_STRUCTURAL_INCOMPLETE_OR_FAIL",py)
        self.assertNotIn("results/IS_PHASE10_11_LEGACY_646_SNP_REPLAY_V1",py)

    def test_r_sensitivity_source_parses(self):
        result=subprocess.run(["Rscript","--vanilla","-e",
            "invisible(parse(file=commandArgs(trailingOnly=TRUE)[1]));cat('SCALAR_R_PARSE_PASS\\n')",
            str(R)],capture_output=True,text=True,timeout=15)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn("SCALAR_R_PARSE_PASS",result.stdout)

    def test_notebook_builder_no_overwrite(self):
        with tempfile.TemporaryDirectory() as t:
            out=Path(t)/"note.ipynb"
            build=load(BUILDER).build
            build(TEMPLATE,R,out)
            data=out.read_bytes()
            with self.assertRaises(FileExistsError):
                build(TEMPLATE,R,out)
            self.assertEqual(data,out.read_bytes())

if __name__=="__main__":
    unittest.main()
