"""Tests for 646 original GTEx SNP-varying QTL-N common-sdY calibration.

The science contract: N can enter coloc's sdY.est via a vector, but does not
replace the supplied per-SNP beta/varbeta or substitute matched QTL cohort LD.
"""
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
RUNNER=ROOT/"scripts/is/is_646_qtl_variant_n_sdy_calibration.R"
BUILDER=ROOT/"scripts/is/build_is_646_variant_n_sdy_colab.py"
NOTEBOOK=ROOT/"notebooks/is/IS_PHASE10_11_646_VARIANT_N_SDY_ONECELL.ipynb"
TEMPLATE=ROOT/"notebooks/is/IS_PHASE10_11_646_QTL_N_SENSITIVITY_ONECELL.ipynb"
EXPECTED_CATALOG="d1e237d3c0ad8281adf3d9fe78f69091505eeeab6d152ae0bc54edd5e964853f"

def load_builder():
    spec=importlib.util.spec_from_file_location("source_builder_sdy",BUILDER)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

class TestISVariantNCommonSDY(unittest.TestCase):
    def test_sdY_runner_R_parses(self):
        r=subprocess.run(["Rscript","--vanilla","-e",
            "invisible(parse(file=commandArgs(trailingOnly=TRUE)[1]));cat('OK_SDY_SCRIPT\\n')",
            str(RUNNER)],capture_output=True,text=True,timeout=15)
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertIn("OK_SDY_SCRIPT",r.stdout)

    def test_estimator_statistical_contract_on_synthetic_case(self):
        code=r"""
          suppressPackageStartupMessages(library(coloc))
          f=getFromNamespace("sdY.est","coloc")
          n=c(90,110,120,140,150,155,160,170,180,200)
          maf=seq(.08,.45,length.out=length(n))
          ysd=1.2
          V=ysd^2/(2*n*maf*(1-maf))
          a=suppressWarnings(f(V,maf,n))
          b=suppressWarnings(f(V,maf,rep(max(n),length(n))))
          stopifnot(is.finite(a),abs(a-ysd)<1e-10,abs(b-a)>1e-3)
          dat=list(type="quant",beta=rep(.05,length(n)),varbeta=V,sdY=a,
                   snp=paste0("rs",seq_along(n)))
          dat1=list(type="quant",beta=rep(.05,length(n)),varbeta=V,sdY=b,
                   snp=paste0("rs",seq_along(n)))
          cc=list(type="cc",beta=rep(.02,length(n)),varbeta=rep(.003,length(n)),
                  snp=dat$snp,N=174686,s=22664/174686)
          abf=suppressMessages(coloc.abf(cc,dat,p1=1e-4,p2=1e-4,p12=1e-5))
          old=suppressMessages(coloc.abf(cc,dat1,p1=1e-4,p2=1e-4,p12=1e-5))
          stopifnot(all(is.finite(abf$summary[grep("PP.H",names(abf$summary))])))
          cat("VECTOR_N_SDY_FORMULA_VALIDATED\n")
        """
        env=os.environ.copy()
        env.setdefault("R_LIBS_USER","/srv/is-analysis/.Rlib")
        r=subprocess.run(["Rscript","--vanilla","-e",code],capture_output=True,
                         text=True,timeout=20,env=env)
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertIn("VECTOR_N_SDY_FORMULA_VALIDATED",r.stdout)

    def test_single_cell_output_and_hash_gates(self):
        note=json.loads(NOTEBOOK.read_text())
        self.assertEqual(len(note["cells"]),1)
        code="".join(note["cells"][0]["source"])
        ast.parse(code)
        self.assertIn(repr(RUNNER.read_text()),code)
        self.assertIn(EXPECTED_CATALOG,code)
        self.assertIn("GWAS_SOURCE_EFFECT_DISAGREES_ACROSS_GTEX_ASSAYS",code)
        self.assertIn("FULL_646_GWAS_NATIVE_ALLELE_CATALOG_VERIFIED_BY_EXACT_SET_AND_EFFECTS",code)
        self.assertIn("IS_PHASE10_11_646_VARIANT_N_SDY_CALIBRATION_V1",code)
        self.assertIn("IS_VARIANT_N_SDY_CALIBRATION_SUMMARY.json",code)
        self.assertIn("IS_646_VARIANT_N_SDY_COLAB_MANIFEST.json",code)
        self.assertNotIn("results/IS_PHASE10_11_646_QTL_N_SCALAR_SENSITIVITY_V1",code)
        self.assertIn("PER_VARIANT_QTL_N_ONLY_IN_COMMON_SDY_ESTIMATE_NO_VARIANT_SE_RECALC",code)
        self.assertIn("if OUTPUT.exists()",code)

    def test_builder_reproducible_and_no_overwrite(self):
        m=load_builder()
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/NOTEBOOK.name
            m.build(TEMPLATE,RUNNER,p)
            self.assertEqual(
              hashlib.sha256(p.read_bytes()).hexdigest(),
              hashlib.sha256(NOTEBOOK.read_bytes()).hexdigest())
            with self.assertRaises(FileExistsError):
                m.build(TEMPLATE,RUNNER,p)
            self.assertTrue(p.exists())

    def test_scientific_gate_script_does_not_recalculate_se(self):
        s=RUNNER.read_text()
        self.assertIn('getFromNamespace("sdY.est","coloc")',s)
        self.assertIn("sdy_est(v,maf,n)",s)
        self.assertIn("sdy_est(v,maf,n[1])",s)
        self.assertIn('e$sdY<-yv',s)
        self.assertIn("BASELINE_SDY_EXPLICIT_REPRODUCTION_FAILURE",s)
        self.assertIn("NOT_ESTABLISHED",s)
        self.assertIn("digits=16",s)
        self.assertNotIn("varbeta=v*n",s)

if __name__=="__main__":
    unittest.main()
