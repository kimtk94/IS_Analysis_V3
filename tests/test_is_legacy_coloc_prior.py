"""Regression tests for IS legacy coloc summary-prior diagnostics and overlay."""
import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/"scripts/is/plot_is_legacy_coloc_prior_sensitivity.py"
JOIN=ROOT/"scripts/is/integrate_is_p0_genetics_cell_evidence.py"
RSCRIPT=ROOT/"scripts/is/reweight_is_legacy_coloc_priors.R"


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


class PriorAuditTests(unittest.TestCase):
    def test_input_guard_for_integrated_gene_overlay(self):
        m=module(JOIN,"is_join")
        with tempfile.TemporaryDirectory() as t:
            d=Path(t)
            with self.assertRaises(FileNotFoundError):
                m.run(d/"a",d/"b",d/"c",d/"d",d/"out")
            out=d/"out"
            out.mkdir()
            (out/"history.txt").write_text("DONOTOVERWRITE")
            with self.assertRaises(FileExistsError):
                m.run(d/"a",d/"b",d/"c",d/"d",out)
            self.assertEqual((out/"history.txt").read_text(),"DONOTOVERWRITE")

    def test_prior_figure_from_synthetic_grid(self):
        m=module(FIG,"is_p12_plot")
        with tempfile.TemporaryDirectory() as t:
            d=Path(t)
            src=d/"grid.tsv"
            with src.open("w",newline="") as h:
                w=csv.writer(h,delimiter="\t")
                w.writerow(["gene_symbol","conditional_p12","PP_H4"])
                for gene in m.FOCUS:
                    for p in (1e-6,3e-6,1e-5,3e-5,1e-4):
                        w.writerow([gene,p,0.5])
            out=d/"figure.svg"
            m.run(src,out)
            ET.parse(out)
            self.assertGreater(out.stat().st_size,1000)
            with self.assertRaises(FileExistsError):
                m.run(src,out)

    @unittest.skipUnless(shutil.which("Rscript"),"Rscript absent")
    def test_r_coloc_posterior_reweight_646_fixture(self):
        # Test requires coloc 5.2.3. Local env opt-in for packages.
        env=os.environ.copy()
        if env.get("IS_COLOC_R_LIB"):
            env["R_LIBS_USER"]=env["IS_COLOC_R_LIB"]
        probe=subprocess.run(
            ["Rscript","--vanilla","-e",
             "cat(if (requireNamespace('coloc',quietly=TRUE) && as.character(packageVersion('coloc'))=='5.2.3') 'YES' else 'NO')"],
            capture_output=True,text=True,env=env,timeout=30
        )
        if probe.stdout.strip()!="YES":
            self.skipTest("coloc 5.2.3 not available in test R environment")
        with tempfile.TemporaryDirectory() as t:
            d=Path(t)
            src=d/"master.tsv"
            cols=["locus","dataset_key","gene_base","gene_symbol","status","nsnps",
                  "qtl_n","PP.H0","PP.H1","PP.H2","PP.H3","PP.H4"]
            with src.open("w",newline="") as h:
                w=csv.DictWriter(h,fieldnames=cols,delimiter="\t")
                w.writeheader()
                for i in range(646):
                    w.writerow(dict(locus="BBJ_IS_L001",dataset_key=f"TISSUE_{i}",
                                    gene_base=f"ENSG_{i:06d}",gene_symbol="FGF5",
                                    status="PASS",nsnps=1694,qtl_n=175,
                                    **{"PP.H0":0.01,"PP.H1":0.1,"PP.H2":0.01,
                                       "PP.H3":0.2,"PP.H4":0.68}))
            out=d/"out"
            proc=subprocess.run(["Rscript","--vanilla",str(RSCRIPT),str(src),str(out)],
                                env=env,text=True,capture_output=True,timeout=90)
            self.assertEqual(proc.returncode,0,msg=proc.stderr[-3000:])
            self.assertIn("GRID_ROWS=3230",proc.stdout)
            with (out/"IS_LEGACY_646_ABF_P12_GRID.tsv").open() as h:
                rec=list(csv.DictReader(h,delimiter="\t"))
            self.assertEqual(len(rec),3230)
            baseline=[r for r in rec if float(r["conditional_p12"])==1e-5]
            self.assertEqual(len(baseline),646)
            self.assertTrue(all(abs(float(r["PP_H4"])-0.68)<1e-10 for r in baseline))
            rerun=subprocess.run(["Rscript","--vanilla",str(RSCRIPT),str(src),str(out)],
                                env=env,text=True,capture_output=True,timeout=90)
            self.assertNotEqual(rerun.returncode,0)
            self.assertIn("refusing overwrite",rerun.stderr.lower())


if __name__=="__main__":
    unittest.main()
