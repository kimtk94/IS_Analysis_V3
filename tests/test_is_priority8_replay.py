"""Regression: eight real-phase source workflows, with entirely synthetic test input."""
import csv
import importlib.util
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/is/audit_is_priority8_snp_vs_summary_prior_grid.py"
RSCRIPTS=(
    ROOT/"scripts/is/replay_is_phase10_abf_priority8_snp_inputs.R",
    ROOT/"scripts/is/audit_is_legacy_crossancestry_maf_priority8.R",
)


def load():
    sp=importlib.util.spec_from_file_location("is_priority8",SCRIPT)
    mod=importlib.util.module_from_spec(sp)
    sp.loader.exec_module(mod)
    return mod


def fixture(folder):
    m=load()
    keys=list(m.GENE_ID.items())
    pri=folder/"priors.tsv"
    rer=folder/"replay.tsv"
    cols=["locus","dataset_key","gene_base","conditional_p12",
          "PP_H0","PP_H1","PP_H2","PP_H3","PP_H4"]
    rcols=["gene","locus","tissue","p12","qtl_n_mode",
           "PP_H0","PP_H1","PP_H2","PP_H3","PP_H4"]
    grid=[1e-6,3e-6,1e-5,3e-5,1e-4]
    def h(p):
        # Synthetic valid posterior vector (does not assert any real-data biological fit).
        v=.3+(p/1e-4)*.12
        return [0.1,0.2,0.1,0.6-v,v]
    with pri.open("w",newline="") as fp,rer.open("w",newline="") as fr:
        wp=csv.DictWriter(fp,fieldnames=cols,delimiter="\t")
        wr=csv.DictWriter(fr,fieldnames=rcols,delimiter="\t")
        wp.writeheader();wr.writeheader()
        for j in range(646):
            if j<8:
                gene,ens=keys[j]
                locus,tissue=f"LOCUS_{gene}",f"GTEx_TISSUE_{gene}"
            else:
                gene,ens="FILLER",f"ENSG_FILLER_{j}"
                locus,tissue=f"LOCUS_FILLER_{j}",f"TISSUE_FILLER_{j}"
            for p in grid:
                d=dict(zip(["PP_H0","PP_H1","PP_H2","PP_H3","PP_H4"],h(p)))
                wp.writerow(dict(locus=locus,dataset_key=tissue,gene_base=ens,conditional_p12=p,**d))
                if j<8:
                    for mode in ("FIRST","MIN","MAX"):
                        wr.writerow(dict(gene=gene,locus=locus,tissue=tissue,p12=p,qtl_n_mode=mode,**d))
    return pri,rer


class TestPriority8(unittest.TestCase):
    def test_full_eight_by_five_grid(self):
        m=load()
        self.assertEqual(len(m.GENE_ID),8)
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);a,b=fixture(p)
            res=m.run(a,b,p/"audit.json")
            self.assertEqual(res["paired_analyses"],40)
            self.assertEqual(res["pairs"],8)
            self.assertLess(res["max_abs_posterior_delta"],1e-12)

    def test_reject_inconsistent_original_baseline_h4(self):
        m=load()
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);a,b=fixture(p)
            lines=b.read_text().splitlines()
            headers=lines[0].split("\t")
            h4=headers.index("PP_H4")
            changed=False
            for i in range(1,len(lines)):
                cells=lines[i].split("\t")
                if cells[headers.index("qtl_n_mode")]=="FIRST" and float(cells[headers.index("p12")])==1e-5:
                    cells[h4]="0.900"
                    lines[i]="\t".join(cells)
                    changed=True
                    break
            assert changed
            content="\n".join(lines)+"\n"
            b.write_text(content)
            with self.assertRaises(ValueError):
                m.run(a,b,p/"result.json")

    def test_refuse_nonempty_output(self):
        m=load()
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);a,b=fixture(p)
            out=p/"keep.json";out.write_text("IMPORTANT")
            with self.assertRaises(FileExistsError):
                m.run(a,b,out)
            self.assertEqual(out.read_text(),"IMPORTANT")

    def test_summary_ledger_refuses_existing_output(self):
        path=ROOT/"scripts/is/summarize_is_priority8_snp_evidence.py"
        sp=importlib.util.spec_from_file_location("is_priority8_summary",path)
        module=importlib.util.module_from_spec(sp)
        sp.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as t:
            d=Path(t)
            output=d/"legacy.json"
            output.write_text("UNCHANGED")
            with self.assertRaises(FileExistsError):
                module.run(d/"r",d/"q",d/"m",d/"a",output)
            self.assertEqual(output.read_text(),"UNCHANGED")

    @unittest.skipUnless(shutil.which("Rscript"),"Rscript absent")
    def test_priority8_r_source_parses(self):
        for s in RSCRIPTS:
            proc=subprocess.run(
                ["Rscript","--vanilla","-e",
                 "invisible(parse(file=commandArgs(TRUE)[1]));cat('R_PARSE_OK\\n')",str(s)],
                capture_output=True,text=True,timeout=25)
            self.assertEqual(proc.returncode,0,msg=proc.stderr)
            self.assertIn("R_PARSE_OK",proc.stdout)


if __name__=="__main__":
    unittest.main()
