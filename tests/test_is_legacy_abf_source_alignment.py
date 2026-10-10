"""Source-provenance and SNP-replay gates for original IS ABF coloc."""
import csv
import importlib.util
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/is/audit_is_legacy_abf_source_alignment.py"

def mod():
    sp=importlib.util.spec_from_file_location("is_abf_alignment",SCRIPT)
    m=importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m

def fixtures(folder):
    idx=folder/"i.tsv"
    master=folder/"m.tsv"
    ann=folder/"a.tsv"
    keys=("locus","dataset_key","gene_base")
    post=["PP.H0","PP.H1","PP.H2","PP.H3","PP.H4"]
    with idx.open("w",newline="") as hi,master.open("w",newline="") as hm,ann.open("w",newline="") as ha:
        wi=csv.DictWriter(hi,delimiter="\t",fieldnames=[*keys,"status","nsnps","qtl_n_first","qtl_n_min","qtl_n_max"])
        wm=csv.DictWriter(hm,delimiter="\t",fieldnames=[*keys,"status","nsnps","qtl_n",*post])
        wa=csv.DictWriter(ha,delimiter="\t",fieldnames=[*keys,"status","nsnps","qtl_n",*post,"tissue_class"])
        wi.writeheader();wm.writeheader();wa.writeheader()
        for j in range(646):
            key=dict(locus=f"L{j%4}",dataset_key=f"D{j%16}",gene_base=f"ENSG{j}")
            wi.writerow(dict(key,status="READY",nsnps=1500,qtl_n_first=175,qtl_n_min=170,qtl_n_max=175))
            row=dict(key,status="PASS",nsnps=1500,qtl_n=175,**dict(zip(post,[.1,.2,.1,.2,.4])))
            wm.writerow(row);wa.writerow(dict(row,tissue_class="brain"))
    return idx,master,ann

class TestABFSourceProvenance(unittest.TestCase):
    def test_646_exact_key_and_hypothesis_match(self):
        m=mod()
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)
            i,master,ann=fixtures(p)
            result=m.run(i,master,ann,p/"audit.json")
            self.assertEqual(result["rows"],646)
            self.assertEqual(result["matching_H0_to_H4_annotation"],646)
            self.assertEqual(result["pairs_with_variable_snp_specific_QTL_N"],646)
            self.assertEqual(result["N_range_median"],5)

    def test_source_modified_h4_is_rejected(self):
        m=mod()
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);i,master,ann=fixtures(p)
            contents=ann.read_text()
            contents=contents.replace("\t0.4\tbrain\n","\t0.41\tbrain\n",1)
            ann.write_text(contents)
            with self.assertRaises(ValueError):
                m.run(i,master,ann,p/"result.json")

    def test_fail_closed_on_existing_output(self):
        m=mod()
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);i,master,ann=fixtures(p)
            old=p/"do_not_overwrite"
            old.write_text("OLD")
            with self.assertRaises(FileExistsError):
                m.run(i,master,ann,old)
            self.assertEqual(old.read_text(),"OLD")

    @unittest.skipUnless(shutil.which("Rscript"),"R not installed")
    def test_legacy_snp_replay_r_syntax(self):
        for name in ("replay_is_phase10_abf_snp_inputs.R","audit_is_legacy_crossancestry_maf.R"):
            p=ROOT/"scripts/is"/name
            proc=subprocess.run(["Rscript","--vanilla","-e",
                "invisible(parse(file=commandArgs(trailingOnly=TRUE)[1]));cat('R_PARSE_OK\\n')",
                str(p)],capture_output=True,text=True,timeout=20)
            self.assertEqual(proc.returncode,0,msg=proc.stderr)
            self.assertIn("R_PARSE_OK",proc.stdout)


if __name__=="__main__":
    unittest.main()
