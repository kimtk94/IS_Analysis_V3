"""Small self-contained matrix fixtures for original IS 646 source LD numeric auditor."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

try:
    import numpy as np
except ImportError:
    np=None

ROOT=Path(__file__).resolve().parents[1]
AUDIT=ROOT/"scripts/is/audit_is_legacy_646_ld_numerical.py"


@unittest.skipUnless(np is not None,"Numerical QC needs numpy")
class TestISReferenceLDDiagnostics(unittest.TestCase):
    def auditor(self):
        spec=importlib.util.spec_from_file_location("ISLegacyLDNumeric",AUDIT)
        module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def fixture(self,t):
        root=Path(t)
        input_dir=root/"inputs"
        ldroot=root/"reference"/"ld_v3"
        genotype=ldroot.parent/"ld"
        for d in (input_dir,ldroot,genotype):
            d.mkdir(parents=True,exist_ok=True)
        locus="BBJ_IS_L001"
        dataset="GTEx_V8__Brain_Cerebellar_Hemisphere"
        ensg="ENSG00000138675"
        file=input_dir/f"{locus}__{dataset}__{ensg}.tsv"
        fields=["locus","dataset_key","gene_base","variant_id","match_key",
                "gwas_maf","gwas_eaf","harmonization"]
        with file.open("w",encoding="utf8") as f:
            f.write("\t".join(fields)+"\n")
            for i in range(101):
                source=f"4:{10000+i}:A:C"
                dest=f"4:{20000+i}:A:C"
                f.write(f"{locus}\t{dataset}\t{ensg}\t{source}\t{dest}\t0.2\t0.2\tMATCHED\n")
        ids=[f"4:{10000+i}:A:C" for i in range(101)]
        (ldroot/f"{locus}.unphased.vcor1.bin.vars").write_text("\n".join(ids)+"\n")
        np.eye(101,dtype="<f8").tofile(ldroot/f"{locus}.unphased.vcor1.bin")
        (ldroot/f"{locus}.log").write_text("PLINK 2.0\n504 samples (founders)\n")
        (genotype/f"{locus}.matched.pvar").write_text(
            "#CHROM\tPOS\tID\tREF\tALT\n" +
            "".join(f"4\t{10000+i}\t{ids[i]}\tA\tC\n" for i in range(101)))
        (genotype/f"{locus}.matched.psam").write_text(
            "#FID\tIID\tSEX\n"+"".join(f"ID{i}\tID{i}\tNA\n" for i in range(504)))
        return input_dir,ldroot

    def test_basic_reference_matrix_qc_pass(self):
        m=self.auditor()
        with tempfile.TemporaryDirectory() as t:
            inp,ld=self.fixture(t)
            out=Path(t)/"report"
            v=m.run(inp,ld,out)
            self.assertEqual(v["source_file_count"],1)
            self.assertEqual(v["n_source_SNP_rows"],101)
            self.assertEqual(v["panels"][0]["sample_count"],504)
            self.assertEqual(v["panels"][0]["min_sampled_eigenvalue"],1.0)
            self.assertTrue((out/"IS_646_CACHE_LD_NUMERIC_AUDIT.json").exists())

    def test_nonempty_output_fails_closed(self):
        m=self.auditor()
        with tempfile.TemporaryDirectory() as t:
            inp,ld=self.fixture(t)
            out=Path(t)/"report"
            out.mkdir()
            sentinel=out/"preserve_original.txt"
            sentinel.write_text("PRESERVE")
            with self.assertRaises(FileExistsError):
                m.run(inp,ld,out)
            self.assertEqual(sentinel.read_text(),"PRESERVE")

    def test_invalid_correlation_matrix_does_not_generate_pass(self):
        m=self.auditor()
        with tempfile.TemporaryDirectory() as t:
            inp,ld=self.fixture(t)
            p=ld/"BBJ_IS_L001.unphased.vcor1.bin"
            mm=np.memmap(p,dtype="<f8",mode="r+",shape=(101,101))
            mm[0,1]=mm[1,0]=1.5
            mm.flush()
            out=Path(t)/"report"
            with self.assertRaises(ValueError):
                m.run(inp,ld,out)
            self.assertFalse(out.exists())


if __name__=="__main__":
    unittest.main()
