import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

FILE=Path(__file__).resolve().parents[1]/"scripts/is/audit_gtex_v8_gencode26_cis_geometry.py"
SPEC=importlib.util.spec_from_file_location("gtex_tss",FILE)
m=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)

def row(name,gene,chrom,start,end,strand):
    return "\t".join(["chr"+chrom,"HAVANA","gene",str(start),str(end),
                      ". ",strand,".",f'gene_id "{gene}.11"; gene_name "{name}";']).replace("\t. \t","\t.\t")

GOOD={
 "FGF5":("ENSG00000138675","4",79800000,79900000,"+"),
 "C4orf22":("ENSG00000197826","4",80000000,80500000,"-"),
 "CALHM2":("ENSG00000138172","10",103000000,103200000,"+"),
 "NEURL":("ENSG00000107954","10",103400000,103500000,"-"),
}
def fixture(rows):
    with tempfile.NamedTemporaryFile(suffix=".gtf.gz") as f:
        body="##gencode version 26\n"+"\n".join(
            row(name,*parts) for name,parts in rows.items())+"\n"
        f.write(gzip.compress(body.encode(),mtime=0))
        f.flush()
        return m.inspect(f.name)

class CisGeometryTests(unittest.TestCase):
    def test_complete_mock_tss_geometry_bounded(self):
        a=fixture(GOOD)
        self.assertTrue(a["four_gene_TSS_windows_covered"])
        self.assertIn("PROXY_INSIDE",a["geometry_gate"])
        self.assertFalse(a["GTEx_original_collapsed_gene_model_independently_verified"])
        self.assertFalse(a["source_specific_GTEx_original_n_variants_independently_verified"])
        self.assertFalse(a["coloc_performed"])

    def test_reverse_strand_uses_gene_end_not_start(self):
        a=fixture(GOOD)
        self.assertEqual(a["gene_boundary_screen"]["C4orf22"]["TSS_GRCh38"],80500000)
        self.assertEqual(a["gene_boundary_screen"]["FGF5"]["TSS_GRCh38"],79800000)

    def test_truncated_interval_rejected(self):
        x=dict(GOOD)
        x["FGF5"]=("ENSG00000138675","4",81350000,81400000,"+")
        result=fixture(x)
        self.assertFalse(result["four_gene_TSS_windows_covered"])
        self.assertIn("TRUNCATES",result["geometry_gate"])

    def test_missing_gene_is_error(self):
        x=dict(GOOD); x.pop("NEURL")
        with self.assertRaises(ValueError):
            fixture(x)

    def test_gene_id_and_chr_gate(self):
        x=dict(GOOD)
        x["CALHM2"]=("ENSG00000138172","4",103000000,103200000,"+")
        result=fixture(x)
        self.assertFalse(result["gene_boundary_screen"]["CALHM2"]["both_ends_inside_predeclared_interval"])

    def test_oversized_file_fails_closed(self):
        with tempfile.NamedTemporaryFile(suffix=".gtf.gz") as f:
            f.write(b"not gzip");f.flush()
            with self.assertRaises(OSError):
                m.inspect(f.name)

    def test_ref_and_magnitude_retained(self):
        z=fixture(GOOD)
        self.assertEqual(z["gene_boundary_screen"]["FGF5"]["GENCODE_v26_gene_TSS_plusminus1Mb_start"],78800000)
        self.assertEqual(z["gene_boundary_screen"]["FGF5"]["previous_screen_region"],"4:79000000-81500000")
        self.assertTrue(z["scientific_status"]=="HOLD")

if __name__=="__main__":
    unittest.main()
