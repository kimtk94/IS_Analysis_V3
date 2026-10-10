"""Unit/regression tests for FASTA-anchored GWAS source allele provenance."""
import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path

SOURCE=(Path(__file__).resolve().parents[1]/
        "scripts/is/audit_gigastroke_ais_refalt_provenance.py")
spec=importlib.util.spec_from_file_location("gigastroke_refalt_prov",SOURCE)
audit=importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)

class TestAlleleProvenance(unittest.TestCase):
    def test_effect_is_reference(self):
        self.assertEqual(audit.normalized_alleles("T","C","T"),("T","C",-1))

    def test_effect_is_alternative(self):
        self.assertEqual(audit.normalized_alleles("C","T","T"),("T","C",1))

    def test_reject_missing_fasta_allele(self):
        self.assertIsNone(audit.normalized_alleles("C","G","T"))

    def test_reject_non_snv(self):
        self.assertIsNone(audit.normalized_alleles("C","CT","C"))

    def test_fai_base_lookup(self):
        with tempfile.TemporaryDirectory() as td:
            fasta=Path(td)/"toy.fa"
            fasta.write_bytes(b">4\nTACG\n>12\nGTAC\n")
            fasta.with_suffix(".fa.fai").write_text("4\t4\t3\t4\t5\n12\t4\t12\t4\t5\n")
            idx=audit.indexed_fasta(fasta)
            with fasta.open("rb") as f:
                self.assertEqual(audit.ref_at(f,idx,"4",1),"T")
                self.assertEqual(audit.ref_at(f,idx,"12",2),"T")

    def test_gzip_awk_scan_targets(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"toy.tsv.gz"
            with gzip.open(p,"wt") as f:
                f.write("chr\tpos\teffect_allele\tother_allele\n")
                f.write("4\t1\tC\tT\n4\t2\tA\tT\n12\t3\tT\tG\n")
            rows=audit.awk_scan(p,1,2,{("4",1)},sample_stride=0)
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]["pos"],"1")

if __name__=="__main__":unittest.main()
