import importlib.util
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

source=Path(__file__).resolve().parents[1]/"scripts/is/audit_tpmi_grch38_chain_ucsc.py"
spec=importlib.util.spec_from_file_location("tpmi_chain_audit",source)
audit=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=audit
spec.loader.exec_module(audit)

class FakeLiftOver:
    def convert_coordinate(self,chrom,pos):
        return [(chrom,pos+100,"+",0)]

class TestChainAndReference(unittest.TestCase):
    def test_chain_uses_zero_based_coordinate(self):
        target={"variant_id":"12:10:G:A","chrom":"12","pos":"10"}
        out=audit.intervals_for_targets([target],FakeLiftOver())
        self.assertEqual(out[0]["pos38"],110)
        self.assertEqual(out[0]["chr38"],"chr12")

    def test_ref_fetch_position(self):
        rows=[{"chr38":"chr12","pos38":101},{"chr38":"chr12","pos38":103}]
        mock={"genome":"hg38","chrom":"chr12","start":100,"end":103,"dna":"GTA"}
        with patch.object(audit,"get_json",return_value=mock):
            got=audit.fetch_reference_bases(rows)
        self.assertEqual(got[("chr12",101)],"G")
        self.assertEqual(got[("chr12",103)],"A")

    def test_ensembl_rsid_matches(self):
        row={"chrom":"12","pos38":111803962,"allele_verified_rsids":"rs671",
             "ref":"G","alt":"A"}
        variation={"name":"rs671","mappings":[
            {"assembly_name":"GRCh38","seq_region_name":"12",
             "start":111803962,"end":111803962,
             "strand":1,"allele_string":"G/A"}]}
        self.assertEqual(audit.check_ensembl_rs(row,variation),
                         "ENSEMBL_RSID_POSITION_ALLELES_PASS")
        variation["mappings"][0]["end"]=111803963
        self.assertEqual(audit.check_ensembl_rs(row,variation),
                         "ENSEMBL_MAPPING_MISMATCH_OR_DUPLICATE")

if __name__=="__main__":
    unittest.main()
