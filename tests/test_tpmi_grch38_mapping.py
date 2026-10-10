import importlib.util
import sys
import unittest
from pathlib import Path

source=Path(__file__).resolve().parents[1]/"scripts/is/build_tpmi_grch38_mapping.py"
spec=importlib.util.spec_from_file_location("tpmi_build",source)
module=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=module
spec.loader.exec_module(module)

class TestTPMIReferenceMapping(unittest.TestCase):
    def setUp(self):
        self.target={"allele_verified_rsids":"rs671","chrom":"12","ref":"G","alt":"A"}
        self.variation={"name":"rs671","mappings":[
            {"assembly_name":"GRCh38","seq_region_name":"12",
             "start":111803962,"end":111803962,"strand":1,"allele_string":"G/A"}]}
        self.sequence={"id":"chromosome:GRCh38:12:111803962:111803962:1","seq":"G"}

    def test_rs671(self):
        status,vid=module.resolve(self.target,self.variation,self.sequence)
        self.assertEqual(status,"VERIFIED_GRCH38_REF_ALT")
        self.assertEqual(vid,"12:111803962:G:A")

    def test_reverse_reference(self):
        self.sequence["seq"]="A"
        status,vid=module.resolve(self.target,self.variation,self.sequence)
        self.assertEqual(vid,"12:111803962:A:G")

    def test_reject_alt_not_present(self):
        self.variation["mappings"][0]["allele_string"]="G/C"
        status,_=module.resolve(self.target,self.variation,self.sequence)
        self.assertEqual(status,"SOURCE_ALLELES_NOT_SUPPORTED_GRCH38")

    def test_reject_wrong_assembly(self):
        self.variation["mappings"][0]["assembly_name"]="GRCh37"
        status,_=module.resolve(self.target,self.variation,self.sequence)
        self.assertEqual(status,"MAPPING_NOT_UNIQUE")

    def test_reject_multimapping(self):
        self.variation["mappings"].append(dict(self.variation["mappings"][0]))
        status,_=module.resolve(self.target,self.variation,self.sequence)
        self.assertEqual(status,"MAPPING_NOT_UNIQUE")

if __name__=="__main__":
    unittest.main()
