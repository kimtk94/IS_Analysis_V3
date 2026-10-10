"""Matrix regression: repaired REF/ALT IDs require original-source provenance."""
import importlib.util
import sys
import unittest
from pathlib import Path

PATH=(Path(__file__).resolve().parents[1]/
      "scripts/is/build_alcohol_ais_verified_variant_matrix.py")
spec=importlib.util.spec_from_file_location("alcohol_ais_matrix",PATH)
module=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=module
spec.loader.exec_module(module)


def fixtures():
    exposure=[]
    provenance=[]
    overlap=[]
    for i in range(11):
        vid=f"4:{100+i}:T:C"
        missing=i==10
        prov_status=("MISSING_BOTH" if missing else
                     "CANONICAL_ID_IS_OTHER_EFFECT_NOT_REF_ALT" if i==0 else
                     "PASS_FASTA_REF_ALT")
        exposure.append({"ID":vid,"reference_ALT":"C","source_beta_ALT":"0.2",
                         "source_se":"0.01","source_original_p":"1e-10"})
        provenance.append({"variant_id_fasta":vid,"status":prov_status,
                           "source_md5_pass":"True","source_canonical_stats_equal":"True",
                           "validated_fasta_variant_id":vid,
                           "original_canonical_id":f"4:{100+i}:C:T" if i==0 else vid,
                           "ais_beta_alt":"-0.3","ais_se":"0.1",
                           "ais_p":"0.02","ais_eaf_alt":"0.25"})
        for ds in ("BBJ_JAPAN_IS","GIGASTROKE_EAS_AIS"):
            status=("MISSING_FROM_DATASET" if missing else
                    "REF_ALT_SWAP_REVIEW" if ds=="GIGASTROKE_EAS_AIS" and i==0 else
                    "MATCH_EXACT")
            overlap.append({"variant_id":vid,"dataset":ds,"locus":"ADH1B",
                            "rsid_verified":f"rs{i}","status":status,
                            "beta_alt":"-0.3" if not missing else "",
                            "se":"0.1" if not missing else "",
                            "p":"0.02" if not missing else "",
                            "eaf_alt":"0.25" if not missing else "",
                            "source_variant_id":vid})
    return overlap,provenance,exposure


class TestSourceAuditedMatrix(unittest.TestCase):
    def test_repaired_legacy_id_is_verified(self):
        overlap,prov,expo=fixtures()
        rows=module.make_matrix(overlap,prov,expo)
        self.assertEqual(len(rows),22)
        first=[r for r in rows if r["rsid"]=="rs0" and
               r["AIS_source"]=="GIGASTROKE_EAS_AIS"][0]
        self.assertEqual(first["source_quality_status"],
                         "VERIFIED_SOURCE_FASTA_ID_REPAIRED_IN_DERIVATIVE")
        self.assertEqual(first["variant_id_GRCh37_REF_ALT"],"4:100:T:C")
        self.assertEqual(first["ais_beta_ALT"],"-0.3")

    def test_fail_on_raw_beta_discordance(self):
        overlap,prov,expo=fixtures()
        prov[0]["ais_beta_alt"]="0.99"
        with self.assertRaises(ValueError):
            module.make_matrix(overlap,prov,expo)

    def test_fail_on_missing_source_audit(self):
        overlap,prov,expo=fixtures()
        prov[0]["source_canonical_stats_equal"]="False"
        with self.assertRaises(ValueError):
            module.make_matrix(overlap,prov,expo)

    def test_fail_on_non_11_unique_inputs(self):
        overlap,prov,expo=fixtures()
        with self.assertRaises(ValueError):
            module.make_matrix(overlap[:-1],prov,expo)


if __name__=="__main__":
    unittest.main()
