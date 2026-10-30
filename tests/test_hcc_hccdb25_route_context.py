"""Offline controls for HCCDB Vue route module discovery."""
import pathlib
import sys
import unittest
from unittest.mock import patch
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
import probe_hccdb25_route_context as p

class RouteTests(unittest.TestCase):
    def test_public_static_route_chunk_is_structure_not_rna(self):
        raw=b'const routes=[{path:"/repository",component:()=>n.e("chunk-9b202").then(n.bind(null,"7199"))},{path:"/download",name:"Download"}]'
        obs=p.context_windows(raw)
        self.assertTrue(any(x["marker"]=="/repository" for x in obs))
        self.assertTrue(any("chunk-9b202" in x["nearby_chunk_names"] for x in obs))
        self.assertTrue(all(not x["context_has_direct_data_file_extension"] for x in obs))
        self.assertNotIn("7199",str(obs))

    def test_empty_script_yields_no_module(self):
        self.assertEqual(p.context_windows(b"window.__my_app=true"),[])

    def test_route_probe_never_claims_quant_matrix(self):
        with patch.object(p,"small_get",return_value=(b'path:"/repository"',{"status":206,"truncated":False})):
            r=p.probe()
        self.assertEqual(r["status"],"PUBLIC_ROUTE_MODULE_INSPECTED")
        self.assertFalse(r["hccdb25_expression_file_verified"])
        self.assertFalse(r["public_matrix_download_performed"])
        self.assertFalse(r["pdc_case_join_verified"])

if __name__=="__main__":
    unittest.main()
