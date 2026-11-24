import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

FILE=Path(__file__).resolve().parents[1]/"scripts/is/audit_gtex_v8_api_collapsed_exons.py"
spec=importlib.util.spec_from_file_location("gtex_api_audit",FILE)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

IDS={
 "FGF5":{"gencodeId":"ENSG00000138675.14","chromosome":"4","strand":"+","regular_TSS":80266599},
 "NEURL":{"gencodeId":"ENSG00000107954.12","chromosome":"10","strand":"+","regular_TSS":103493979},
 "C4orf22":{"gencodeId":"ENSG00000197826.10","chromosome":"4","strand":"+","regular_TSS":80335720},
 "CALHM2":{"gencodeId":"ENSG00000138172.13","chromosome":"10","strand":"-","regular_TSS":103452402},
}
def example(url):
    if "datasetId=gtex_v8" not in url:raise ValueError("wrong release")
    if "00138675" in url:
        item={"start":80266599,"end":80290000,"chromosome":"chr4","strand":"+"}
    elif "00107954" in url:
        item={"start":103493979,"end":103510000,"chromosome":"chr10","strand":"+"}
    elif "00197826" in url:
        item={"start":80335720,"end":80400000,"chromosome":"chr4","strand":"+"}
    elif "00138172" in url:
        item={"start":103440000,"end":103452402,"chromosome":"chr10","strand":"-"}
    else:raise ValueError("wrong ID")
    return json.dumps({"data":[item],"paging_info":{"total":1}}).encode()

class TestOfficialAPI(unittest.TestCase):
    def test_gene_id_version_required(self):
        with self.assertRaises(ValueError):
            m.api_call("ENSG00000138675",fetcher=example)

    def test_request_v8_not_default_v10(self):
        q=m.api_call("ENSG00000138675.14",fetcher=example)
        self.assertIn("datasetId=gtex_v8",q["url"])
        self.assertEqual(q["response"]["data"][0]["start"],80266599)

    def test_api_four_genes_positive_and_outside(self):
        o=m.audit(ids=IDS,fetcher=example)
        self.assertTrue(o["all_four_gene_exon_models_parsed"])
        self.assertTrue(o["two_missing_QTD_variants_outside_API_collapsed_TSS"])
        self.assertEqual(o["by_gene"]["FGF5"]["outside_bp"],19)
        self.assertEqual(o["by_gene"]["NEURL"]["outside_bp"],229)
        self.assertFalse(o["FastQTL_original_tested_gene_BED_TSS_verified"])
        self.assertEqual(o["scientific_status"],"HOLD")

    def test_missing_strand_is_hold(self):
        def no_strand(url):
            result=json.loads(example(url))
            for item in result["data"]:item.pop("strand")
            return json.dumps(result).encode()
        o=m.audit(ids=IDS,fetcher=no_strand)
        self.assertFalse(o["all_four_gene_exon_models_parsed"])
        self.assertIn("HOLD",o["api_gate"])

    def test_wrong_chr_is_hold(self):
        def bad(url):
            d=json.loads(example(url))
            if "00138675" in url:d["data"][0]["chromosome"]="chr17"
            return json.dumps(d).encode()
        o=m.audit(ids=IDS,fetcher=bad)
        self.assertFalse(o["all_four_gene_exon_models_parsed"])

    def test_inside_cis_boundary_does_not_pass(self):
        def changed(url):
            d=json.loads(example(url))
            if "00138675" in url:d["data"][0]["start"]=80300000
            return json.dumps(d).encode()
        o=m.audit(ids=IDS,fetcher=changed)
        self.assertFalse(o["two_missing_QTD_variants_outside_API_collapsed_TSS"])
        self.assertEqual(o["api_gate"],"HOLD_API_SOURCE_INCOMPLETE_OR_DIFFERENT_BOUNDARIES")

    def test_api_failure_does_not_imply_original_qtl_absence(self):
        def fail(url):raise OSError("network")
        o=m.audit(ids=IDS,fetcher=fail)
        self.assertFalse(o["all_four_gene_exon_models_parsed"])
        self.assertFalse(o["source_variant_missingness_proven"])
        self.assertFalse(o["coloc_performed"])

    def test_gtf_extracts_true_versioned_gene_ids(self):
        lines=[]
        for name,vals in IDS.items():
            gene=vals["gencodeId"]
            chrom=vals["chromosome"]
            strand=vals["strand"]
            pos=vals["regular_TSS"]
            start,end=(pos,pos+100) if strand=="+" else (pos-100,pos)
            lines.append("\t".join(["chr"+chrom,"TEST","gene",str(start),str(end),".",strand,".",
                        'gene_id "'+gene+'"; gene_name "'+name+'";']))
        with tempfile.NamedTemporaryFile(suffix=".gz") as f:
            f.write(gzip.compress(("\n".join(lines)+"\n").encode()))
            f.flush()
            ids=m.versioned_ids(f.name)
        self.assertEqual(ids["FGF5"]["gencodeId"],IDS["FGF5"]["gencodeId"])
        self.assertEqual(ids["CALHM2"]["strand"],"-")

    def test_payload_size_cap(self):
        with self.assertRaises(ValueError):
            m.api_call("ENSG00000138675.14",fetcher=lambda _:b"0"*(m.MAX_API+1))

if __name__=="__main__":
    unittest.main()
