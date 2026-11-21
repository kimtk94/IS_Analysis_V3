"""Offline tests for deterministic independent GRCh37 44-sentinel allele gate."""
from __future__ import annotations
import hashlib
import io
import json
import unittest
from scripts.ckm_stroke.audit_bmi_grch37_reference_sentinels import (
    classify_reference_base,pick_sentinels,inspect_batch,ENSEMBL_URL,
)


def make_snp(i,ch):
    return f"rs{ch*100+i}"


def fake_rows(chromosomes=(1,2),n=5):
    for ch in chromosomes:
        for i in range(n):
            yield (i+1, [str(ch),str(1000+10*i),make_snp(i,ch),"A","C",
                         "0.2","0.1","0.03","0.01","600000"])
    yield (0,["__SOURCE_COMPLETE_SHA256_PASS__","46754937"])


def fake_candidate(chromosomes=(1,2),n=5):
    return {make_snp(i,ch):("A","C") for ch in chromosomes for i in range(n)}


class Response(io.BytesIO):
    def __init__(self, data):
        super().__init__(data)
        self.status=200
    def geturl(self):
        return ENSEMBL_URL


def valid_provider(request,timeout):
    self=valid_provider
    self.last_request=request
    data=json.loads(request.data)
    output=[]
    for x in data["regions"]:
        output.append({
            "query":x,
            "id":"chromosome:GRCh37:"+x.split(":")[0]+":"+x.split(":")[1].split("..")[0]+":"+x.split(":")[1].split("..")[-1]+":1",
            "seq":"A",
        })
    return Response(json.dumps(output).encode())


class SentinelTests(unittest.TestCase):
    def test_allele_reference_classifications(self):
        self.assertEqual(classify_reference_base("A","C","A"),"reference_base_in_forward_allele_pair")
        self.assertEqual(classify_reference_base("A","C","G"),"reverse_complement_only")
        self.assertEqual(classify_reference_base("A","C","T"),"reverse_complement_only")
        self.assertEqual(classify_reference_base("A","C","N"),"unknown_or_nonacgt_reference")

    def test_sampling_deterministic_and_per_chromosome(self):
        cand=fake_candidate()
        s,meta=pick_sentinels(cand,fake_rows(),chromosomes=(1,2),per_chr=2,strict_original=False)
        self.assertEqual(len(s),4)
        self.assertEqual({x["chr"] for x in s},{1,2})
        self.assertEqual(s,[*pick_sentinels(cand,fake_rows(),chromosomes=(1,2),per_chr=2,strict_original=False)[0]])
        self.assertEqual(meta["matched_candidate_rsids_in_original"],10)

    def test_missing_full_source_sha_rejected(self):
        rows=list(fake_rows())
        rows.pop()
        with self.assertRaisesRegex(ValueError,"sentinel"):
            pick_sentinels(fake_candidate(),iter(rows),chromosomes=(1,2),per_chr=2,strict_original=False)

    def test_discordant_candidate_alleles_rejected(self):
        candidates=fake_candidate()
        candidates["rs100"]=("G","C")
        with self.assertRaisesRegex(ValueError,"allele orientation differs"):
            pick_sentinels(candidates,fake_rows(),chromosomes=(1,2),per_chr=2,strict_original=False)

    def test_44_region_api_batch_and_forward_bases(self):
        panel=[{
            "chr":ch,"pos":1000+i,"rsid":f"rs{ch*10+i}","A1":"A","A2":"C",
            "rank_sha256":hashlib.sha256(str((ch,i)).encode()).hexdigest()
        } for ch in range(1,23) for i in range(2)]
        rows,counts=inspect_batch(panel,provider=valid_provider)
        self.assertEqual(len(rows),44)
        self.assertEqual(counts["reference_base_in_forward_allele_pair"],44)
        self.assertEqual(len(json.loads(valid_provider.last_request.data)["regions"]),44)
        self.assertEqual(valid_provider.last_request.get_method(),"POST")

    def test_api_reference_wrong_build_rejected(self):
        panel=[{"chr":1,"pos":1000,"A1":"A","A2":"C","rsid":"rs10","rank_sha256":"0"*64}]
        def bad_provider(req,timeout):
            return Response(json.dumps([{"query":"1:1000..1000:1",
                "id":"chromosome:GRCh38:1:1000:1000:1","seq":"A"}]).encode())
        with self.assertRaisesRegex(ValueError,"GRCh37"):
            inspect_batch(panel,provider=bad_provider)

    def test_duplicate_sentinel_position_rejected(self):
        a={"chr":1,"pos":1000,"A1":"A","A2":"C","rsid":"rs1","rank_sha256":"0"*64}
        with self.assertRaisesRegex(ValueError,"duplicate sentinel"):
            inspect_batch([a,dict(a)],provider=valid_provider)

    def test_unexpected_allele_does_not_pass(self):
        panel=[{"chr":1,"pos":1000,"A1":"A","A2":"C","rsid":"rs10","rank_sha256":"0"*64}]
        def provider(req,timeout):
            return Response(json.dumps([{"query":"1:1000..1000:1",
                 "id":"chromosome:GRCh37:1:1000:1000:1","seq":"G"}]).encode())
        rows,counts=inspect_batch(panel,provider=provider)
        self.assertEqual(counts["reverse_complement_only"],1)


if __name__=="__main__":
    unittest.main()
