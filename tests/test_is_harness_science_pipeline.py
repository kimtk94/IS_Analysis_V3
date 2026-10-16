"""Offline synthetic science gates; never use live genotype files in CI."""
import importlib.util
import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

FILE = Path(__file__).resolve().parents[1] / "scripts/is/harness_science_pipeline.py"
spec = importlib.util.spec_from_file_location("is_science_stage_under_test",FILE)
science = importlib.util.module_from_spec(spec)
spec.loader.exec_module(science)


def tsv(path,header,rows):
    with path.open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,delimiter="\t",fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)


class ScienceGateTests(unittest.TestCase):
    def test_variant_keys_not_evidence_of_causality(self):
        self.assertEqual(science.parse_variant("chr1:123:a:g"),("1",123,"A","G"))
        self.assertEqual(science.parse_variant("1:123:A:G"),science.parse_variant("chr1:123:A:G"))
        self.assertIsNone(science.parse_variant("rs12345"))
        self.assertIsNone(science.parse_variant("1:22:A:A"))
        self.assertTrue(science.palindromic(science.parse_variant("1:22:A:T")))
        self.assertFalse(science.palindromic(science.parse_variant("1:22:A:C")))
        self.assertIsNone(science.numeric_af("1.3"))
        self.assertIsNone(science.numeric_af("NA"))

    def test_independent_key_af_and_strand_qc_counts(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/"variants.tsv"
            header=["region_id","QTL_dataset","gene_id","mapped_variant_GRCh37",
                    "BBJ_exact_SNP","EAS_AIS_exact_SNP","QTL_alt_eaf_ac_div_an",
                    "BBJ_alt_eaf","EAS_AIS_alt_eaf"]
            values=[
             {"region_id":"R1","QTL_dataset":"GTEx_v8","gene_id":"ENSG0001",
              "mapped_variant_GRCh37":"chr1:15:A:C","BBJ_exact_SNP":"1:15:A:C",
              "EAS_AIS_exact_SNP":"1:15:A:C","QTL_alt_eaf_ac_div_an":"0.30",
              "BBJ_alt_eaf":"0.32","EAS_AIS_alt_eaf":"0.31"},
             {"region_id":"R1","QTL_dataset":"GTEx_v8","gene_id":"ENSG0001",
              "mapped_variant_GRCh37":"chr1:18:A:T","BBJ_exact_SNP":"1:18:A:T",
              "EAS_AIS_exact_SNP":"1:18:T:A","QTL_alt_eaf_ac_div_an":"0.30",
              "BBJ_alt_eaf":"0.32","EAS_AIS_alt_eaf":"0.90"},
             {"region_id":"R2","QTL_dataset":"QTL2","gene_id":"ENSG0002",
              "mapped_variant_GRCh37":"rs723","BBJ_exact_SNP":"1:19:A:G",
              "EAS_AIS_exact_SNP":"1:19:A:G","QTL_alt_eaf_ac_div_an":"",
              "BBJ_alt_eaf":"","EAS_AIS_alt_eaf":""},
            ]
            tsv(path,header,values)
            obj=science.independent_p2_qc(path,max_rows=10)
            c=obj["counters"]
            self.assertEqual(c["total_rows"],3)
            self.assertEqual(c["mapped_variant_parseable"],2)
            self.assertEqual(c["dual_allele_key_equal"],1)
            self.assertEqual(c["nonpal_dual_allele_key_equal"],1)
            self.assertEqual(c["qtl_eas_af_delta_gt_0_2_all_rows_diagnostic_only"],1)
            self.assertEqual(c.get("qtl_eas_af_delta_gt_0_2_source_dual_review",0),0)
            self.assertEqual(obj["fasta_status"],"BLOCKED_NO_APPROVED_GRCH37_FASTA")
            self.assertNotIn("coloc",obj["schema"])

    def test_boolean_exact_source_flags_are_not_independent_snp_ids(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/"boolean_flags.tsv"
            columns=["region_id","QTL_dataset","gene_id","mapped_variant_GRCh37",
                     "BBJ_exact_SNP","EAS_AIS_exact_SNP","QTL_alt_eaf_ac_div_an",
                     "BBJ_alt_eaf","EAS_AIS_alt_eaf"]
            tsv(path,columns,[
                {"region_id":"R1","QTL_dataset":"Q1","gene_id":"ENSG1",
                 "mapped_variant_GRCh37":"chr1:125:A:C","BBJ_exact_SNP":"True",
                 "EAS_AIS_exact_SNP":"True","QTL_alt_eaf_ac_div_an":"0.4",
                 "BBJ_alt_eaf":"0.4","EAS_AIS_alt_eaf":"0.41"},
                {"region_id":"R1","QTL_dataset":"Q1","gene_id":"ENSG1",
                 "mapped_variant_GRCh37":"chr1:130:C:G","BBJ_exact_SNP":"False",
                 "EAS_AIS_exact_SNP":"True","QTL_alt_eaf_ac_div_an":"0.4",
                 "BBJ_alt_eaf":"0.42","EAS_AIS_alt_eaf":"0.43"}
            ])
            report=science.independent_p2_qc(path)
            totals=report["counters"]
            self.assertEqual(totals["dual_source_reported_exact_snp"],1)
            self.assertEqual(totals.get("dual_allele_key_equal",0),0)
            self.assertEqual(totals["bbj_SOURCE_FALSE_FLAG"],1)
            self.assertEqual(totals["eas_SOURCE_TRUE_FLAG"],2)
            self.assertEqual(totals["mapped_variant_parseable"],2)
            p1={"rows":[{"gene_id":"ENSG1","region_id":"R1",
                          "QTL_dataset":"Q1","p1_status":"BLOCKED_NO_VERIFIED_COHORT_MATCHED_LD"}]}
            p3={"rows":[{"gene_id":"ENSG1","status":"BLOCKED"}]}
            p4=science.join_p4_candidates({"candidates":[{
                "gene_id":"ENSG1","gene_symbol":"GENE1",
                "nominal_min_p":0.01,"nominal_variant_rows":2
            }]},report,p1,p3)
            self.assertEqual(p4["candidates"][0]["evidence_tier"],"SOURCE_MATCH_FLAG_ONLY")

    def test_p1_ld_cohort_and_sha_must_match(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            ld=root/"verified_ld.tsv"
            ld.write_text("snp\tchr\nrs1\t1\n")
            digest=hashlib.sha256(ld.read_bytes()).hexdigest()
            template={c:"" for c in science.P1_COLUMNS}
            template.update({"region_id":"R1","QTL_dataset":"Q1","gene_id":"ENSG1",
                "qtl_source_cohort":"CohortA","source_ancestry":"EAS",
                "ld_source_cohort":"CohortA","ld_ancestry":"EAS",
                "ld_matrix_path":str(ld),"ld_sha256":digest,"verified_by":"PI",
                "verified_at":"2026-10-16","evidence_url":"https://example.org/evidence"})
            approve,counts=science.verified_p1_registry([template],root=root)
            self.assertEqual(len(approve),1)
            self.assertEqual(counts["PROVENANCE_REVIEWED_LD_HASH_MATCH"],1)
            bad={**template,"ld_source_cohort":"UnmatchedSource"}
            approve,counts=science.verified_p1_registry([bad],root=root)
            self.assertEqual(len(approve),0)
            self.assertEqual(counts["LD_COHORT_OR_ANCESTRY_MISMATCH"],1)
            bad={**template,"ld_sha256":"0"*64}
            approve,counts=science.verified_p1_registry([bad],root=root)
            self.assertEqual(len(approve),0)
            self.assertEqual(counts["LD_SHA256_MISMATCH"],1)

    def test_p3_no_coloc_without_independently_verified_fasta_and_inputs(self):
        p1={"rows":[{"region_id":"R1","QTL_dataset":"Q1","gene_id":"ENSG1"}]}
        p2={"groups":[{"region_id":"R1","QTL_dataset":"Q1","gene_id":"ENSG1",
                       "dual_allele_key_equal":20,"fasta_and_nonpal_allele_key_equal":20}]}
        gate=science.evaluate_p3_gate(p1,p2,{},[],p2_fasta_approved=False)
        self.assertEqual(gate["coloc_executed"],0)
        self.assertIn("BLOCKED_COHORT_MATCHED_LD",gate["rows"][0]["status"])
        self.assertIn("BLOCKED_INDEPENDENT_GRCH37_FASTA_ALLELE",gate["rows"][0]["status"])
        self.assertIn("BLOCKED_MISSING_REVIEWED_COLOC_INPUT",gate["rows"][0]["status"])

    def test_p4_independent_snp_evidence_is_not_a_causal_gene(self):
        cand={"candidates":[{"gene_id":"ENSG1","gene_symbol":"GENE1","nominal_min_p":0.001,"nominal_variant_rows":4},
                            {"gene_id":"ENSG2","gene_symbol":"GENE2","nominal_min_p":0.01,"nominal_variant_rows":2}]}
        p2={"counters":{"mapped_variant_parseable":8},
            "groups":[{"gene_id":"ENSG1","dual_allele_key_equal":5,"nonpal_dual_allele_key_equal":4}]}
        p1={"rows":[{"gene_id":"ENSG1","p1_status":"BLOCKED_NO_VERIFIED_COHORT_MATCHED_LD"}]}
        p3={"rows":[{"gene_id":"ENSG1","status":"BLOCKED","coloc_executed":False}]}
        out=science.join_p4_candidates(cand,p2,p1,p3)
        self.assertEqual(out["candidate_count"],2)
        self.assertEqual(out["tier_counts"]["HAS_INDEPENDENT_SNP_ID_OVERLAP"],1)
        self.assertEqual(out["candidates"][0]["p3_coloc_executed"],False)

    def test_missing_all_registries_generates_blocked_pipeline_not_fake_pass(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            output=root/"evidence"
            variant=root/"variants.tsv"
            ready=root/"ready.tsv"
            tsv(variant,["region_id","QTL_dataset","gene_id","mapped_variant_GRCh37",
                "BBJ_exact_SNP","EAS_AIS_exact_SNP","QTL_alt_eaf_ac_div_an",
                "BBJ_alt_eaf","EAS_AIS_alt_eaf"],[
                {"region_id":"R1","QTL_dataset":"Q1","gene_id":"ENSG1",
                 "mapped_variant_GRCh37":"1:1:A:C","BBJ_exact_SNP":"1:1:A:C",
                 "EAS_AIS_exact_SNP":"1:1:A:C","QTL_alt_eaf_ac_div_an":"0.1",
                 "BBJ_alt_eaf":"0.1","EAS_AIS_alt_eaf":"0.11"}])
            tsv(ready,["region_id","QTL_dataset","gene_id"],[
                {"region_id":"R1","QTL_dataset":"Q1","gene_id":"ENSG1"}])
            detail={"candidates":[{"gene_id":"ENSG1","gene_symbol":"GENE1",
                    "nominal_min_p":0.0001,"nominal_variant_rows":2}]}
            summary=science.advance(variant,ready,detail,output=output,
                registry=root/"missing1.tsv",reference=root/"missingref.json",
                p3_registry=root/"missing3.tsv",root=root)
            self.assertEqual(summary["p1"]["approved_matched_ld"],0)
            self.assertEqual(summary["p1"]["source_datasets"],1)
            self.assertEqual(summary["p3"]["coloc_executed"],0)
            self.assertEqual(summary["p4"]["candidate_count"],1)
            self.assertEqual(summary["p2"]["nonpalindromic_records"],1)
            self.assertEqual(summary["scientific_result"],"NO_NEW_COLOC_OR_CAUSAL_GENE_EVIDENCE")
            for name in ("p1_cohort_ld_review.json","p2_independent_variant_qc.json",
                         "p3_coloc_gates.json","p4_candidate_cross_evidence.json",
                         "science_pipeline_summary.json","p1_source_catalog.json"):
                self.assertTrue((output/"science"/name).is_file())


if __name__ == "__main__":
    unittest.main()
