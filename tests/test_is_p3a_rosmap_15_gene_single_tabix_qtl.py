"""Synthetic QTL source audit tests; no 4GB download or live EBI request in tests."""
from __future__ import annotations

import gzip
import hashlib
import importlib
import io
import json
import unittest
from pathlib import Path
from unittest.mock import patch

M=importlib.import_module("scripts.is.audit_p3a_rosmap_15_gene_single_tabix_qtl")


def manifest():
    return json.loads(Path(
        "docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json"
    ).read_text(encoding="utf-8"))


def candidates():
    return [g for g in manifest()["groups"] if g["QTL_dataset"]==M.DATASET]


def metadata():
    ret={}
    for i,g in enumerate(candidates()):
        start=102_000_000 + i*10_000
        ret[g["gene_id"]]={"chromosome":"10", "phenotype_pos":start,
                           "gene_interval_start":start, "gene_interval_end":start+999,
                           "strand":"+"}
    return ret


def region():
    return M.one_combined_region(metadata(),set(metadata()))


def source_row(g, pos=None, alt="G", an="1120", maf="0.22",
               beta="0.04", se="0.01", rsid="rs123", ref="C",
               chrom="10", pvalue="0.005", m_id=None):
    pos=pos or metadata()[g]["phenotype_pos"]+500
    raw={"molecular_trait_id":m_id or g,"chromosome":chrom,
         "position":str(pos),"ref":ref,"alt":alt,
         "variant":f"chr10_{pos}_{ref}_{alt}",
         "ma_samples":"240","maf":maf,"pvalue":pvalue,
         "beta":beta,"se":se,"type":"SNP","ac":"200","an":an,"r2":"0.99",
         "molecular_trait_object_id":g,"gene_id":g,
         "median_tpm":"1.1","rsid":rsid}
    return "\t".join(raw[x] for x in M.COLUMNS)+"\n"


def synthetic_rows(include_all=True):
    return [source_row(g["gene_id"]) for g in
            (candidates() if include_all else candidates()[:1])]


def metadata_file():
    lines=["\t".join(["phenotype_id","gene_id","phenotype_pos","chromosome",
                       "gene_start","gene_end","strand","gene_name"])]
    for i,g in enumerate(candidates()):
        pos=102_000_000 + i*10_000
        lines.append(f"{g['gene_id']}\t{g['gene_id']}\t{pos}\t10\t{pos}\t{pos+999}\t+\t{g['gene_symbol']}")
    return gzip.compress(("\n".join(lines)+"\n").encode())


class OriginalQTL15GeneOneRegionalTabixTests(unittest.TestCase):
    def test_manifest_contains_exactly_15_rosmap_targets(self):
        self.assertEqual(len(candidates()),15)
        self.assertEqual(len({r["gene_id"] for r in candidates()}),15)

    def test_2023_gene_metadata_publisher_checksum_guard(self):
        raw=metadata_file()
        h=hashlib.md5(raw,usedforsecurity=False).hexdigest()
        with patch.object(M,"Z_MD5",h):
            result=M.source_gene_metadata(raw,set(metadata()))
        self.assertEqual(len(result),15)

    def test_publisher_checksum_drift_fails_closed(self):
        with self.assertRaises(ValueError):
            M.source_gene_metadata(metadata_file(),set(metadata()))

    def test_original_gene_coordinate_lineage_not_approved_R7(self):
        r=region()
        self.assertEqual(r["chromosome"],"10")
        self.assertTrue(r["all_15_historical_phenotype_windows_contained"])
        self.assertFalse(r["R7_actual_qtlmap_run_cis_parameter_independently_verified"])

    def test_exactly_one_broad_region_not_15_queries(self):
        r=region()
        self.assertEqual(r["start"],101_000_000)
        self.assertEqual(r["end"],103_140_000)
        self.assertLess(r["span_bp"],M.MAX_REGION_BP)

    def test_incorrect_or_missing_gene_position_rejected(self):
        m=metadata();m.pop(next(iter(m)))
        with self.assertRaises(ValueError):
            M.one_combined_region(m,set(metadata()))

    def test_too_wide_union_fails_closed(self):
        m=metadata();m[next(iter(m))]["phenotype_pos"]=10_000_000
        with self.assertRaises(ValueError):
            M.one_combined_region(m,set(metadata()))

    def test_historical_phenotype_pos_not_within_gene_span_rejected(self):
        raw=metadata_file().replace(b"\n",b"\r\n")
        with patch.object(M,"Z_MD5",
                          hashlib.md5(raw,usedforsecurity=False).hexdigest()):
            self.assertEqual(len(M.source_gene_metadata(raw,set(metadata()))),15)

    def test_one_batch_recovers_all_15_with_expected_an(self):
        x=M.aggregate_tabix_stream(synthetic_rows(),manifest(),region(),metadata())
        self.assertEqual(x["n_ROSMAP_genes_with_original_any_SNP_QTL"],15)
        self.assertEqual(x["sum_gene_unique_full_variant_ids"],15)
        self.assertEqual(x["n_genes_all_AN_div2_560"],15)
        self.assertEqual(x["n_original_genes_approved_fullcis_for_formal_Coloc"],0)
        self.assertTrue(all(not g["coloc_abf_eligible"] for g in x["groups"]))

    def test_tabix_input_header_is_ignored_if_returned(self):
        rows=["\t".join(M.COLUMNS)+"\n"]+synthetic_rows()
        x=M.aggregate_tabix_stream(rows,manifest(),region(),metadata())
        self.assertEqual(x["n_ROSMAP_genes_with_original_any_SNP_QTL"],15)

    def test_missing_14_genes_are_scientific_hold_not_fabricated(self):
        x=M.aggregate_tabix_stream(synthetic_rows(include_all=False),
                                   manifest(),region(),metadata())
        self.assertEqual(x["n_ROSMAP_genes_with_missing_source_SNP_QTL"],14)
        self.assertEqual(x["n_genes_all_AN_div2_560"],1)

    def test_duplicate_rsid_same_effect_collapses(self):
        g=candidates()[0]["gene_id"]
        rows=synthetic_rows()+[source_row(g,rsid="rs456")]
        x=M.aggregate_tabix_stream(rows,manifest(),region(),metadata())
        self.assertEqual(x["sum_gene_unique_full_variant_ids"],15)
        self.assertEqual(x["sum_rsid_duplicates_removed_for_audit"],1)

    def test_duplicate_rsid_inconsistent_effect_blocks(self):
        g=candidates()[0]["gene_id"]
        rows=synthetic_rows()+[source_row(g,rsid="rs456",beta="0.09")]
        with self.assertRaises(ValueError):
            M.aggregate_tabix_stream(rows,manifest(),region(),metadata())

    def test_multiallelic_same_site_is_reported(self):
        g=candidates()[0]["gene_id"]
        rows=synthetic_rows()+[source_row(g,alt="A",rsid="rs789")]
        x=M.aggregate_tabix_stream(rows,manifest(),region(),metadata())
        self.assertEqual(x["sum_gene_unique_full_variant_ids"],16)
        self.assertEqual(x["groups"][0]["multi_alt_positions_present"],1)

    def test_an_variation_does_not_approve_n(self):
        g=candidates()[0]["gene_id"]
        rows=synthetic_rows()+[source_row(g,pos=metadata()[g]["phenotype_pos"]+600,
                                           an="1118")]
        x=M.aggregate_tabix_stream(rows,manifest(),region(),metadata())
        self.assertEqual(x["groups"][0]["AN_distinct_values"],2)
        self.assertFalse(x["groups"][0]["source_genewise_regression_N_certified"])

    def test_malformed_an_or_maf_fails(self):
        g=candidates()[0]["gene_id"]
        for field,invalid in [("an","1119"),("an","0"),("maf","1.1"),
                              ("maf","nan"),("se","0"),("se","nan"),
                              ("beta","nan"),("pvalue","2")]:
            with self.subTest(field=field,invalid=invalid):
                changed=source_row(g,**{field:invalid})
                with self.assertRaises(ValueError):
                    M.aggregate_tabix_stream([changed],manifest(),region(),metadata())

    def test_variant_id_alleles_and_position_must_match(self):
        g=candidates()[0]["gene_id"]
        row=source_row(g).replace("chr10_102000500_C_G",
                                   "chr10_102000500_T_G")
        with self.assertRaises(ValueError):
            M.aggregate_tabix_stream([row],manifest(),region(),metadata())

    def test_wrong_chromosome_and_out_of_range_variant_fails(self):
        g=candidates()[0]["gene_id"]
        for row in [source_row(g,chrom="12"),
                    source_row(g,pos=500_000_000),
                    source_row(g,pos=metadata()[g]["phenotype_pos"]+1_000_001)]:
            with self.assertRaises(ValueError):
                M.aggregate_tabix_stream([row],manifest(),region(),metadata())

    def test_19_column_schema_cannot_silently_drift(self):
        bad=source_row(candidates()[0]["gene_id"]).rstrip("\n").split("\t")[:-1]
        with self.assertRaises(ValueError):
            M.aggregate_tabix_stream(["\t".join(bad)+"\n"],
                                     manifest(),region(),metadata())

    def test_fails_when_max_row_budget_exceeded(self):
        g=candidates()[0]["gene_id"]
        with patch.object(M,"MAX_SOURCE_LINES",0):
            with self.assertRaises(ValueError):
                M.aggregate_tabix_stream([source_row(g)],manifest(),region(),metadata())

    def test_qtl_effects_and_sample_ids_not_exported(self):
        x=M.aggregate_tabix_stream(synthetic_rows(),manifest(),region(),metadata())
        s=json.dumps(x)
        self.assertNotIn("rs123",s)
        self.assertNotIn("chr10_102000500_C_G",s)
        self.assertNotIn('"beta": 0.04',s)
        self.assertFalse(x["raw_SNP_ids_effects_or_donor_ids_exported"])

    def test_legacy_eGene_counts_not_used_as_real_source_approval(self):
        x=M.aggregate_tabix_stream(synthetic_rows(),manifest(),region(),metadata())
        self.assertTrue(all(
            v["scientific_status"].startswith("HOLD") for v in x["groups"]))
        self.assertEqual(x["n_original_genes_approved_fullcis_for_formal_Coloc"],0)


if __name__=="__main__":
    unittest.main()
