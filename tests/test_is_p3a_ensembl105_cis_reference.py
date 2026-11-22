"""Synthetic acceptance tests: official Ensembl 105 reference ≠ approved true cis."""
import copy
import gzip
import importlib
import json
import unittest
from pathlib import Path

M=importlib.import_module("scripts.is.audit_p3a_ensembl105_cis_reference")
MANIFEST=Path("docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json")


def mock():
    m=json.loads(MANIFEST.read_text(encoding="utf8"))
    q={"schema":"is-p3a-36-original-eqtl-catalogue-indexed-chr10-audit-v1",
       "technical_status":"PASS_3_OFFICIAL_INDEXED_QTL_CONTEXT_REPORTS",
       "scientific_status":"HOLD","groups":[],
       "gene_full_cis_window_coverage_verified":False,
       "coloc_abf_eligible":0,
       "converted_GRCh38_queried_context":{
            **M.REGION,"is_lower_bound_of_full_cis_with_guaranteed_coverage":False}}
    pos={}
    for g in m["groups"]:
        pos[g["gene_id"]]={"gene_id":g["gene_id"],
             "chromosome":"10","phenotype_pos_GRCh38":103430000,
             "gene_name":g["gene_symbol"],
             "gene_strand":"+","original_phenotype_id":g["gene_id"]}
        q["groups"].append({"region_id":g["region_id"],
            "QTL_dataset":g["QTL_dataset"],"gene_id":g["gene_id"],
            "GRCh38_context_gene_associations_seen":2900,
            "GRCh38_context_valid_nonpal_unique":2300})
    return m,q,pos


class PositionTests(unittest.TestCase):
    def test_36_keys_20_annotation_gene_reference(self):
        m,q,p=mock()
        d=M.audit(m,q,p,{"publisher_MD5":M.EXPECTED_MD5})
        self.assertEqual(d["n_real_groups"],36)
        self.assertEqual(d["n_unique_ref_position_found"],20)
        self.assertEqual(d["n_reference_cis_entirely_in_prior_context"],20)
        self.assertEqual(d["n_group_coloc_abf_eligible"],0)
        self.assertTrue(all(not x["coloc_abf_eligible"] for x in d["groups"]))

    def test_right_side_missing_is_detected_in_all_tissues(self):
        m,q,p=mock()
        first=next(iter(p))
        p[first]["phenotype_pos_GRCh38"]=103600000
        d=M.audit(m,q,p,{})
        self.assertEqual(d["n_reference_cis_not_fully_covered"],1)
        self.assertEqual(d["genes"][first]["missing_right_bp"],139407)
        self.assertEqual(d["n_group_coloc_abf_eligible"],0)

    def test_left_side_missing(self):
        m,q,p=mock()
        first=next(iter(p))
        p[first]["phenotype_pos_GRCh38"]=103200000
        d=M.audit(m,q,p,{})
        self.assertEqual(d["genes"][first]["missing_left_bp"],203919)
        self.assertFalse(d["genes"][first]["reference_expected_window_fully_inside_common_context"])

    def test_gene_missing_held_not_imputed(self):
        m,q,p=mock()
        p.pop(next(iter(p)))
        d=M.audit(m,q,p,{})
        self.assertEqual(d["technical_status"],"HOLD_MISSING_GENE_REFERENCE")
        self.assertEqual(d["n_unique_ref_position_found"],19)
        self.assertEqual(d["n_group_coloc_abf_eligible"],0)

    def test_gene_cis_cannot_be_approved_from_reference_only(self):
        m,q,p=mock()
        d=M.audit(m,q,p,{})
        self.assertFalse(d["full_gene_specific_nominal_tested_cis_proven"])
        self.assertTrue(d["reference_105_positions_cannot_approve_dataset_specific_TSS"])
        self.assertTrue(all(x["source_study_exact_phenotype_metadata_matched"] is False
                            for x in d["groups"]))

    def test_wrong_scientific_status_blocks(self):
        m,q,p=mock()
        q["scientific_status"]="ELIGIBLE"
        with self.assertRaises(ValueError):M.audit(m,q,p,{})

    def test_faked_complete_cis_status_blocks(self):
        m,q,p=mock()
        q["gene_full_cis_window_coverage_verified"]=True
        with self.assertRaises(ValueError):M.audit(m,q,p,{})

    def test_qtl_unknown_source_group_blocks(self):
        m,q,p=mock()
        q["groups"][0]["QTL_dataset"]="QTDxxxxxx"
        with self.assertRaises(ValueError):M.audit(m,q,p,{})

    def test_wrong_manifest_gene_blocks(self):
        m,q,p=mock()
        m["groups"][0]["gene_id"]="ENSG00009999999"
        with self.assertRaises(ValueError):M.audit(m,q,p,{})

    def test_official_phenotype_gzip_gene_pos_extraction(self):
        header="phenotype_id\tgene_id\tchromosome\tphenotype_pos\tgene_name\tstrand\n"
        lines=(header+
            "\n".join(f"ENSG{j:011d}\tENSG{j:011d}\t10\t103430000\tG{j}\t+"
                      for j in range(351))+"\n")
        d=M.load_genes(gzip.compress(lines.encode()))
        self.assertEqual(len(d),351)
        self.assertEqual(d["ENSG00000000000"]["phenotype_pos_GRCh38"],103430000)

    def test_reject_duplicate_conflicting_phenotypes(self):
        header="phenotype_id\tgene_id\tchromosome\tphenotype_pos\n"
        v=header+"\n".join(f"ENSG{j:011d}\tENSG{j:011d}\t10\t103430000"
                               for j in range(351))
        v+="\nENSG00000000000\tENSG00000000000\t10\t104000000\n"
        with self.assertRaises(ValueError):M.load_genes(gzip.compress(v.encode()))

    def test_published_checksum_is_pinned(self):
        self.assertEqual(M.EXPECTED_MD5,"b999b21dcffeecbfe0bf886a5a91997b")
        self.assertTrue(M.URL.endswith("gene_counts_Ensembl_105_phenotype_metadata.tsv.gz?download=1"))


if __name__=="__main__":
    unittest.main()
