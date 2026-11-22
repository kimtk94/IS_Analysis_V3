"""Synthetic negative gates for official three-study P3A molecular QTL metadata."""
import copy
import importlib
import json
import unittest
from pathlib import Path

M = importlib.import_module("scripts.is.audit_official_p3a_qtl_catalogue_provenance")
MANIFEST = Path("docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json")

COLS = ("# Column names of the nominal eQTL summary statistics files (*.all.tsv.gz)\n"
        "variant - Based on GRCh38 coordinates\n"
        "ref - GRCh38 reference allele\n"
        "**alt** - GRCh38 alternative allele (also the effect allele)\n"
        "**beta** - Regression coefficient\n"
        "**se** - Standard error\n"
        "**maf** - Minor allele frequency\n")

PHEAD = "study_id\tdataset_id\tstudy_label\tsample_group\ttissue_id\t" \
        "tissue_label\tcondition_label\tsample_size\tquant_method\tftp_path\t" \
        "ftp_cs_path\tftp_lbf_path\n"
MHEAD = "study_id\tdataset_id\tstudy_label\tsample_group\ttissue_id\t" \
        "tissue_label\tcondition_label\tsample_size\tquant_method\tpmid\tstudy_type\n"


def fixture():
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    p = [PHEAD]
    m = [MHEAD]
    for qtl, x in M.EXPECTED.items():
        tissue = {"QTD000131":"artery (aorta)", "QTD000171":"brain (cortex)",
                  "QTD000434":"brain (DLPFC)"}[qtl]
        ftp = ("ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
               + x["study_id"] + "/" + qtl + "/" + qtl + ".all.tsv.gz")
        common = [x["study_id"], qtl, x["study_label"], x["sample_group"],
                  x["tissue_id"],tissue, "naive",x["sample_size"],x["quant_method"]]
        p.append("\t".join(common+[ftp, "unused_cs","unused_lbf"])+"\n")
        m.append("\t".join(common+[x["pmid"],"bulk"])+"\n")
    return manifest, {"paths":"".join(p),"meta":"".join(m),"columns":COLS}


class QTLSourceTest(unittest.TestCase):
    def test_exact_official_datasets_and_real_36_groups(self):
        result=M.validate(*fixture())
        self.assertEqual(result["n_original_assay_groups"],36)
        self.assertEqual(result["n_tissue_datasets"],3)
        self.assertEqual(result["n_distinct_genes"],20)
        self.assertEqual(result["coloc_abf_eligible"],0)
        self.assertEqual(result["qtl_reference_build"],"GRCh38")
        self.assertFalse(result["liftover_complete"])
        self.assertTrue(all(not x["gene_effective_n_approved"] for x in result["group_rows"]))
        self.assertEqual([r["nominal_catalogue_sample_size"] for r in
                          result["QTL_datasets"].values()],[387,205,560])

    def test_study_sample_size_cannot_drift(self):
        manifest,sources=fixture()
        sources["meta"]=sources["meta"].replace("\t387\tge\t","\t388\tge\t")
        with self.assertRaises(ValueError):
            M.validate(manifest,sources)

    def test_wrong_population_source_rejected(self):
        manifest,sources=fixture()
        sources["meta"]=sources["meta"].replace("\tROSMAP\t","\tGTEx\t")
        with self.assertRaises(ValueError):
            M.validate(manifest,sources)

    def test_wrong_tissue_rejected(self):
        manifest,sources=fixture()
        sources["paths"]=sources["paths"].replace("\tbrain_cortex\t","\tbrain_other\t")
        with self.assertRaises(ValueError):
            M.validate(manifest,sources)

    def test_wrong_QTL_effect_build_rejected(self):
        manifest,sources=fixture()
        sources["columns"]=sources["columns"].replace("GRCh38","GRCh37")
        with self.assertRaises(ValueError):
            M.validate(manifest,sources)

    def test_wrong_nominal_source_endpoint_rejected(self):
        manifest,sources=fixture()
        sources["paths"]=sources["paths"].replace("QTD000434.all.tsv.gz",
                                                 "QTD000434.cc.tsv.gz")
        with self.assertRaises(ValueError):
            M.validate(manifest,sources)

    def test_missing_study_rejected(self):
        manifest,sources=fixture()
        sources["meta"]="\n".join(x for x in sources["meta"].split("\n")
                                   if "QTD000171" not in x)
        with self.assertRaises(ValueError):
            M.validate(manifest,sources)

    def test_extra_assay_cannot_be_approved(self):
        manifest,sources=fixture()
        manifest["groups"][0]["approved_for_new_coloc"]=True
        with self.assertRaises(ValueError):
            M.validate(manifest,sources)

    def test_wrong_QTL_study_membership_rejected(self):
        manifest,sources=fixture()
        manifest["groups"][0]["QTL_dataset"]="QTD000171"
        with self.assertRaises(ValueError):
            M.validate(manifest,sources)

    def test_git_blob_known_identity(self):
        raw=b"fixture\n"
        self.assertEqual(len(M.git_blob(raw)),40)
        self.assertNotEqual(M.git_blob(raw),M.git_blob(raw+b"changed"))

    def test_http_head_does_not_imply_data_download(self):
        manifest,sources=fixture()
        result=M.validate(manifest,sources)
        class FakeResponse:
            status=200
            headers={"Content-Length":"9999"}
            def __enter__(self):return self
            def __exit__(self,*args):return None
        calls=[]
        def fake_head(req,timeout=24):
            self.assertEqual(req.get_method(),"HEAD")
            calls.append(req.full_url)
            return FakeResponse()
        M.probe_head(result,fake_head)
        self.assertEqual(len(calls),6)
        self.assertEqual(sum(url.endswith(".tbi") for url in calls),3)
        self.assertFalse(result["any_source_fullcis_associations_downloaded"])
        self.assertEqual(result["coloc_abf_eligible"],0)
        self.assertTrue(all(not row["source_file_variant_rows_verified"] for row
                            in result["QTL_datasets"].values()))


if __name__ == "__main__":
    unittest.main()
