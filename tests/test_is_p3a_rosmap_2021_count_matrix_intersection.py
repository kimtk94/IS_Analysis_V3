"""ROSMAP original 2021 matrix source/privacy gates: synthetic fixtures only."""
import csv
import gzip
import importlib
import io
import json
import unittest

M = importlib.import_module("scripts.is.audit_p3a_rosmap_2021_count_matrix_intersection")


def manifest():
    genes = [(f"ENSG{i:011}", f"GENE{i}") for i in range(20)]
    rows = []
    for q, indices in (("QTD000434", range(15)), ("QTD000131", range(10)),
                       ("QTD000171", range(9, 20))):
        rows += [{"region_id": "IS_XDATA_G0015", "QTL_dataset": q,
                  "gene_id": genes[i][0], "gene_symbol": genes[i][1]} for i in indices]
    return {"schema": "is-p3a-36-real-group-keys-from-legacy-393-v1", "groups": rows}


def metadata():
    out = io.StringIO()
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(("pseudo_id", "qtl_group", "rna_qc_passed", "genotype_qc_passed"))
    for i in range(581):
        writer.writerow((f"PRIVATE_SYNTHETIC_{i:03}", "brain_naive",
                         "True", "True" if i < 569 else "False"))
    return out.getvalue().encode()


def matrix(ncols=570, bad_row=False, duplicate=False):
    cols = [f"PRIVATE_SYNTHETIC_{i:03}" for i in range(ncols)]
    if duplicate:
        cols[1] = cols[0]
    out = io.StringIO()
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["gene_id"] + cols)
    for i in range(20):
        values = [f"ENSG{i:011}"] + ["3"] * ncols
        if bad_row and i == 7:
            values.pop()
        writer.writerow(values)
    return gzip.compress(out.getvalue().encode())


def r7(n=560):
    return ("study_id\tdataset_id\tstudy_label\tsample_group\tsample_size\n"
            f"QTS000025\tQTD000434\tROSMAP\tbrain_naive\t{n}\n").encode()


class OriginalSourceTests(unittest.TestCase):
    def test_expected_original_matrix_MD5(self):
        self.assertEqual(M.COUNT_NAME, "ROSMAP_merged_counts_pseudo.tsv.gz")
        self.assertEqual(M.COUNT_BYTES, 25484860)
        self.assertEqual(M.COUNT_MD5, "d215dc801f687e15d5b938356d80d183")

    def test_2021_intersection_and_36_group_status_no_identifiers(self):
        result = M.assemble(manifest(), metadata(), matrix(), r7(), {"MD5":"synthetic"})
        a = result["2021_matrix_aggregates"]
        self.assertEqual(a["2021_matrix_metadata_pseudo_ID_column_count"], 570)
        self.assertEqual(a["2021_matrix_and_RNA_genotype_QC_passed_intersection_count"], 569)
        self.assertEqual(a["2021_matrix_and_QC_failed_metadata_column_count"], 1)
        self.assertEqual(a["2021_QC_passed_metadata_samples_absent_from_matrix_count"], 0)
        self.assertEqual(a["original_15_ROSMAP_gene_ids_in_2021_count_matrix"], 15)
        self.assertEqual(result["n_gene_tissue_groups"], 36)
        self.assertEqual(result["n_formal_coloc_abf_executed"], 0)
        self.assertTrue(all(not x["coloc_eligible"] for x in result["groups"]))
        self.assertNotIn("PRIVATE_SYNTHETIC_", json.dumps(result))

    def test_even_560_matrix_matches_not_final_release7_proof(self):
        result = M.assemble(manifest(), metadata(), matrix(560), r7(), {})
        a = result["2021_matrix_aggregates"]
        self.assertTrue(a["2021_intersection_count_equals_Release7_560_DESCRIPTIVE"])
        self.assertEqual(a["2021_QC_passed_metadata_samples_absent_from_matrix_count"], 9)
        self.assertFalse(a["any_9_missing_donors_attributed_to_Release7_exclusion"])
        self.assertFalse(a["gene_specific_original_normalized_sdY_measured"])

    def test_metadata_duplicate_rejected(self):
        altered = metadata().replace(b"PRIVATE_SYNTHETIC_001", b"PRIVATE_SYNTHETIC_000")
        with self.assertRaises(ValueError):
            M.sample_sets(altered)

    def test_metadata_unknown_QC_rejected(self):
        changed = metadata().replace(b"\tTrue\tFalse\n", b"\tUNKNOWN\tFalse\n", 1)
        with self.assertRaises(ValueError):
            M.sample_sets(changed)

    def test_matrix_duplicate_columns_rejected(self):
        ids, passed = M.sample_sets(metadata())
        with self.assertRaises(ValueError):
            M.count_matrix_aggregates(matrix(duplicate=True), ids, passed,
                                      M.manifest_groups(manifest()))

    def test_matrix_no_matching_sample_id_rejected(self):
        ids, passed = M.sample_sets(metadata())
        with self.assertRaises(ValueError):
            M.count_matrix_aggregates(gzip.compress(b"gene_id\tunknown\nENSG00000000000\t4\n"),
                                      ids, passed, M.manifest_groups(manifest()))

    def test_matrix_malformed_row_rejected(self):
        ids, passed = M.sample_sets(metadata())
        with self.assertRaises(ValueError):
            M.count_matrix_aggregates(matrix(bad_row=True), ids, passed,
                                      M.manifest_groups(manifest()))

    def test_release7_wrong_N_fails(self):
        with self.assertRaises(ValueError):
            M.release7(r7(569))

    def test_release7_wrong_dataset_fails(self):
        with self.assertRaises(ValueError):
            M.release7(r7().replace(b"QTD000434", b"QTD000431"))

    def test_manifest_duplicate_gene_in_tissue_fails(self):
        m = manifest()
        m["groups"][0]["gene_id"] = m["groups"][1]["gene_id"]
        with self.assertRaises(ValueError):
            M.manifest_groups(m)

    def test_corrupt_gzip_rejected(self):
        ids, passed = M.sample_sets(metadata())
        with self.assertRaises((OSError, EOFError)):
            M.count_matrix_aggregates(b"not compressed sample data", ids, passed,
                                      M.manifest_groups(manifest()))


if __name__ == "__main__":
    unittest.main()
