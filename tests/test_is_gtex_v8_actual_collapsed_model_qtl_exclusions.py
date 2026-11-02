import importlib.util
from pathlib import Path
import tempfile
import unittest

PATH=Path(__file__).resolve().parents[1]/"scripts/is/audit_gtex_v8_actual_collapsed_model_qtl_exclusions.py"
spec=importlib.util.spec_from_file_location("gtex_v8_collapsed",PATH)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def ln(chrom,start,end,strand,gene):
    return "\t".join((chrom,"gtex","exon",str(start),str(end),".",strand,".",
                      'gene_id "'+gene+'.11"; transcript_id "'+gene+'.11";'))

def run_rows(lines):
    with tempfile.NamedTemporaryFile(suffix=".gtf") as f:
        f.write(("\n".join(lines)+"\n").encode())
        f.flush()
        return m.inspect_gtf(f.name)

FGF5="ENSG00000138675"
NEURL="ENSG00000107954"
GOOD=[
    ln("chr4",80266599,80266700,"+",FGF5),
    ln("chr4",80267000,80267500,"+",FGF5),
    ln("chr10",103493979,103495000,"+",NEURL),
    ln("chr10",103497000,103500000,"+",NEURL),
]

class CollapsedModelTests(unittest.TestCase):
    def test_both_two_missing_qtd_snps_outside_gtEx_v8_window(self):
        result=run_rows(GOOD)
        self.assertEqual(result["source_model_QTL_boundary_gate"],
            "ACTUAL_GTEXV8_COLLAPSED_GENE_MODEL_SUPPORTS_TWO_CIS_EXCLUSIONS")
        self.assertTrue(result["two_target_SNPs_outside_actual_GTEx_v8_collapsed_cis_window"])
        self.assertEqual(result["genes"]["FGF5"]["QTD_only_position_bp_after_upper"],19)
        self.assertEqual(result["genes"]["NEURL"]["QTD_only_position_bp_before_lower"],229)
        self.assertFalse(result["coloc_performed"])

    def test_reverse_strand_uses_max_end_for_tss(self):
        lines=[GOOD[0],GOOD[1],ln("chr10",103490000,103493979,"-",NEURL)]
        out=run_rows(lines)
        self.assertEqual(out["genes"]["NEURL"]["GTEx_collapsed_TSS_GRCh38"],103493979)

    def test_gene_model_inside_cis_does_not_falsely_pass(self):
        rows=[ln("chr4",80300000,80300400,"+",FGF5),GOOD[2]]
        out=run_rows(rows)
        self.assertFalse(out["two_target_SNPs_outside_actual_GTEx_v8_collapsed_cis_window"])
        self.assertEqual(out["source_model_QTL_boundary_gate"],
            "HOLD_COLLAPSED_TSS_DISAGREES_WITH_QTL_EXCLUSION_HYPOTHESIS")

    def test_missing_gene_is_fail_closed(self):
        with self.assertRaises(ValueError):
            run_rows(GOOD[:2])

    def test_gene_strand_conflict_fails(self):
        with self.assertRaises(ValueError):
            run_rows([GOOD[0],ln("chr4",80267000,80267500,"-",FGF5),GOOD[2]])

    def test_public_source_file_size_capped(self):
        with tempfile.NamedTemporaryFile(suffix=".gtf") as f:
            f.write(b"X" * 101);f.flush()
            old=m.MAX_SOURCE_BYTES
            try:
                m.MAX_SOURCE_BYTES=100
                with self.assertRaises(ValueError):
                    m.inspect_gtf(f.name)
            finally:
                m.MAX_SOURCE_BYTES=old

    def test_no_statistical_promotion_from_boundary_match(self):
        out=run_rows(GOOD)
        self.assertFalse(out["exact_full_GTEx_original_file_MD5_recomputed"])
        self.assertFalse(out["QTL_GWAS_LD_N_approved"])
        self.assertTrue(out["not_evidence_of_independent_QTL_replication"])
        self.assertEqual(out["scientific_status"],"HOLD")

if __name__=="__main__":
    unittest.main()
