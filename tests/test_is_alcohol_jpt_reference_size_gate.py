"""Synthetic fail-closed panel-size LD research scientific gate tests."""
import csv,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_alcohol_jpt_eas_reference_size_gates import audit
def tab(path,rows):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader();w.writerows(rows)
def setup(root):
    filtering=[];fixed=[];ctrl=[]
    for locus,n,original_s,jpt_s in [("ADH1B",637,.0133,.082),
                                        ("ALDH2",574,.0908,.194)]:
        p=root/locus; q=p/"JPT104_sensitivity";q.mkdir(parents=True)
        (q/"input_qc.json").write_text(json.dumps({
          "status":"JPT104_LD_SOURCE_RECONSTRUCTED_EAS_BASELINE_CROSSCHECK_PASS",
          "max_reconstructed_existing_EAS_signed_LD_absolute_difference":0,
          "JPT_n":104}))
        (p/"input_qc.json").write_text(json.dumps({"reference_EAS_n":504,"variants":n}))
        (p/"ALCOHOL_SUSIE_DIAGNOSTIC_AND_GATE.json").write_text(json.dumps({
          "median_sample_n":154570,"LD_mismatch_s":original_s}))
        for pop,snp,s in (("EAS504",n,original_s),("JPT104",n,jpt_s)):
            fixed.append({"locus":locus,"reference":pop,
                  "used_snp_count":snp,"source_snp_count":snp,
                  "source_genomic_variants_identical_across_references":"TRUE",
                  "LD_mismatch_s":s})
            for threshold in (.10,.05,.03):
                filtering.append({"locus":locus,"panel":pop,"max_abs_af_difference":threshold})
        for i in range(1,13):
            ctrl.append({"locus":locus,"replicate_index":i,"seed":20261010,
                         "diagnostic_s":.10+(i-6)*.011 if locus=="ADH1B"
                              else .22+(i-6)*.013,
                         "MR_or_causal_finemapping_performed":"FALSE"})
    tab(root/"ALCOHOL_ALDH2_ADH1B_EAS504_VS_JPT104_LD_MISMATCH.tsv",filtering)
    tab(root/"ALCOHOL_EAS504_JPT104_FIXED_IDENTICAL_SNP_LD_MISMATCH.tsv",fixed)
    tab(root/"ALCOHOL_EAS104_RANDOM_DOWNSAMPLING_LD_MISMATCH.tsv",ctrl)
    (root/"ALCOHOL_JPT104_EAS104_PANEL_SIZE_CONTROL_SUMMARY.json").write_text(
       json.dumps({"repeats":12}))
class ScienceGate(unittest.TestCase):
    def test_only_descriptive_ld_and_no_causal_iv(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);setup(root)
            r=audit(root)
            self.assertTrue(r["ALDH2_finemap_blocked"])
            self.assertFalse(r["susie_finite_reference_correction_applied"])
            self.assertFalse(r["quantitative_two_sample_MR_or_mediation_performed"])
            self.assertEqual(r["comparison_basis"],
                  "EXACT_IDENTICAL_SNP_LIST_NO_ANCESTRY_SPECIFIC_AF_REMOVAL")
    def test_variant_set_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);setup(root)
            p=root/"ALCOHOL_EAS504_JPT104_FIXED_IDENTICAL_SNP_LD_MISMATCH.tsv"
            with p.open() as f:rows=list(csv.DictReader(f,delimiter="\t"))
            rows[-1]["used_snp_count"]="571"
            tab(p,rows)
            with self.assertRaisesRegex(ValueError,"Different JPT/EAS SNP"):
                audit(root)
if __name__=="__main__":unittest.main()
