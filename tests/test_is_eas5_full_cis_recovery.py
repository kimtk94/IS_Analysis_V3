"""Full cis tested-universe source completeness and AIS case fraction regression."""
import csv,gzip,importlib.util,sys,tempfile,unittest
from pathlib import Path
root=Path(__file__).resolve().parents[1]/"scripts/is"
def load(p,name):
 spec=importlib.util.spec_from_file_location(name,root/p)
 m=importlib.util.module_from_spec(spec)
 sys.modules[spec.name]=m
 spec.loader.exec_module(m)
 return m
merge=load("merge_is_eas5_full_cis_readiness.py","merge_cis_eas5")
ext=load("extend_is_eas5_public_qtl_complete_cis.py","ext_cis_eas5")
case=load("audit_is_gigastroke_eas_ais_case_fraction.py","case_EAS_AIS")
TESTED={
 "IS_XDATA_G0015":(37,41,43),
 "IS_XDATA_G0023":(9,7,6),
 "IS_XDATA_G0007":(3,4,4),
 "IS_XDATA_G0025":(51,58,58),
 "IS_XDATA_G0018":(25,25,22),
}
def fake():
 before=[];after=[]
 studies=("QTD000131","QTD000171","QTD000434")
 for reg,size in merge.REGIONS.items():
  for i,did in enumerate(studies):
   for j in range(size):
    tested=j<TESTED[reg][i]
    is_baseline_complete=reg=="IS_XDATA_G0015" and j<(3 if i<2 else 4)
    old_state="NO_EBI_PERMUTED_GENE_RECORD_CANNOT_ASSESS" if not tested else (
      "FULL_CIS_EQUIVALENCE_VS_EBI_NVARIANTS_COUNT_ONLY" if is_baseline_complete else
      "REGIONAL_SLICE_INCOMPLETE_CIS_NOT_SAFE_FOR_FULL_LOCUS_COLOC")
    new_state=merge.FULL if tested else merge.UNTESTED
    kw={"region_id":reg,"QTL_dataset":did,"gene_id":f"ENSG{reg[-4:]}{j:07d}",
      "gene_symbol":f"SYMBOL{j}"}
    before.append(dict(kw,
       archived_ABF_assays="0",
       nominal_regional_distinct_SNP_INDEL_variants=str(100 if is_baseline_complete else 50 if tested else 0),
       ebi_permuted_full_cis_n_variants="100" if tested else "",
       ebi_source_test_status=old_state))
    after.append(dict(kw,
       historical_ABF_tested="False",
       distinct_extended_cis_variant_ids="100" if tested else "0",
       full_cis_official_ebi_permuted_n_variants="100" if tested else "",
       official_permutation_lead_in_merged_nominal_QTL="True" if tested else "",
       cis_completeness_status=new_state))
 return before,after,{"case_fraction_status":"SECONDARY_PUBLISHED_RECONCILED_NOT_ORIGINAL_COHORT_LEVEL",
  "secondary_published_cases":"19032",
  "secondary_published_controls":"237242",
  "secondary_fraction_cases":str(19032/256274)}
class TestCISRecovery(unittest.TestCase):
 def test_720_after_recovers_383_and_preserves_327_unknown(self):
  rows,sub=merge.audit(*fake())
  self.assertEqual(len(rows),720)
  self.assertEqual(sum(int(x["new_full_cis_count_matched"]) for x in sub),393)
  self.assertEqual(sum(int(x["newly_restored_full_cis_count_matched"]) for x in sub),383)
  self.assertEqual(sum(int(x["no_permutation_record"]) for x in sub),327)
  self.assertTrue(all(x["coloc_recomputed"] is False for x in rows))
 def test_false_complete_status_fails_closed(self):
  b,a,c=fake()
  a[0]["distinct_extended_cis_variant_ids"]="99"
  with self.assertRaises(ValueError):merge.audit(b,a,c)
 def test_permutation_lead_variant_genomic_tuple(self):
  with tempfile.TemporaryDirectory() as t:
   p=Path(t)/"permuted.tsv.gz"
   with gzip.open(p,"wt") as f:
    f.write("molecular_trait_object_id\tvariant\n")
    f.write("ENSG123\tchr10_100300_A_G\n")
   leads=ext.official_permuted_lead_variants(p,{"ENSG123"})
   self.assertEqual(leads["ENSG123"],("10",100300,"A","G"))
 def test_EAS_AIS_case_fraction_consistent_with_EBI_yaml(self):
  with tempfile.TemporaryDirectory() as t:
   p=Path(t)/"GCST90104545-meta.yaml"
   p.write_text("genome_assembly: GRCh37\nsamples:\n  - sample_ancestry_category:\n      - East Asian\n    sample_size: 256274\n")
   row=case.check(p,"GCST90104545")
   self.assertAlmostEqual(row["secondary_fraction_cases"],19032/256274)
   self.assertFalse(row["coloc_casecontrol_input_from_original_cohort_table_attested"])
   p.write_text("genome_assembly: GRCh37\nsamples:\n  - sample_ancestry_category:\n      - East Asian\n    sample_size: 999\n")
   with self.assertRaises(ValueError):case.check(p,"GCST90104545")
if __name__=="__main__":unittest.main()
