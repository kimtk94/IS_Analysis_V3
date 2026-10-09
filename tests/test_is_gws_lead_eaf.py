"""EAS source/reference ALT allele frequency comparison synthetic regression."""
import csv,gzip,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_gws_lead_eaf import run
class EAFTest(unittest.TestCase):
 def test_frequency_comparison_and_missingness(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp)/"results";root.mkdir()
   qc=root/"reference_frequency_qc";qc.mkdir()
   g=[]
   for i in range(7):
    group=f"G{i}"
    pos=100+i
    g.append(dict(group_id=group,lead_variant=f"2:{pos}:A:G",has_gws="1"))
    with (qc/f"{group}.afreq").open("w") as f:
     f.write("#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n")
     f.write(f"2\t2:{pos}:A:G\tA\tG\t0.20\t1008\n")
   with (root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv").open("w") as f:
    w=csv.DictWriter(f,fieldnames=list(g[0]),delimiter="\t");w.writeheader();w.writerows(g)
   with gzip.open(root/"IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz","wt") as f:
    cols=["group_id","pos","dataset","alt_effect_eaf","variant_id","qc_status","ancestry"]
    w=csv.DictWriter(f,fieldnames=cols,delimiter="\t");w.writeheader()
    for i in range(7):
     w.writerow(dict(group_id=f"G{i}",pos=100+i,dataset="EAS",
        alt_effect_eaf=".5" if i==0 else ".22",
        variant_id=f"2:{100+i}:A:G",qc_status="MATCH_ALT_EFFECT",ancestry="EAS"))
   report=run(root,Path(tmp)/"unused_new",Path(tmp)/"unused_old")
   self.assertEqual(report["lead_groups_reference_freq"],7)
   self.assertEqual(report["lead_gwas_frequency_pairs"],7)
   self.assertEqual(report["difference_above_0_10"],1)
   self.assertEqual(report["difference_at_most_0_10"],6)
if __name__=="__main__":unittest.main()
