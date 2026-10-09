"""G0022 provenance and matched-SNP subset audit (no unverified N transfer)."""
import csv,gzip,json,sys,tempfile,unittest
from pathlib import Path
import numpy as np
import yaml

SCRIPTS=Path(__file__).resolve().parents[1]/"scripts/is"
sys.path.insert(0,str(SCRIPTS))
from prepare_g0022_effective_n import neff,build as build_neff
from prepare_g0022_matched_as_ais_ld import run as build_match

class G0022DiagnosticsTest(unittest.TestCase):
 def test_neff_formula_rejects_wrong_sample_counts(self):
  val=neff(19032,237242,256274)
  self.assertAlmostEqual(val,70474.41011,places=3)
  self.assertAlmostEqual(neff(27413,237242,264655),98294.2313,places=3)
  with self.assertRaisesRegex(ValueError,"inconsistent"):
   neff(19032,237242,260000)
 def test_external_effective_n_source_matches_yaml(self):
  with tempfile.TemporaryDirectory() as td:
   base=Path(td);meta=base/"meta";meta.mkdir()
   for acc,total in [("GCST90104544",264655),("GCST90104545",256274)]:
    source={"samples":[{"sample_ancestry_category":["East Asian"],"sample_size":total}]}
    (meta/f"{acc}_buildGRCh37.tsv.gz-meta.yaml").write_text(yaml.safe_dump(source))
   r=build_neff(meta,base/"out")
   self.assertEqual(r["working_n"],{"AS":98294,"AIS":70474})
   self.assertIn("APPROXIMATE",r["status"])
 def test_same_ld_subset_and_order(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td)/"root";root.mkdir()
   gwas=Path(td)/"allele_matched.tsv.gz"
   template=[]
   for suffix in ("signal1_111629389","signal2_112930475"):
    f=root/("GCST90104544_AS_"+suffix);f.mkdir()
    vs=[]
    for i in range(130):
     pos=111000000+i if suffix.startswith("signal1") else 112900000+i
     vs.append({"variant_id":f"12:{pos}:A:G","gwas_z":str(i/10),
       "reference_alt_af":"0.30"})
     if i<120:
      template.append({"group_id":"IS_XDATA_G0022","dataset":"GCST90104545",
        "qc_status":"MATCH_ALT_EFFECT","variant_id":f"12:{pos}:A:G",
        "alt_effect_beta":str(0.1+0.001*i),"se":"0.1","alt_effect_eaf":"0.31"})
    with (f/"variants.tsv").open("w") as h:
     w=csv.DictWriter(h,delimiter="\t",fieldnames=list(vs[0]))
     w.writeheader();w.writerows(vs)
    np.savetxt(f/"signed_ld.tsv.gz",np.eye(len(vs)),delimiter="\t")
   with gzip.open(gwas,"wt") as h:
    w=csv.DictWriter(h,delimiter="\t",fieldnames=list(template[0]))
    w.writeheader();w.writerows(template)
   rows=build_match(root,gwas)
   self.assertEqual(len(rows),2)
   self.assertTrue(all(r["shared_study_variant_pairs"]==120 for r in rows))
   for r in rows:
    out=root/"matched_as_ais_ld_sensitivity"/r["window"]
    ld=np.loadtxt(out/"ld.tsv.gz",delimiter="\t")
    self.assertEqual(ld.shape,(120,120))
    self.assertEqual(len(np.loadtxt(out/"ais_z.tsv")),120)
    self.assertIn("NOT_INDEPENDENT",r["proof_level"])
if __name__=="__main__": unittest.main()
