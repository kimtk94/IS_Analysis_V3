import json,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class ResourcesTest(unittest.TestCase):
 def test_shipped_catalog_and_manifest_provenance(self):
  import csv
  from masteromics.resources import DEFAULT,validate,plan
  catalog=validate(json.loads(DEFAULT.read_text()));indexed={r['dataset_id']:r for r in catalog['resources']}
  self.assertEqual(len(indexed),26)
  with (ROOT/'data/metadata/ckd_public_download_manifest.tsv').open() as f:
   for old in csv.DictReader(f,delimiter='\t'):
    if old['access']=='public':
     self.assertEqual(indexed[old['data_id']]['artifact_url'],old['url'])
     self.assertEqual(indexed[old['data_id']]['verification_scope'],'existing_manifest_not_probed')
  result=plan(catalog,list(indexed))
  self.assertEqual(result['downloads_executed'],0)
  self.assertTrue(all(i['status']=='BLOCKED' for i in result['items']))
 def test_access_cannot_be_overridden_and_pin_gates(self):
  from masteromics.resources import DEFAULT,plan
  catalog=json.loads(DEFAULT.read_text())
  b={'artifact_url':'https://example.org/fixture.tsv','file_name':'fixture.tsv','release':'fixture-v1',
     'sha256':'a'*64,'license_reviewed':True,'genome_build':'GRCh37','ancestry':'EUR','access_type':'open'}
  result=plan(catalog,['finngen_gwas','koges_cohort','ukb_ppp_pqtl'],{k:b for k in ['finngen_gwas','koges_cohort','ukb_ppp_pqtl']})
  self.assertIn('AUTHENTICATION_REQUIRED',result['items'][0]['blocks'])
  self.assertIn('DATA_ACCESS_APPROVAL_REQUIRED',result['items'][1]['blocks'])
  self.assertEqual(result['items'][2]['status'],'READY_FOR_ACQUISITION')
  self.assertEqual(result['items'][2]['scientific_status'],'NOT_REVIEWED')
  for field,gate in [('sha256','PIN_SHA256'),('release','PIN_RELEASE'),('genome_build','PIN_GENOME_BUILD'),('ancestry','PIN_ANCESTRY'),('license_reviewed','REVIEW_LICENSE_OR_DUA')]:
   bad=dict(b);bad[field]=None
   self.assertIn(gate,plan(catalog,['ukb_ppp_pqtl'],{'ukb_ppp_pqtl':bad})['items'][0]['blocks'])
 def test_catalog_and_selection_fail_closed(self):
  import copy
  from masteromics.resources import DEFAULT,validate,plan
  catalog=json.loads(DEFAULT.read_text())
  with self.assertRaisesRegex(ValueError,'Unknown'):plan(catalog,['unknown'])
  with self.assertRaises(ValueError):plan(catalog,[])
  with self.assertRaises(ValueError):plan(catalog,['gigastroke','gigastroke'])
  with self.assertRaises(ValueError):plan(catalog,['gigastroke'],{'megastroke':{}})
  bad=copy.deepcopy(catalog);bad['resources'].append(bad['resources'][0])
  with self.assertRaises(ValueError):validate(bad)
  bad=copy.deepcopy(catalog);bad['resources'][0]['source_url']='https://user:secret@example.org/data'
  with self.assertRaises(ValueError):validate(bad)
  with self.assertRaises(ValueError):plan(catalog,['gigastroke'],{'gigastroke':{'artifact_url':'https://example.org/file','file_name':'../escape'}})
  bad=copy.deepcopy(catalog);r=next(r for r in bad['resources'] if r['dataset_id']=='allen_human_brain');r['species']='Mus musculus'
  with self.assertRaisesRegex(ValueError,'human category'):validate(bad)
 def test_cli_works_without_scientific_environment(self):
  import subprocess
  result=subprocess.run([sys.executable,'-S','-m','masteromics','resources','validate'],cwd=ROOT,capture_output=True,text=True)
  self.assertEqual(result.returncode,0,result.stderr)
  result=subprocess.run([sys.executable,'-S','-m','masteromics','resources','plan','--ids','koges_cohort'],cwd=ROOT,capture_output=True,text=True)
  self.assertEqual(result.returncode,2,result.stderr)
  self.assertEqual(json.loads(result.stdout)['downloads_executed'],0)
if __name__=='__main__':unittest.main()
