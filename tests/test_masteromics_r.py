"""Actual scientific backend tests. Required (no skips) in the R CI job."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
import numpy as np
import pandas as pd
from scipy.stats import norm
from masteromics import science
from masteromics.compile import compile_project
from masteromics.engine import atomic_json,execute,sha
ROOT=Path(__file__).resolve().parents[1]
HAS_R=shutil.which('Rscript') is not None
REQUIRE_R=os.environ.get('MASTEROMICS_REQUIRE_R')=='1'

class ScientificBackendTest(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  if not HAS_R:
   if REQUIRE_R:raise RuntimeError('Required R backend unavailable')
   raise unittest.SkipTest('R integration runs in dedicated CI')
  subprocess.run(['Rscript','-e','for(p in c("TwoSampleMR","coloc","susieR","jsonlite","lme4","survival")) stopifnot(requireNamespace(p,quietly=TRUE))'],check=True)
 def bridge(self,*args):
  return subprocess.run(['Rscript',str(ROOT/'masteromics/analysis.R'),*map(str,args)],capture_output=True,text=True,timeout=180)
 def test_median_egger_and_eligibility(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);x=np.array([.05,.08,.11,.14,.17])
   d=pd.DataFrame({'key':[f's{i}' for i in range(5)],'beta_x':x,'beta_y':.5*x+np.array([-.0005,.0005,-.0003,.0002,.0001]),'se_x':.002,'se_y':.003})
   science.write(d,p/'h.tsv');science.write(d[['key']],p/'i.tsv')
   result=self.bridge('mr',p/'h.tsv',p/'i.tsv',p/'result.tsv');self.assertEqual(result.returncode,0,result.stdout+result.stderr)
   output=science.table(p/'result.tsv');self.assertEqual(set(output.status),{'SUCCESS'})
   self.assertTrue(np.allclose(output.beta,.5,atol=.03));self.assertTrue((output.se>0).all())
   science.write(d[['key']].iloc[:2],p/'i.tsv')
   result=self.bridge('mr',p/'h.tsv',p/'i.tsv',p/'small.tsv');self.assertEqual(result.returncode,0,result.stderr)
   self.assertEqual(set(science.table(p/'small.tsv').status),{'NOT_RUN_INSUFFICIENT_INSTRUMENTS'})
 def test_coloc_susie_full_dag_and_resume(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);n=120;idx=np.arange(n);r=.8**np.abs(idx[:,None]-idx[None,:]);keys=[f'GRCh38:1:{1000+i}:AG' for i in idx]
   registry={'datasets':{}}
   for did,effect in [('exposure',.3),('outcome',.2)]:
    b=effect*r[:,60]
    df=pd.DataFrame({'chr':'1','pos':1000+idx,'effect_allele':'A','other_allele':'G','beta':b,'se':.01,'pval':2*norm.sf(abs(b/.01)),'eaf':.25,'n':10000})
    science.write(df,p/f'{did}.tsv')
    registry['datasets'][did]={'id':did,'path':str(p/f'{did}.tsv'),'sha256':sha(p/f'{did}.tsv'),'columns':{c:c for c in science.REQ},'build':'GRCh38','ancestry':'EUR','effect_scale':'beta','type':'quant','sdY':1}
   pd.DataFrame(r,index=keys,columns=keys).to_csv(p/'ld.tsv',sep='\t')
   science.write(pd.DataFrame({'key':keys,'effect_allele':'A'}),p/'variants.tsv')
   ld={'matrix':str(p/'ld.tsv'),'variants':str(p/'variants.tsv'),'build':'GRCh38','ancestry':'EUR'}
   recipe={'project':'rtruth','output_root':str(p/'run'),'units':[{'id':'shared','gene':'TEST','exposure':'exposure','outcome':'outcome','build':'GRCh38','chr':'1','start':1000,'end':1119,'exposure_ld':ld,'outcome_ld':ld}]}
   atomic_json(p/'registry.json',registry);atomic_json(p/'recipe.json',recipe);compile_project(p/'recipe.json',p/'registry.json',p/'dag.json')
   self.assertEqual(execute(p/'dag.json',ROOT,jobs=2),0,(p/'run/run_summary.json').read_text())
   coloc=science.table(p/'run/shared_coloc/coloc.tsv');self.assertEqual(len(coloc),3);self.assertGreater(coloc['PP.H4.abf'].min(),.9)
   susie=science.table(p/'run/shared_susie/susie.tsv');self.assertEqual(set(susie.status),{'SUCCESS'});self.assertGreater(susie['PP.H4.abf'].max(),.9)
   self.assertAlmostEqual(science.table(p/'run/shared_mr/mr.tsv').iloc[0].beta,2/3,places=5)
   self.assertEqual(execute(p/'dag.json',ROOT,jobs=2),0)
   states=json.loads((p/'run/run_summary.json').read_text())['stages'];self.assertTrue(all(s['status']=='SKIPPED_EXISTING_VALID' for s in states.values()))
   # Independently reject a permuted LD axis in the R boundary.
   malformed=pd.DataFrame(r,index=keys,columns=keys).iloc[::-1];malformed.to_csv(p/'badld.tsv',sep='\t')
   result=self.bridge('susie',p/'run/shared_harmonize/regional.tsv',p/'compiled_inputs/unit_shared.json',p/'badld.tsv',p/'run/shared_ld/yld.tsv',p/'should_not_exist.tsv')
   self.assertNotEqual(result.returncode,0);self.assertIn('LD order mismatch',result.stderr)
 def test_longitudinal_and_incident_models(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);rng=np.random.default_rng(20261003);n=400
   scores=rng.normal(size=n);ages=rng.uniform(40,70,n);sex=rng.integers(0,2,n);pcs=rng.normal(size=n)
   intercept=rng.normal(0,7,n);slope=rng.normal(0,1.5,n);rows=[]
   for i in range(n):
    for time in range(5):
     rows.append({'id':f's{i}','time':time,'egfr':105-.3*ages[i]+intercept[i]+time*(-1-.7*scores[i]+slope[i])+rng.normal(0,1.5),'score':scores[i],'age':ages[i],'sex':sex[i],'PC1':pcs[i]})
   science.write(pd.DataFrame(rows),p/'long.tsv')
   result=self.bridge('cohort',p/'long.tsv',p/'mixed.tsv','PC1');self.assertEqual(result.returncode,0,result.stderr)
   mixed=science.table(p/'mixed.tsv');beta=mixed.set_index('term').loc['time:score','beta'];self.assertLess(abs(beta+.7),.3)
   eventtime=rng.exponential(10/np.exp(.6*scores));censor=10.;events=(eventtime<=censor).astype(int)
   d=pd.DataFrame({'id':[f's{i}' for i in range(n)],'followup':np.minimum(eventtime,censor),'event':events,'baseline_ckd':0,'score':scores,'age':ages,'sex':sex,'PC1':pcs})
   science.write(d,p/'incident.tsv');result=self.bridge('incident',p/'incident.tsv',p/'cox.tsv','PC1');self.assertEqual(result.returncode,0,result.stderr)
   cox=science.table(p/'cox.tsv').set_index('term');self.assertGreater(cox.loc['score','hr'],1);self.assertLess(abs(cox.loc['score','beta']-.6),.3)
if __name__=='__main__':unittest.main()
