import json,tempfile,unittest,sys
from pathlib import Path
import numpy as np
import pandas as pd
from masteromics.engine import execute,atomic_json,sha,validate,check_graph
from masteromics import science
from masteromics.compile import compile_project
ROOT=Path(__file__).resolve().parents[1]
class CoreTest(unittest.TestCase):
 def test_checkpoint(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);inp=p/'input';inp.write_text('first');out=p/'result.tsv';cfg=p/'cfg.json';script=p/'worker.py'
   script.write_text("import sys\nfrom pathlib import Path\nPath(sys.argv[2]).write_text('value\\n'+Path(sys.argv[1]).read_text()+'\\n')\n")
   atomic_json(cfg,{'project':'fixture','output_root':str(p/'run'),'stages':[{'id':'one','argv':[sys.executable,str(script),str(inp),'@out0'],'inputs':[str(inp)],'outputs':[{'path':str(out),'kind':'tsv','columns':['value']}] }]})
   self.assertEqual(execute(cfg,ROOT),0);self.assertEqual(execute(cfg,ROOT),0)
   self.assertEqual(json.loads((p/'run/run_summary.json').read_text())['stages']['one']['status'],'SKIPPED_EXISTING_VALID')
   out.write_text('value\ncorrupted\n');self.assertEqual(execute(cfg,ROOT),0);self.assertIn('first',out.read_text())
   inp.write_text('second');self.assertEqual(execute(cfg,ROOT),0);self.assertIn('second',out.read_text())
   script.write_text("import sys\nfrom pathlib import Path\nPath(sys.argv[2]).write_text('value\\n')\n")
   self.assertEqual(execute(cfg,ROOT),1);self.assertIn('second',out.read_text())
 def test_header_graph(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'empty.tsv';p.write_text('key\n')
   with self.assertRaises(ValueError):validate(p,{'kind':'tsv','columns':['key']})
  with self.assertRaises(ValueError):check_graph([{'id':'a','depends':['b'],'outputs':[{}]}])
 def test_pipeline(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);registry={'datasets':{}}
   for did,b in [('x',.1),('y',.05)]:
    frame=pd.DataFrame({'chr':['1','1'],'pos':[100,200],'effect_allele':['A','C'],'other_allele':['G','T'],'beta':[b,b*2],'se':[.01,.01],'pval':[1e-10,1e-11],'eaf':[.2,.3],'n':[10000,10000]})
    frame.to_csv(p/f'{did}.tsv',sep='\t',index=False)
    registry['datasets'][did]={'id':did,'path':str(p/f'{did}.tsv'),'sha256':sha(p/f'{did}.tsv'),'columns':{c:c for c in science.REQ},'build':'GRCh38','ancestry':'EUR','effect_scale':'beta','type':'quant','sdY':1}
   keys=['GRCh38:1:100:AG','GRCh38:1:200:CT']
   pd.DataFrame(np.eye(2),index=keys,columns=keys).to_csv(p/'ld.tsv',sep='\t')
   pd.DataFrame({'key':keys,'effect_allele':['A','C']}).to_csv(p/'variants.tsv',sep='\t',index=False)
   ld={'matrix':str(p/'ld.tsv'),'variants':str(p/'variants.tsv'),'build':'GRCh38','ancestry':'EUR'}
   recipe={'project':'synthetic','output_root':str(p/'run'),'advanced_mr':False,'regional_analysis':False,'units':[{'id':'truth','gene':'TEST','exposure':'x','outcome':'y','build':'GRCh38','chr':'1','start':1,'end':1000,'exposure_ld':ld,'outcome_ld':ld}]}
   atomic_json(p/'registry.json',registry);atomic_json(p/'recipe.json',recipe);compile_project(p/'recipe.json',p/'registry.json',p/'dag.json')
   self.assertEqual(execute(p/'dag.json',ROOT,jobs=2),0)
   self.assertTrue(np.allclose(science.table(p/'run/MASTER_MR_EVIDENCE.tsv').beta,.5))
   self.assertEqual(execute(p/'dag.json',ROOT,jobs=2),0)
   science.normalize(registry['datasets']['x'],p/'norm.tsv');d=science.table(p/'norm.tsv')
   with self.assertRaises(ValueError):science.region(d,{'build':'GRCh37','chr':'1','start':1,'end':1000})
   with self.assertRaises(ValueError):science.load_ld({**ld,'ancestry':'EAS'},d)
 def test_swap_palindrome(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);d=pd.DataFrame({'chr':['1']*2,'pos':[10,20],'effect_allele':['A','A'],'other_allele':['G','T'],'beta':[.1,.1],'se':[.01]*2,'pval':[1e-9]*2,'eaf':[.2]*2,'n':[1000]*2})
   d.to_csv(p/'x.tsv',sep='\t',index=False)
   spec={'id':'x','path':str(p/'x.tsv'),'columns':{c:c for c in science.REQ},'build':'GRCh38','ancestry':'EUR','effect_scale':'beta'}
   science.normalize(spec,p/'xn.tsv');y=d.copy();y.effect_allele=d.other_allele;y.other_allele=d.effect_allele;y.beta=-d.beta*.5;y.to_csv(p/'y.tsv',sep='\t',index=False)
   science.normalize({**spec,'path':str(p/'y.tsv'),'id':'y'},p/'yn.tsv')
   science.harmonize(p/'xn.tsv',p/'yn.tsv',{'build':'GRCh38','chr':'1','start':1,'end':100},p/'h.tsv',p/'qc.json')
   h=science.table(p/'h.tsv');self.assertEqual(len(h),1);self.assertAlmostEqual(h.iloc[0].beta_y,.05)
class ScoreTest(unittest.TestCase):
 def test_allele_alignment_and_incomplete(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)
   w=pd.DataFrame({'key':['a','b'],'effect_allele':['A','C'],'other_allele':['G','T'],'build':['GRCh38']*2,'beta':[.1,.2]})
   d=pd.DataFrame({'id':['001','001'],'key':['a','b'],'effect_allele':['G','C'],'other_allele':['A','T'],'build':['GRCh38']*2,'dosage':[2.,1.]})
   science.write(w,p/'w.tsv');science.write(d,p/'d.tsv');science.score(p/'d.tsv',p/'w.tsv',p/'score.tsv')
   out=science.table(p/'score.tsv');self.assertEqual(out.iloc[0].id,'001');self.assertAlmostEqual(out.iloc[0].score,.2)
   science.write(d.iloc[:1],p/'d.tsv')
   with self.assertRaises(ValueError):science.score(p/'d.tsv',p/'w.tsv',p/'score.tsv')
if __name__=='__main__':unittest.main()
