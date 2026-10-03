"""CKD rebuild from local source archives, using reviewed legacy preparation adapters.

Reference LD is reused explicitly. Downloads, genotype LD rebuilding, IS and
cohort annotation are outside this migration run, never silently marked complete.
"""
import argparse
import csv
import fcntl
import gzip
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from .engine import atomic_json,sha,code_digest
from . import regression

ROOT=Path(__file__).resolve().parents[1]
SOURCES=[
 ('EUR','eGFRcrea','stanzick_egfr','outcome/EUR/ckdgen_stanzick2021_egfr_eur/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz'),
 ('EUR','CKD','wuttke_ckd','outcome/EUR/ckdgen_wuttke2019_ckd_eur/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz'),
 ('EUR','BUN','wuttke_bun','support/EUR/ckdgen_wuttke2019_bun_eur/BUN_overall_EA_YL_20171108_METAL1_nstud24.dbgap.txt.gz'),
 ('EUR','eGFRcys','gorski_egfrcys','support/EUR/ckdgen_gorski2017_egfrcys_eur/CKDGen_1000Genomes_DiscoveryMeta_eGFRcys_overall.csv.gz'),
 ('EUR','UACR','teumer_uacr','support/EUR/ckdgen_teumer2019_uacr_eur/formatted_20180517-UACR_overall-EA-nstud_18-SumMac_400.tbl.rsid.gz'),
 ('EAS','eGFRcrea','chen_egfr','eas/eGFR.gz'),('EAS','BUN','chen_bun','eas/BUN.gz')]

def table(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:return list(csv.DictReader(f,delimiter='\t'))

def compare(old,new,keys):
    a=table(old);b=table(new)
    def indexed(rows):
        out={}
        for r in rows:
            k=tuple(r[x] for x in keys)
            if k in out:raise ValueError('Duplicate comparison key: '+str(k))
            out[k]=r
        return out
    aa=indexed(a);bb=indexed(b);differences=[]
    for k in sorted(aa.keys()|bb.keys()):
        if k not in aa or k not in bb:
            differences.append({'key':k,'field':'row_presence'});continue
        if set(aa[k])!=set(bb[k]):differences.append({'key':k,'field':'schema'});continue
        for field,x in aa[k].items():
            y=bb[k][field]
            if x==y:continue
            try:
                fx,fy=float(x),float(y)
                same=(math.isnan(fx) and math.isnan(fy)) or (math.isfinite(fx) and math.isfinite(fy) and math.isclose(fx,fy,rel_tol=1e-6,abs_tol=1e-12))
            except (ValueError,TypeError):same=False
            if not same:differences.append({'key':k,'field':field,'legacy':x,'rebuilt':y})
    return {'status':'PASS' if not differences else 'DIFFERENCE','legacy_rows':len(a),'rebuilt_rows':len(b),'difference_count':len(differences),'examples':differences[:30],'legacy_sha256':sha(old),'rebuilt_sha256':sha(new)}

class Stages:
    def __init__(self,out):self.out=out;self.states={}
    def run(self,name,inputs,action):
        target=self.out/name;record=self.out/'checkpoints'/(name+'.json')
        h={str(p):sha(p) for p in inputs}
        scripts={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'scripts').glob('*.py'))}
        signature={'inputs':h,'code':code_digest(ROOT),'adapters':scripts,'python':sys.version}
        if record.exists() and target.exists():
            old=json.loads(record.read_text())
            if old['signature']==signature and old['outputs'] and all((target/p).is_file() and sha(target/p)==digest for p,digest in old['outputs'].items()):
                self.states[name]='SKIPPED_VALID';return target
        if target.exists():raise ValueError('Changed inputs/output: use a new --out instead of overwriting '+str(target))
        tmp=Path(tempfile.mkdtemp(prefix=name+'_',dir=self.out/'work'))
        try:
            action(tmp)
            files=[p for p in tmp.rglob('*') if p.is_file()]
            if not files:raise ValueError('No stage outputs')
            for p in files:
                if p.name.endswith(('.tsv','.tsv.gz')) and not table(p):raise ValueError('Header-only output: '+str(p))
            outputs={str(p.relative_to(tmp)):sha(p) for p in files}
            os.replace(tmp,target);atomic_json(record,{'signature':signature,'outputs':outputs})
            self.states[name]='SUCCESS';return target
        except BaseException:
            self.states[name]='FAILED';raise
        finally:
            if tmp.exists():shutil.rmtree(tmp)
    def command(self,name,args):
        with (self.out/'logs'/(name+'.log')).open('w') as log:
            r=subprocess.run([str(x) for x in args],stdout=log,stderr=subprocess.STDOUT)
        if r.returncode:raise RuntimeError(name+' failed; see '+str(self.out/'logs'/(name+'.log')))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('/srv/is-analysis'))
    p.add_argument('--out',type=Path,required=True);p.add_argument('--raw-root',type=Path)
    p.add_argument('--pqtl-root',type=Path);p.add_argument('--coordinates',type=Path);p.add_argument('--chain',type=Path)
    p.add_argument('--ld-root',type=Path);p.add_argument('--baseline',type=Path)
    p.add_argument('--rebuild-ld',action='store_true',help='Regenerate LD with reviewed bcftools/plink2 adapter in isolated work root; may download reference chromosomes')
    p.add_argument('--plan',action='store_true');p.add_argument('--max-ld-variants',type=int,default=8200)
    a=p.parse_args();a.raw_root=a.raw_root or a.root/'data/ckd/rawdata';a.pqtl_root=a.pqtl_root or a.root/'data/ckd/stage2_pqtl'
    a.coordinates=a.coordinates or a.root/'data/reference/gene_coordinates_hg38.tsv';a.chain=a.chain or a.root/'data/reference/hg38ToHg19.over.chain'
    a.ld_root=a.ld_root or a.root/'data/ckd/stage2c_ld/ld';a.baseline=a.baseline or a.root/'results/ckd'
    protected=[a.raw_root,a.pqtl_root,a.coordinates.parent,a.ld_root,a.baseline,ROOT]
    if any(a.out.resolve()==x.resolve() or x.resolve() in a.out.resolve().parents for x in protected):p.error('--out must be outside source, baseline and code roots')
    sources=[a.raw_root/'instruments/Sun2023.xlsx',*[a.raw_root/x[3] for x in SOURCES],a.coordinates,a.chain,ROOT/'data/metadata/ukb_ppp_download_manifest.tsv']
    missing=[str(x) for x in sources if not x.is_file()]
    print(json.dumps({'scope':'local raw statistics rebuild; reviewed legacy preparation; '+('regenerated LD' if a.rebuild_ld else 'cached reference LD')+'; CKD only','missing':missing,'raw_root':str(a.raw_root),'pqtl_root':str(a.pqtl_root),'ld_root':str(a.ld_root),'output':str(a.out)},indent=2),flush=True)
    if a.plan:return 1 if missing else 0
    if missing:p.error('Missing source files; pass explicit path overrides')
    for package in ['openpyxl','pyliftover']:__import__(package)
    subprocess.run(['Rscript','-e','for(p in c("jsonlite","coloc","susieR")) stopifnot(requireNamespace(p,quietly=TRUE))'],check=True)
    a.out.mkdir(parents=True,exist_ok=True)
    for name in ['work','logs','checkpoints']:(a.out/name).mkdir(exist_ok=True)
    lock=(a.out/'rebuild.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    cfg={k:str(v) for k,v in vars(a).items()};manifest=a.out/'REBUILD_CONFIG.json'
    if manifest.exists() and json.loads(manifest.read_text())!=cfg:p.error('Resume configuration differs; use a new output root')
    atomic_json(manifest,cfg);runner=Stages(a.out);status='INCOMPLETE';errors=[];comparisons={};rr=[];provenance={}
    try:
        def ready_action(t):
            inst=t/'UKBPPP_ST16_cis_independent.tsv.gz'
            runner.command('st16',[sys.executable,ROOT/'scripts/extract_ukbppp_st16_cis.py','--xlsx',sources[0],'--output',inst,'--summary',t/'instruments.summary.json'])
            for ancestry,phenotype,schema,relative in SOURCES:
                runner.command('extract_'+ancestry+'_'+phenotype,[sys.executable,ROOT/'scripts/extract_ckd_matched_outcomes.py','--schema',schema,'--input',a.raw_root/relative,'--instruments',inst,'--output',t/ancestry/(phenotype+'.tsv.gz'),'--summary',t/ancestry/(phenotype+'.summary.json'),'--phenotype',phenotype,'--ancestry',ancestry])
        ready=runner.run('ready',sources[:8],ready_action)
        stage1=runner.run('stage1',list(ready.rglob('*.tsv.gz')),lambda t:runner.command('stage1',[sys.executable,ROOT/'scripts/run_ckd_stage1_mr.py','--ready-root',ready,'--output-root',t]))
        stage2=runner.run('stage2',[stage1/'stage1_protein_summary.tsv',stage1/'instrument_qc.tsv.gz',sources[-1]],lambda t:runner.command('stage2',[sys.executable,ROOT/'scripts/run_ckd_stage2_candidates.py','--stage1-root',stage1,'--download-manifest',sources[-1],'--output-root',t]))
        genes=sorted({r['gene_symbol'] for r in table(stage2/'stage2_candidates.tsv')})
        # Archives are resolved by the same reviewed adapter, not guessed.
        import importlib.util
        spec=importlib.util.spec_from_file_location('ckd_prep',ROOT/'scripts/prepare_ckd_stage2b_coloc.py');adapter=importlib.util.module_from_spec(spec);spec.loader.exec_module(adapter)
        archives=[adapter.archive_for_gene(a.pqtl_root,g) for g in genes]
        regional=runner.run('stage2b_coloc',[stage2/'stage2_candidates.tsv',a.coordinates,a.chain,sources[1],*archives],lambda t:runner.command('regional',[sys.executable,ROOT/'scripts/prepare_ckd_stage2b_coloc.py','--candidates',stage2/'stage2_candidates.tsv','--pqtl-root',a.pqtl_root,'--gene-coordinates',a.coordinates,'--chain',a.chain,'--outcome-egfr',sources[1],'--output-root',t]))
        for ancestry,phenotype,_,_ in SOURCES:
            for sub,key in [('harmonized',['protein_id','rsid']),('mr',['protein_id'])]:
                rel=Path(sub)/ancestry/(phenotype+('.tsv.gz' if sub=='harmonized' else '.tsv'))
                comparisons[str(rel)]=compare(a.baseline/'stage1'/rel,stage1/rel,key)
        comparisons['candidates']=compare(a.baseline/'stage2/stage2_candidates.tsv',stage2/'stage2_candidates.tsv',['protein_id'])
        for gene in genes:comparisons['regional:'+gene]=compare(a.baseline/'stage2b_coloc/coloc_input'/(gene+'.tsv.gz'),regional/'coloc_input'/(gene+'.tsv.gz'),['snp'])
        # Prepare LD-aligned source files from freshly rebuilt full regions;
        # reference metadata/order and binary LD are deliberately reused.
        if a.rebuild_ld:
            for tool in ['bcftools','plink2','curl']:
                if not shutil.which(tool):raise ValueError('Missing LD regeneration tool: '+tool)
            rebuilt_ld=runner.run('ld_rebuild',[stage2/'stage2_candidates.tsv',*list((regional/'coloc_input').glob('*.tsv.gz'))],lambda t:runner.command('ld_rebuild',[sys.executable,ROOT/'scripts/prepare_ckd_stage2c_ld.py','--candidates',stage2/'stage2_candidates.tsv','--stage2b-input',regional/'coloc_input','--work-root',t/'reference','--output-root',t/'analysis']))
            a.ld_root=rebuilt_ld/'reference/ld';susie=rebuilt_ld/'analysis/susie_input'
        else:
            susie=a.out/'stage2c_susie/susie_input';susie.mkdir(parents=True,exist_ok=True)
            for gene in genes:shutil.copyfile(regional/'coloc_input'/(gene+'.tsv.gz'),susie/(gene+'.tsv.gz'))
        rr=[];policy=[]
        # Numerical MR uses rebuilt harmonization and historical estimates.
        for ancestry,phenotype,_,_ in SOURCES:
            mirror=Path(tempfile.mkdtemp(dir=a.out/'work'))
            try:
                (mirror/'harmonized').symlink_to(stage1/'harmonized',target_is_directory=True);(mirror/'mr').symlink_to(a.baseline/'stage1/mr',target_is_directory=True)
                rows,_=regression.replay_mr(mirror,a.out/'central_mr'/ancestry/phenotype,genes,1e-8,ancestry,phenotype);rr+=rows
            except Exception as e:errors.append({'stage':'central_mr','ancestry':ancestry,'phenotype':phenotype,'error':str(e)})
            finally:shutil.rmtree(mirror)
        for gene in genes:
            try:
                rows,delta,pp=regression.replay_coloc(gene,regional/'coloc_input'/(gene+'.tsv.gz'),a.baseline/'stage2b_coloc/coloc_results/STAGE2B_COLOC_DEFAULT.tsv',a.out/'central_coloc'/gene,1e-4);rr+=rows;policy+=delta;provenance['coloc:'+gene]=pp
                rows,pp=regression.replay_susie(gene,susie/(gene+'.tsv.gz'),a.baseline/'stage2c_susie/susie_results/STAGE2C_SUSIE_DEFAULT.tsv',a.ld_root,a.out/'central_susie'/gene,1e-4,a.max_ld_variants);rr+=rows;provenance['susie:'+gene]=pp
            except Exception as e:errors.append({'stage':'central_regional','gene':gene,'error':str(e)})
        import pandas as pd
        if rr:regression.science.write(pd.DataFrame(rr),a.out/'REBUILD_NUMERICAL_COMPARISON.tsv')
        if policy:regression.science.write(pd.DataFrame(policy),a.out/'POLICY_DELTA.tsv')
        status='INCOMPLETE' if errors or any(r['status']=='UNRESOLVED' for r in rr) else 'DIFFERENCE' if any(x['status']=='DIFFERENCE' for x in comparisons.values()) or any(r['status']=='DIFFERENCE' for r in rr) else 'PASS'
    except (Exception,SystemExit) as e:errors.append({'stage':'preparation','error':str(e)})
    atomic_json(a.out/'REBUILD_SUMMARY.json',{'status':status,'stages':runner.states,'errors':errors,'numerical_metrics':len(rr),'pass_metrics':sum(r['status']=='PASS' for r in rr),'difference_metrics':sum(r['status']=='DIFFERENCE' for r in rr),'unresolved_metrics':sum(r['status']=='UNRESOLVED' for r in rr),'provenance':provenance,'preparation_comparisons':comparisons,'scope':'CKD local source-statistics rebuild using reviewed legacy preparation adapters; source downloads, IS and cohort/tissue annotation excluded','LD_mode':'regenerated' if a.rebuild_ld else 'reused reference LD','R_library':os.environ.get('R_LIBS_USER','default'),'code_sha256':code_digest(ROOT)})
    print(json.dumps({'status':status,'errors':errors,'output':str(a.out)},indent=2));return 0 if status=='PASS' else 1
if __name__=='__main__':sys.exit(main())
