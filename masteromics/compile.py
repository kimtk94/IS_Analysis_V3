"""Compile project recipes into one validated DAG; identical code for every phenotype."""
import json
from pathlib import Path
from .engine import atomic_json,check_graph


def compile_project(recipe_path, registry_path, output):
    recipe=json.loads(Path(recipe_path).read_text()); registry=json.loads(Path(registry_path).read_text())
    mode=recipe.get('mode')
    if mode=='frozen_locus_first_migration':
        raise ValueError('Canonical ischemic-stroke config is a read-only locus-first migration contract; use `python -m masteromics is-adapter`, not generic pQTL run')
    if mode=='legacy_pqtl_reference':
        raise ValueError('Legacy ischemic-stroke pQTL recipe is reference-only and cannot be compiled as the canonical thesis pipeline')
    project=recipe['project']; variables=recipe.get('variables',{})
    stages=[]; used=set(); mr_outputs=[]; evidence_outputs=[]
    def stage(id,deps,argv,inputs,outputs):
        stages.append({'id':id,'depends':deps,'argv':argv,'inputs':inputs,'outputs':outputs})
    def contract(path,columns=None,keys=None):
        return {'path':path,'kind':'json','keys':keys} if keys else {'path':path,'kind':'tsv','columns':columns,'min_rows':1}
    def result(id,file):return '{results}/'+id+'/'+file
    py=['{python}','-m','masteromics']
    for u in recipe['units']:
        uid=u['id']
        if not uid.replace('_','').isalnum(): raise ValueError('Unsafe unit id')
        for did in [u['exposure'],u['outcome']]:
            if did in used:continue
            if not did.replace('_','').isalnum():raise ValueError('Unsafe dataset id')
            spec=registry['datasets'][did]
            if spec['id']!=did:raise ValueError('Registry identity mismatch')
            meta=Path(output).parent/'compiled_inputs'/f'{did}.json'
            atomic_json(meta,spec); metapath=str(meta.resolve())
            acq='acquire_'+did; norm='normalize_'+did
            raw=spec['path']; marker=result(acq,'integrity.json'); normalized=result(norm,'full.tsv')
            stage(acq,[],py+['acquire',metapath,'@out0'],[metapath]+([raw] if Path(raw).is_file() else []),[contract(marker,keys=['path','sha256'])])
            stage(norm,[acq],py+['normalize',metapath,'@out0'],[metapath,raw,marker],[contract(normalized,['key','beta','se','pval','build','ancestry'])])
            used.add(did)
        xs=registry['datasets'][u['exposure']];ys=registry['datasets'][u['outcome']]
        if xs['build']!=u['build'] or ys['build']!=u['build']:raise ValueError('Prepare validated same-build inputs before compilation')
        unit={**u,'exposure':xs,'outcome':ys,'min_regional_snps':u.get('min_regional_snps',100),
            'p1':u.get('p1',1e-4),'p2':u.get('p2',1e-4),'p12':u.get('p12',1e-5),
            'p12_grid':u.get('p12_grid',[1e-6,1e-5,1e-4]),'maxit':u.get('maxit',1000)}
        meta=Path(output).parent/'compiled_inputs'/f'unit_{uid}.json';atomic_json(meta,unit);meta=str(meta.resolve())
        exposure=result('normalize_'+u['exposure'],'full.tsv');outcome=result('normalize_'+u['outcome'],'full.tsv')
        instrument=uid+'_instruments';harm=uid+'_harmonize';mr=uid+'_mr';advanced=uid+'_sensitivity';coloc=uid+'_coloc';ld=uid+'_ld';susie=uid+'_susie'
        inst=result(instrument,'instruments.tsv'); matched=result(harm,'regional.tsv')
        ld_inputs=[u['exposure_ld']['matrix'],u['exposure_ld']['variants']]
        stage(instrument,['normalize_'+u['exposure']],py+['instruments',meta,exposure,'@out0','@out1'],[meta,exposure]+ld_inputs,
              [contract(inst,['key','f_statistic']),contract(result(instrument,'qc.json'),keys=['regional_rows','clumped_rows'])])
        hqc=result(harm,'qc.json')
        stage(harm,['normalize_'+u['exposure'],'normalize_'+u['outcome']],py+['harmonize',meta,exposure,outcome,'@out0','@out1'],[meta,exposure,outcome],
              [contract(matched,['key','beta_x','beta_y']),contract(hqc,keys=['matched','retained'])])
        mro=result(mr,'mr.tsv');mr_outputs.append(mro)
        stage(mr,[instrument,harm],py+['mr',matched,inst,'@out0'],[matched,inst],[{**contract(mro,['method','beta','se','pval']), 'finite_columns':['beta','se','pval'],'probability_columns':['pval']}])
        if recipe.get('advanced_mr',True):
            adv=result(advanced,'sensitivity.tsv');mr_outputs.append(adv)
            stage(advanced,[instrument,harm],['Rscript','{code}/masteromics/analysis.R','mr',matched,inst,'@out0'],[matched,inst],
                  [contract(adv,['method','status','intercept_pval'])])
        if recipe.get('regional_analysis',True):
            co=result(coloc,'coloc.tsv');evidence_outputs.append(co)
            stage(coloc,[harm],py+['regional_gate',meta,matched,hqc,'@out0'],[meta,matched,hqc],[{**contract(co,['p12','PP.H4.abf']), 'finite_columns':['PP.H4.abf'],'probability_columns':['PP.H4.abf']}])
            xld=result(ld,'xld.tsv');yld=result(ld,'yld.tsv')
            stage(ld,[harm],py+['ld',meta,matched,'@out0','@out1'],[meta,matched]+ld_inputs+[u['outcome_ld']['matrix'],u['outcome_ld']['variants']],
                  [contract(xld,[]),contract(yld,[])])
            so=result(susie,'susie.tsv');evidence_outputs.append(so)
            stage(susie,[coloc,ld],['Rscript','{code}/masteromics/analysis.R','susie',matched,meta,xld,yld,'@out0'],[meta,matched,xld,yld],
                  [contract(so,['status','PP.H4.abf'])])
        if recipe.get('annotation'):
            an=uid+'_annotation';ao=result(an,'annotation.tsv');evidence_outputs.append(ao)
            stage(an,[],py+['annotate',recipe['annotation'],u['gene'],'@out0'],[recipe['annotation']],[contract(ao,['gene_symbol','cell_type','expression','source'])])
    if recipe.get('cohort'):
        c=recipe['cohort'];co=result('cohort','longitudinal.tsv');evidence_outputs.append(co)
        cohortdeps=[]; panel=c.get('panel')
        if c.get('dosages'):
            scoreout=result('score','scores.tsv')
            stage('score',[],py+['score',c['dosages'],c['weights'],'@out0'],[c['dosages'],c['weights']],[contract(scoreout,['id','score'])])
            panel=result('cohort_panel','panel.tsv')
            stage('cohort_panel',['score'],py+['cohort_panel',c['phenotypes'],scoreout,'@out0'],[c['phenotypes'],scoreout],[contract(panel,['id','time','egfr','score'])])
            cohortdeps=['cohort_panel']
        if not panel:raise ValueError('Cohort needs a validated scored panel or dosages + weights + phenotypes')
        stage('cohort',cohortdeps,['Rscript','{code}/masteromics/analysis.R','cohort',panel,'@out0',','.join(c['pcs'])],[panel],
              [contract(co,['term','beta','se','model'])])
        if c.get('incident_panel'):
            coxout=result('incident','incident_ckd.tsv');evidence_outputs.append(coxout)
            stage('incident',[],['Rscript','{code}/masteromics/analysis.R','incident',c['incident_panel'],'@out0',','.join(c['pcs'])],[c['incident_panel']],[contract(coxout,['term','beta','se','pval','hr','ph_pval'])])
    # All scientifically requested stages must finish before the final evidence table.
    stage('report',[s['id'] for s in stages],py+['report','@out0']+mr_outputs,[*mr_outputs,*evidence_outputs],
          [contract('{results}/MASTER_MR_EVIDENCE.tsv',['method','beta','pval','qval','source_result'])])
    stage('evidence',[s['id'] for s in stages],py+['evidence','@out0']+mr_outputs+evidence_outputs,mr_outputs+evidence_outputs,[contract('{results}/MASTER_CANDIDATE_EVIDENCE.json',keys=['analyses','interpretation'])])
    check_graph(stages)
    atomic_json(output,{'project':project,'variables':variables,'output_root':recipe['output_root'],'stages':stages})
    return output
