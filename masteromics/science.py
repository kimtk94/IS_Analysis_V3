"""Explicit dataset adapters, signed LD QC, conservative harmonization and MR."""
from pathlib import Path
import hashlib
import json
import urllib.request
import os
import numpy as np
import pandas as pd
from scipy import stats
from .engine import sha

REQ=['chr','pos','effect_allele','other_allele','beta','se','pval','eaf','n']


def table(path): return pd.read_csv(path,sep='\t',dtype={'chr':str,'id':str,'key':str})

def write(df,path): df.to_csv(path,sep='\t',index=False,na_rep='NA')


def acquire(spec,path):
    from .acquisition import stage_file,verify
    source=Path(spec['path'])
    if source.is_file():
        # Reviewed local analysis inputs may live in read-only legacy directories.
        verify(source,spec.get('sha256'),spec.get('expected_bytes'),source.name.endswith(('.gz','.bgz')))
    else:
        stage_file(source,spec.get('sha256'),url=spec.get('url'),
                   identity={'dataset_id':spec.get('id'),'build':spec.get('build'),'ancestry':spec.get('ancestry')},
                   expected_bytes=spec.get('expected_bytes'))
    Path(path).write_text(json.dumps({'path':str(source.resolve()),'sha256':spec['sha256']}))


def normalize(spec,path):
    if spec.get('sha256') and sha(spec['path']) != spec['sha256']: raise ValueError('Dataset identity changed after acquisition')
    frame=pd.read_csv(spec['path'],sep=spec.get('separator','\t'),dtype=str)
    cols=spec['columns']
    if not set(REQ).issubset(cols): raise ValueError('Adapter must explicitly map all canonical columns')
    if not set(cols.values()).issubset(frame): raise ValueError('Source column missing')
    df=frame[list(cols.values())].rename(columns={v:k for k,v in cols.items()})
    df['chr']=df['chr'].str.replace(r'^chr','',regex=True)
    for c in ['pos','beta','se','pval','eaf','n']: df[c]=pd.to_numeric(df[c],errors='raise')
    for c in ['effect_allele','other_allele']: df[c]=df[c].str.upper()
    if spec.get('effect_scale') != 'beta': raise ValueError('Declare beta/logOR scale as beta; OR requires adapter conversion')
    if spec['build'] not in ['GRCh37','GRCh38']: raise ValueError('Unknown build')
    valid=np.isfinite(df[['pos','beta','se','pval','eaf','n']]).all(axis=1)
    valid &= (df.pos>0)&(df.pos==np.floor(df.pos))&(df.se>0)&df.pval.between(0,1,inclusive='right')&df.eaf.between(0,1,inclusive='neither')&(df.n>2)
    valid &= df.effect_allele.isin(list('ACGT'))&df.other_allele.isin(list('ACGT'))&(df.effect_allele!=df.other_allele)
    valid &= df.chr.isin([str(i) for i in range(1,23)]+['X','Y'])
    if not valid.all(): raise ValueError(f'Invalid summary-stat rows: {int((~valid).sum())}; resolve at adapter')
    df['pos']=df.pos.astype(int)
    df['build']=spec['build']; df['ancestry']=spec['ancestry']; df['dataset_id']=spec['id']
    df['key']=df['build']+':'+df.chr+':'+df.pos.astype(str)+':'+df.apply(lambda r:''.join(sorted([r.effect_allele,r.other_allele])),axis=1)
    if df.key.duplicated().any(): raise ValueError('Duplicate variant identities in one trait/assay; split registry dataset')
    if df.empty: raise ValueError('Empty source')
    df['source_sha256']=sha(spec['path']); df['source_file']=str(spec['path'])
    for field in ['trait_id','gene_symbol','protein_id','assay_id','variant_id','rsid']:
        if field not in df: df[field]=spec.get(field,'')
    write(df,path)


def region(df,u):
    if not isinstance(u['start'],int) or not isinstance(u['end'],int) or not 0<u['start']<=u['end']: raise ValueError('Define validated positive integer cis coordinates')
    if set(df.build)!={u['build']}: raise ValueError('Build mismatch: liftover and reference validation required first')
    return df[(df.chr==str(u['chr']))&df.pos.between(u['start'],u['end'])].copy()


def load_ld(spec,df):
    if set(df.ancestry)!={spec['ancestry']} or set(df.build)!={spec['build']}: raise ValueError('LD ancestry/build mismatch')
    matrix=pd.read_csv(spec['matrix'],sep='\t',index_col=0)
    alleles=table(spec['variants'])
    if not {'key','effect_allele'}.issubset(alleles): raise ValueError('Signed LD requires counted alleles')
    if alleles.key.duplicated().any() or matrix.index.duplicated().any(): raise ValueError('Duplicate LD variant')
    if list(matrix.index)!=list(matrix.columns): raise ValueError('LD axis order mismatch')
    if len(matrix)>spec.get('max_dense_variants',10000): raise ValueError('Dense LD memory limit')
    r=matrix.to_numpy(dtype=float)
    if not np.isfinite(r).all() or not np.allclose(r,r.T,atol=1e-6) or not np.allclose(np.diag(r),1,atol=1e-6) or np.max(np.abs(r))>1.000001: raise ValueError('Invalid signed LD correlation matrix')
    if np.linalg.eigvalsh(r).min() < -1e-6: raise ValueError('LD matrix not positive semidefinite')
    if not set(df.key).issubset(matrix.index) or not set(df.key).issubset(alleles.key): raise ValueError('LD coverage incomplete; no silent variant dropping')
    a=alleles.set_index('key').loc[df.key,'effect_allele'].to_numpy()
    if not np.all((a==df.effect_allele.to_numpy())|(a==df.other_allele.to_numpy())): raise ValueError('LD counted allele mismatch')
    sign=np.where(a==df.effect_allele.to_numpy(),1,-1)
    r=matrix.loc[df.key,df.key].to_numpy()*np.outer(sign,sign)
    return r


def instruments(exposure,u,ld,path,qc):
    df=region(table(exposure),u)
    if set(df.ancestry)!={ld['ancestry']}: raise ValueError('Exposure LD ancestry mismatch')
    nraw=len(df)
    df=df[(df.pval<=u.get('p_threshold',5e-8))&((df.beta/df.se)**2>=u.get('f_threshold',10))]
    if df.empty: raise ValueError('No strong cis instruments; record biological-empty separately before scaling')
    df=df.sort_values(['pval','key']).reset_index(drop=True); r=load_ld(ld,df)
    selected=[]
    for i in range(len(df)):
        if all(r[i,j]**2<=u.get('clump_r2',0.001) for j in selected): selected.append(i)
    kept=df.iloc[selected].copy(); kept['f_statistic']=(kept.beta/kept.se)**2; write(kept,path)
    Path(qc).write_text(json.dumps({'regional_rows':nraw,'strong_rows':len(df),'clumped_rows':len(kept),'ld_panel':ld}))


def harmonize(exposure,outcome,u,path,qc):
    x=region(table(exposure),u); y=region(table(outcome),u)
    if x.empty or y.empty: raise ValueError('Empty regional inputs')
    merged=x.merge(y,on='key',suffixes=('_x','_y'),validate='one_to_one')
    # Conservative palindromic policy: strand cannot be inferred reliably across populations.
    pair=merged.effect_allele_x+merged.other_allele_x
    pal=pair.isin(['AT','TA','CG','GC'])
    direct=(merged.effect_allele_x==merged.effect_allele_y)&(merged.other_allele_x==merged.other_allele_y)
    swap=(merged.effect_allele_x==merged.other_allele_y)&(merged.other_allele_x==merged.effect_allele_y)
    ok=~pal&(direct|swap)
    df=merged.loc[ok].copy()
    df['beta_y']=df.beta_y*np.where(swap.loc[ok],-1,1)
    df['eaf_y']=np.where(swap.loc[ok],1-df.eaf_y,df.eaf_y)
    if df.empty: raise ValueError('No nonambiguous matched variants')
    # Output columns sufficient for MR and full-region coloc; instruments never replace region.
    write(df,path)
    Path(qc).write_text(json.dumps({'exposure_region':len(x),'outcome_region':len(y),'matched':len(merged),'retained':len(df),'palindromic_drop':int(pal.sum()),'policy':'drop_all_palindromic','minimum_overlap_fraction':u.get('min_overlap_fraction',0.5)}))


def mr(harmonized,instrument_path,path):
    df=table(harmonized); selected=set(table(instrument_path).key); df=df[df.key.isin(selected)]
    if df.empty: raise ValueError('No harmonized instruments')
    x=df.beta_x.to_numpy(); y=df.beta_y.to_numpy(); sy=df.se_y.to_numpy(); k=len(df)
    if np.any(x==0): raise ValueError('Zero exposure effect')
    w=1/sy**2; b=np.sum(w*x*y)/np.sum(w*x*x); se=np.sqrt(1/np.sum(w*x*x))
    q=float(np.sum(w*(y-b*x)**2)); rows=[]
    def row(method,b,se,p,**extras): return dict(method=method,nsnp=k,beta=float(b),se=float(se),pval=float(p),**extras)
    rows.append(row('wald' if k==1 else 'ivw_fixed',b,se,2*stats.norm.sf(abs(b/se)),q=q,q_pval=float(stats.chi2.sf(q,k-1)) if k>1 else np.nan))
    if k>1:
        sem=se*np.sqrt(max(1,q/(k-1))); rows.append(row('ivw_multiplicative',b,sem,2*stats.norm.sf(abs(b/sem))))
    # Official package methods (median/Egger) are handled by R bridge, not approximate reimplementations.
    write(pd.DataFrame(rows),path)


def report(files,path):
    rows=[]
    for name in files:
        df=table(name); df['source_result']=name; rows.append(df)
    df=pd.concat(rows,ignore_index=True)
    if 'pval' in df:
        # BH within method across configured units/outcomes, recorded in report.
        df['qval']=np.nan
        for _,idx in df.groupby('method').groups.items():
            idx=[i for i in idx if pd.notna(df.loc[i,'pval'])]
            vals=df.loc[idx,'pval'].to_numpy(); order=np.argsort(vals); n=len(vals)
            if not n: continue
            adj=np.minimum.accumulate((vals[order]*n/np.arange(1,n+1))[::-1])[::-1]
            mapped=np.empty(n); mapped[order]=np.minimum(adj,1); df.loc[idx,'qval']=mapped
    write(df,path)


def regional_ld(harmonized,xspec,yspec,xout,yout):
    df=table(harmonized)
    for suffix,spec,out in [('x',xspec,xout),('y',yspec,yout)]:
        data=pd.DataFrame({'key':df.key,'build':df['build_'+suffix], 'ancestry':df['ancestry_'+suffix],
            'effect_allele':df.effect_allele_x,'other_allele':df.other_allele_x})
        r=load_ld(spec,data)
        pd.DataFrame(r,index=df.key,columns=df.key).to_csv(out,sep='\t',index=True)


def annotate(source,gene,path):
    df=table(source)
    if not {'gene_symbol','cell_type','expression','source'}.issubset(df): raise ValueError('Annotation schema')
    df=df[df.gene_symbol==gene].copy()
    df['expression']=pd.to_numeric(df.expression,errors='raise')
    if df.empty or not np.isfinite(df.expression).all() or (df.expression<0).any(): raise ValueError('Missing/invalid annotation')
    write(df.sort_values('expression',ascending=False),path)


def evidence(files,path):
    records=[]
    for filename in files:
        df=table(filename)
        unit=Path(filename).parent.name
        record={'unit_stage':unit,'source':filename,'sha256':sha(filename),'rows':len(df)}
        records.append({**record,'results':json.loads(df.to_json(orient='records'))})
    Path(path).write_text(json.dumps({'analyses':records,'interpretation':'Posterior evidence is separate from MR significance; no automatic causal proof'},indent=2,allow_nan=False))


def score(dosages,weights,path):
    d=table(dosages);w=table(weights)
    common={'key','effect_allele','other_allele','build'}
    if not common.union({'id','dosage'}).issubset(d) or not common.union({'beta'}).issubset(w):raise ValueError('Score input schema missing')
    if w.key.duplicated().any() or d[['id','key']].duplicated().any():raise ValueError('Duplicate score variant/sample')
    if set(w.build)!=set(d.build) or len(set(w.build))!=1:raise ValueError('Score genome build mismatch')
    d=d[d.key.isin(w.key)].merge(w,on='key',suffixes=('_d','_w'),validate='many_to_one')
    if d.empty or not (d.groupby('id').size()==len(w)).all():raise ValueError('Every participant must have every weighted variant; define imputation upstream')
    direct=(d.effect_allele_d==d.effect_allele_w)&(d.other_allele_d==d.other_allele_w)
    swap=(d.effect_allele_d==d.other_allele_w)&(d.other_allele_d==d.effect_allele_w)
    if not (direct|swap).all():raise ValueError('Score allele mismatch; strand harmonization required upstream')
    d['dosage']=pd.to_numeric(d.dosage,errors='raise');d['beta']=pd.to_numeric(d.beta,errors='raise')
    if not np.isfinite(d[['dosage','beta']]).all().all() or not d.dosage.between(0,2).all():raise ValueError('Invalid dosage/weights')
    d['weighted']=np.where(swap,2-d.dosage,d.dosage)*d.beta
    result=d.groupby('id',as_index=False).weighted.sum().rename(columns={'weighted':'score'});write(result,path)


def cohort_panel(phenotypes,scores,path):
    d=table(phenotypes);s=table(scores)
    if not {'id','time','egfr','age','sex'}.issubset(d) or not {'id','score'}.issubset(s):raise ValueError('Phenotype/score schema')
    if 'score' in d:raise ValueError('Input phenotype already contains score; avoid accidental overrides')
    if s.id.duplicated().any():raise ValueError('Duplicate score participants')
    if not set(d.id).issubset(s.id):raise ValueError('Missing scores; define sample exclusion upstream')
    write(d.merge(s,on='id',validate='many_to_one'),path)
