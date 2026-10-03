"""Replay frozen legacy inputs through central numerical modules, without refiltering.

This establishes method parity, not end-to-end raw-data equivalence. Optional
policy-delta coloc repeats after the central conservative palindromic filter.
"""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import re
import subprocess
import sys
import platform
import numpy as np
import pandas as pd
from . import science
from .engine import atomic_json,sha,code_digest

CORE=Path(__file__).resolve().parent
GENES=['ACP1','CPVL','F12','GSTA1','GSTA3','HLA-E','INHBC','SDCCAG8','UMOD']


def close(old,new,atol,rtol=1e-6):
    return bool(np.isfinite(float(old)) and np.isfinite(float(new)) and np.isclose(float(old),float(new),atol=atol,rtol=rtol))


def metric(stage,gene,field,old,new,atol,rtol=1e-6,comparison='SAME_INPUT'):
    valid=close(old,new,atol,rtol)
    return {'stage':stage,'gene_symbol':gene,'metric':field,'legacy':old,'central':new,
        'absolute_delta':abs(float(new)-float(old)),'status':'PASS' if valid else 'DIFFERENCE',
        'atol':atol,'rtol':rtol,'comparison':comparison}


def bridge(mode,data,meta,out,ld=None):
    args=['Rscript',str(CORE/'analysis.R'),mode,str(data),str(meta)]
    args+=([str(ld),str(ld),str(out)] if mode=='susie' else [str(out)])
    result=subprocess.run(args,capture_output=True,text=True,timeout=3600)
    (out.parent/(out.name+'.log')).write_text(result.stdout+'\n'+result.stderr)
    if result.returncode:raise RuntimeError(f'{mode} failed; see {out.name}.log: '+result.stderr[-1000:])


def legacy_rows(path):
    d=science.table(path)
    req=['snp','pos37','beta_pqtl','se_pqtl','maf_pqtl','n_pqtl','beta_outcome_aligned','se_outcome','maf_outcome','n_outcome']
    if not set(req).issubset(d):raise ValueError('Legacy regional schema missing columns')
    valid=np.isfinite(d[['beta_pqtl','se_pqtl','beta_outcome_aligned','se_outcome','maf_pqtl','maf_outcome']]).all(axis=1)
    valid &= (d.se_pqtl>0)&(d.se_outcome>0)&(d.maf_pqtl>0)&(d.maf_pqtl<=.5)&(d.maf_outcome>0)&(d.maf_outcome<=.5)
    return d[valid].drop_duplicates('snp',keep='first').copy()


def canonical_replay(d):
    # Already-aligned matched regional SNPs, same order as frozen legacy input.
    x=pd.DataFrame({'key':d.snp.astype(str),'pos_x':d.pos37.astype(int),
        'beta_x':d.beta_pqtl,'se_x':d.se_pqtl,'eaf_x':d.maf_pqtl,
        'n_x':int(round(d.n_pqtl.median())),
        'beta_y':d.beta_outcome_aligned,'se_y':d.se_outcome,'eaf_y':d.maf_outcome,
        'n_y':int(round(d.n_outcome.median()))})
    if x.empty or x.key.duplicated().any():raise ValueError('Empty/duplicate replay inputs')
    return x


def replay_meta():
    return {'exposure':{'type':'quant','sdY':1},'outcome':{'type':'quant','sdY_estimation':'coloc'},
        'p1':1e-4,'p2':1e-4,'p12':1e-5,'p12_grid':[1e-6,1e-5,1e-4],
        'min_regional_snps':50,'maxit':1000}


def replay_mr(stage1,out,genes,atol,ancestry='EUR',phenotype='eGFRcrea'):
    inp=stage1/'harmonized'/ancestry/(phenotype+'.tsv.gz')
    baseline=stage1/'mr'/ancestry/(phenotype+'.tsv')
    harmonized=science.table(inp);base=science.table(baseline)
    req={'protein_id','gene_symbol','rsid','beta_exposure','se_exposure','beta_outcome','se_outcome'}
    if not req.issubset(harmonized):raise ValueError('Stage1 harmonization schema mismatch')
    missing=set(genes)-set(base.gene_symbol)
    if missing:raise ValueError('Requested genes absent from MR baseline: '+','.join(sorted(missing)))
    rows=[]
    for _,old in base[base.gene_symbol.isin(genes)].iterrows():
        d=harmonized[harmonized.protein_id==old.protein_id].copy()
        uid=re.sub('[^A-Za-z0-9_]','_',str(old.gene_symbol)+'_'+str(old.protein_id))
        target=out/uid;target.mkdir(parents=True,exist_ok=True)
        data=pd.DataFrame({'key':d.rsid,'beta_x':d.beta_exposure,'se_x':d.se_exposure,'beta_y':d.beta_outcome,'se_y':d.se_outcome})
        if data.empty or data.key.duplicated().any():raise ValueError('Missing/duplicate frozen MR instruments: '+uid)
        science.write(data,target/'harmonized.tsv');science.write(data[['key']],target/'all_instruments.tsv')
        science.mr(target/'harmonized.tsv',target/'all_instruments.tsv',target/'ivw.tsv')
        fixed=science.table(target/'ivw.tsv').iloc[0]
        for oldcol,newcol in [('ivw_beta_sensitivity','beta'),('ivw_se_sensitivity','se'),('ivw_p_sensitivity','pval'),('ivw_q_sensitivity','q')]:
            rows.append(metric('MR_IVW',uid,oldcol,old[oldcol],fixed[newcol],atol))
        anchor=str(old.get('anchor_rsid',old.top_rsid));selected=data[data.key==anchor]
        if len(selected)!=1:raise ValueError('Frozen anchor unavailable/ambiguous: '+uid)
        science.write(selected[['key']],target/'anchor.tsv');science.mr(target/'harmonized.tsv',target/'anchor.tsv',target/'wald.tsv')
        wald=science.table(target/'wald.tsv').iloc[0]
        for oldcol,newcol in [('wald_beta','beta'),('wald_se','se'),('wald_p','pval')]:
            rows.append(metric('MR_WALD',uid,oldcol,old[oldcol],wald[newcol],atol))
        rows.append(metric('MR_IVW',uid,'nsnp',old.n_harmonized_instruments,len(data),0,0))
    if not rows:raise ValueError('No selected legacy MR results')
    return rows,{'harmonized':sha(inp),'baseline':sha(baseline),'ancestry':ancestry,'phenotype':phenotype}


def replay_coloc(gene,source,baseline,out,atol,policy_delta=True):
    d=legacy_rows(source);out.mkdir(parents=True,exist_ok=True)
    data=canonical_replay(d);science.write(data,out/'same_input.tsv');atomic_json(out/'metadata.json',replay_meta())
    bridge('coloc',out/'same_input.tsv',out/'metadata.json',out/'coloc.tsv')
    new=science.table(out/'coloc.tsv');old=science.table(baseline)
    old=old[old.gene_symbol==gene]
    if old.empty:raise ValueError('No baseline coloc row: '+gene)
    rows=[]
    for _,o in old.iterrows():
        prior=o.p12 if 'p12' in o else 1e-5
        match=new[np.isclose(new.p12,float(prior),rtol=0,atol=1e-12)]
        if len(match)!=1:raise ValueError('Prior baseline mismatch')
        n=match.iloc[0]
        for posterior in ['H0','H1','H2','H3','H4']:
            field='PP.'+poster
            if field in o:rows.append(metric('COLOC',gene,field,o[field],n[field+'.abf'],atol))
        if 'nsnps' in o:rows.append(metric('COLOC',gene,'nsnps',o.nsnps,n.nsnps,0,0))
    policy=[]
    if policy_delta:
        if not {'allele0_pqtl','allele1_pqtl'}.issubset(d):raise ValueError('Policy comparison requires allele columns')
        pal=(d.allele0_pqtl+d.allele1_pqtl).isin(['AT','TA','CG','GC'])
        strict=data.loc[~pal].copy()
        record={'gene_symbol':gene,'legacy_snps':len(data),'retained_snps':len(strict),'palindromic_removed':int(pal.sum()),'status':'NOT_RUN_INSUFFICIENT_SNPS'}
        if len(strict)>=50:
            science.write(strict,out/'policy_delta.tsv');bridge('coloc',out/'policy_delta.tsv',out/'metadata.json',out/'policy_coloc.tsv')
            sf=science.table(out/'policy_coloc.tsv');a=new[np.isclose(new.p12,1e-5,atol=1e-12)].iloc[0];b=sf[np.isclose(sf.p12,1e-5,atol=1e-12)].iloc[0]
            record.update(status='POLICY_COMPARISON',old_h4=a['PP.H4.abf'],strict_h4=b['PP.H4.abf'],delta_h4=b['PP.H4.abf']-a['PP.H4.abf'])
        policy.append(record)
    return rows,policy,{'input':sha(source),'baseline':sha(baseline),'filtered_rows':len(d),'same_input':sha(out/'same_input.tsv'),'N_policy':'median rounded, as legacy R caller','sdY_policy':replay_meta()}


def replay_susie(gene,source,baseline,ldroot,out,atol,max_variants=8000):
    d=legacy_rows(source);meta=ldroot/(gene+'.ld.meta.tsv');order=ldroot/(gene+'.ld.vars');binary=ldroot/(gene+'.ld.bin')
    m=science.table(meta);ids=[s for s in order.read_text().splitlines() if s]
    if list(m.ref_id)!=ids or m.snp.duplicated().any():raise ValueError('LD metadata/order mismatch')
    if len(m)>max_variants:raise ValueError('Dense LD safety limit exceeded; increase --max-ld-variants explicitly')
    if binary.stat().st_size!=len(m)*len(m)*4:raise ValueError('LD binary size mismatch')
    d=d.set_index('snp').loc[m.snp].reset_index()
    r=np.fromfile(binary,dtype='<f4').reshape(len(m),len(m)).astype(float)
    r=(r+r.T)/2;sign=m.effect_vs_ld_major_sign.to_numpy(dtype=float)
    if not np.isin(sign,[-1,1]).all():raise ValueError('LD sign invalid')
    r*=np.outer(sign,sign);np.fill_diagonal(r,1)
    if not np.isfinite(r).all() or np.abs(r).max()>1.00001:raise ValueError('LD numeric range')
    data=canonical_replay(d);out.mkdir(parents=True,exist_ok=True)
    science.write(data,out/'same_input.tsv');atomic_json(out/'metadata.json',replay_meta())
    pd.DataFrame(r,index=data.key,columns=data.key).to_csv(out/'signed_ld.tsv',sep='\t')
    bridge('susie',out/'same_input.tsv',out/'metadata.json',out/'susie.tsv',out/'signed_ld.tsv')
    fit=science.table(out/'susie.tsv');old=science.table(baseline);old=old[old.gene_symbol==gene]
    if len(old)!=1:raise ValueError('SuSiE baseline absent/ambiguous')
    o=old.iloc[0]
    if set(fit.status)!={'SUCCESS'}:
        return [{'stage':'SUSIE','gene_symbol':gene,'metric':'resolution','status':'UNRESOLVED','comparison':'SAME_INPUT'}],{'input':sha(source),'baseline':sha(baseline),'ld':sha(binary),'ld_meta':sha(meta),'ld_order':sha(order)}
    rows=[metric('SUSIE',gene,'max_PP.H4',o['max_PP.H4'],fit['PP.H4.abf'].max(),atol),metric('SUSIE',gene,'nsnps',o.nsnps,len(data),0,0)]
    for oldfield,newfield in [('pqtl_credible_sets','exposure_credible_sets'),('egfr_credible_sets','outcome_credible_sets')]:
        if oldfield in o:rows.append(metric('SUSIE',gene,oldfield,o[oldfield],fit.iloc[0][newfield],0,0))
    return rows,{'input':sha(source),'baseline':sha(baseline),'ld':sha(binary),'ld_meta':sha(meta),'ld_order':sha(order),'n_snps':len(data)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path('/srv/is-analysis'))
    parser.add_argument('--project',choices=['ckd','ischemic_stroke'],default='ckd')
    parser.add_argument('--stage1-root',type=Path);parser.add_argument('--stage2b-root',type=Path);parser.add_argument('--stage2c-root',type=Path);parser.add_argument('--ld-root',type=Path)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--genes',nargs='+',default=GENES)
    parser.add_argument('--stages',nargs='+',choices=['mr','coloc','susie'],default=['mr','coloc'])
    parser.add_argument('--ancestry',choices=['EUR','EAS'],default='EUR');parser.add_argument('--phenotype',default='eGFRcrea')
    parser.add_argument('--mr-atol',type=float,default=1e-8);parser.add_argument('--posterior-atol',type=float,default=1e-4)
    parser.add_argument('--max-ld-variants',type=int,default=8000);parser.add_argument('--no-policy-delta',action='store_true')
    a=parser.parse_args()
    if a.project=='ischemic_stroke':parser.error('This adapter validates the legacy CKD schema only; IS requires its own explicit baseline adapter')
    if a.ancestry!='EUR' and any(s in a.stages for s in ['coloc','susie']):parser.error('Legacy regional/LD adapter is EUR only; EAS replay currently supports MR only')
    base=a.root/'results'/a.project
    a.stage1_root=a.stage1_root or base/'stage1';a.stage2b_root=a.stage2b_root or base/'stage2b_coloc';a.stage2c_root=a.stage2c_root or base/'stage2c_susie';a.ld_root=a.ld_root or a.root/'data'/a.project/'stage2c_ld/ld'
    if a.out.exists():parser.error('Use a new output directory per regression run; no historical output overwrites')
    sources=[a.stage1_root,a.stage2b_root,a.stage2c_root,a.ld_root]
    if any(a.out.resolve()==p.resolve() or p.resolve() in a.out.resolve().parents for p in sources):parser.error('Output must be outside all legacy input roots')
    a.out.mkdir(parents=True);rows=[];policy=[];provenance={};errors=[]
    for stage in a.stages:
        if stage=='mr':
            try:
                rr,pp=replay_mr(a.stage1_root,a.out/'mr',a.genes,a.mr_atol,a.ancestry,a.phenotype);rows+=rr;provenance['mr']=pp
            except Exception as e:errors.append({'stage':stage,'error':str(e)})
        else:
            for gene in a.genes:
                try:
                    if stage=='coloc':
                        baseline=a.stage2b_root/'coloc_results/STAGE2B_COLOC_DEFAULT.tsv'
                        rr,dd,pp=replay_coloc(gene,a.stage2b_root/'coloc_input'/(gene+'.tsv.gz'),baseline,a.out/stage/gene,a.posterior_atol,not a.no_policy_delta);rows+=rr;policy+=dd;provenance[stage+':'+gene]=pp
                    else:
                        rr,pp=replay_susie(gene,a.stage2c_root/'susie_input'/(gene+'.tsv.gz'),a.stage2c_root/'susie_results/STAGE2C_SUSIE_DEFAULT.tsv',a.ld_root,a.out/stage/gene,a.posterior_atol,a.max_ld_variants);rows+=rr;provenance[stage+':'+gene]=pp
                except Exception as e:errors.append({'stage':stage,'gene_symbol':gene,'error':str(e)})
    if rows:science.write(pd.DataFrame(rows),a.out/'REGRESSION_COMPARISON.tsv')
    if policy:science.write(pd.DataFrame(policy),a.out/'POLICY_DELTA.tsv')
    statuses=[r['status'] for r in rows]
    status='INCOMPLETE' if errors or 'UNRESOLVED' in statuses else 'DIFFERENCE' if 'DIFFERENCE' in statuses else 'PASS' if rows else 'NO_RESULTS'
    environment={'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__}
    try:
        r=subprocess.run(['Rscript','-e','sessionInfo()'],capture_output=True,text=True,timeout=30)
        (a.out/'R-session-info.txt').write_text(r.stdout+'\n'+r.stderr);environment['R_session_exit']=r.returncode
    except (OSError,subprocess.TimeoutExpired) as e:environment['R_session_error']=str(e)
    atomic_json(a.out/'REGRESSION_SUMMARY.json',{'project':a.project,'environment':environment,'status':status,'stages':a.stages,'genes':a.genes,'metrics':len(rows),'pass_metrics':statuses.count('PASS'),'difference_metrics':statuses.count('DIFFERENCE'),'errors':errors,'provenance':provenance,'code_sha256':code_digest(CORE.parent),'timestamp':datetime.now(timezone.utc).isoformat(),'scope':'Frozen prepared-input numerical parity; not a full raw-data rebuild. Policy deltas do not count as numerical failures.'})
    print(json.dumps({'status':status,'metrics':len(rows),'errors':errors,'output':str(a.out)},indent=2));return 0 if status=='PASS' else 1
if __name__=='__main__':sys.exit(main())
