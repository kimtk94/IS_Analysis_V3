#!/usr/bin/env python3
from pathlib import Path
import csv,tarfile,gzip,math,json

ROOT=Path('/srv/is-analysis')
WINDOWS=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B2_GENE_CIS_WINDOWS.tsv'
PREFLIGHT=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B_CIS_FILE_PREFLIGHT.tsv'
OUTDIR=ROOT/'results/metabolic_resilience/stage3_full_pgwas/cis_marginal'
SUMMARY=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B3_CIS_EXTRACTION_SUMMARY.tsv'
JSON=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B3_CIS_EXTRACTION_SUMMARY.json'
OUTDIR.mkdir(parents=True,exist_ok=True)

windows={}
with WINDOWS.open() as f:
    for r in csv.DictReader(f,delimiter='\t'):
        if r['assembly']=='GRCh37' and all(r.get(x) for x in ['chrom','gene_start','gene_end','cis_start','cis_end']): windows[r['protein_id']]=r
members={}
with PREFLIGHT.open() as f:
    for r in csv.DictReader(f,delimiter='\t'):
        if r['status']=='PASS': members[r['protein_id']]=r

LOG_GWS=-math.log10(5e-8); LOG_SUN=-math.log10(1.7e-11); summary=[]
for pid in sorted(set(windows)&set(members)):
    w=windows[pid]; m=members[pid]; gene=w['gene_symbol']; chrom=int(w['chrom']); gs=int(w['gene_start']); ge=int(w['gene_end']); cs=int(w['cis_start']); ce=int(w['cis_end'])
    rows=[]; scanned=0
    with tarfile.open(m['tar_path'],'r') as tar:
        raw=tar.extractfile(tar.getmember(m['cis_member']))
        with gzip.GzipFile(fileobj=raw,mode='rb') as gz:
            header=gz.readline().decode(errors='replace').strip().split(); H={x:i for i,x in enumerate(header)}
            for b in gz:
                scanned+=1; p=b.decode(errors='replace').strip().split()
                if len(p)<len(header): continue
                vid=p[H['ID']]; tok=vid.split(':')
                if len(tok)<4: continue
                try: c=int(tok[0]); pos=int(tok[1])
                except: continue
                if c!=chrom or not(cs<=pos<=ce): continue
                try:
                    beta=float(p[H['BETA']]); se=float(p[H['SE']]); eaf=float(p[H['A1FREQ']]); info=float(p[H['INFO']]); n=float(p[H['N']]); chisq=float(p[H['CHISQ']]); logp=float(p[H['LOG10P']])
                except: continue
                if se<=0: continue
                F=(beta/se)**2
                rows.append({'protein_id':pid,'gene_symbol':gene,'chrom_hg19':c,'pos_hg19':pos,'variant_id':vid,'allele0':p[H['ALLELE0']],
                             'allele1':p[H['ALLELE1']],'effect_allele':p[H['ALLELE1']],'other_allele':p[H['ALLELE0']],'eaf':eaf,'info':info,'n':n,
                             'test':p[H['TEST']],'beta':beta,'se':se,'chisq':chisq,'f_stat':F,'log10p':logp,'p_approx':10**(-logp) if logp<300 else 0.0,
                             'info_ge_0p8':int(info>=0.8),'p_lt_5e_8':int(logp>=LOG_GWS),'p_lt_1p7e_11':int(logp>=LOG_SUN),'f_gt_10':int(F>10),
                             'gene_start':gs,'gene_end':ge,'cis_start':cs,'cis_end':ce})
    rows.sort(key=lambda r:r['log10p'],reverse=True)
    outfile=OUTDIR/f"{gene}_{pid.split(':')[2]}_cis1Mb_marginal.tsv.gz"
    if rows:
        with gzip.open(outfile,'wt',newline='') as f:
            wr=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t',lineterminator='\n');wr.writeheader();wr.writerows(rows)
    n_gws=sum(r['info_ge_0p8']==1 and r['p_lt_5e_8']==1 and r['f_gt_10']==1 for r in rows)
    n_sun=sum(r['info_ge_0p8']==1 and r['p_lt_1p7e_11']==1 and r['f_gt_10']==1 for r in rows)
    strongest=rows[0] if rows else None
    summary.append({'protein_id':pid,'gene_symbol':gene,'chrom':chrom,'gene_start':gs,'gene_end':ge,'cis_start':cs,'cis_end':ce,'chromosome_rows_scanned':scanned,
                    'cis_variants':len(rows),'info_ge_0p8':sum(r['info_ge_0p8']==1 for r in rows),'gws_F10_info':n_gws,'sun_threshold_F10_info':n_sun,
                    'strongest_variant':strongest['variant_id'] if strongest else '','strongest_beta':strongest['beta'] if strongest else '',
                    'strongest_se':strongest['se'] if strongest else '','strongest_F':strongest['f_stat'] if strongest else '',
                    'strongest_LOG10P':strongest['log10p'] if strongest else '','strongest_INFO':strongest['info'] if strongest else '',
                    'output_file':str(outfile),'status':'PASS' if rows else 'FAIL_NO_CIS_VARIANTS'})
with SUMMARY.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=summary[0].keys(),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summary)
js={'targets':len(summary),'pass':sum(x['status']=='PASS' for x in summary),'all_have_GWS_instruments':sum(x['gws_F10_info']>0 for x in summary),
    'all_have_Sun_threshold_instruments':sum(x['sun_threshold_F10_info']>0 for x in summary),'coordinate_rule':'GRCh37/hg19 position embedded in UKB-PPP variant ID',
    'cis_definition':'Ensembl GRCh37 target gene interval ±1 Mb','effect_allele':'ALLELE1','exposure_effect':'marginal BETA',
    'primary_preclump_threshold':'P<5e-8, F>10, INFO>=0.8','stringent_preclump_threshold':'P<1.7e-11, F>10, INFO>=0.8'}
JSON.write_text(json.dumps(js,indent=2)+'\n');print(json.dumps(js,indent=2))
