#!/usr/bin/env python3
from pathlib import Path
import csv, gzip, tarfile, json

ROOT = Path('/srv/is-analysis')
EXP = ROOT/'results/metabolic_resilience/stage2_gwas/exposure_marginal/UKBPPP_ST9_STRONGEST_CIS_PER_PROTEIN.tsv.gz'
PREFLIGHT = ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B_CIS_FILE_PREFLIGHT.tsv'
OUT = ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B1_COORDINATE_ANCHOR_AUDIT.tsv'
SUMMARY = ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B1_COORDINATE_SUMMARY.json'

def ff(x):
    try: return float(x)
    except: return None

def ii(x):
    try: return int(float(x))
    except: return None

def ss(x): return '' if x is None else str(x).strip()

anchors = {}
with gzip.open(EXP,'rt',encoding='utf-8',newline='') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        anchors[r['protein_id']] = {
            'gene_symbol': r['gene_symbol'],
            'variant_id_hg19': r['variant_id_hg19'],
            'rsid': r['rsid'],
            'chr_hg19': ii(r['chr_hg19']),
            'pos_hg19': ii(r['pos_hg19']),
            'effect_allele_st9': ss(r['effect_allele']).upper(),
            'other_allele_st9': ss(r['other_allele']).upper(),
            'eaf_st9': ff(r['eaf']),
            'beta_st9': ff(r['beta_exposure']),
        }

members = {}
with PREFLIGHT.open() as f:
    for r in csv.DictReader(f,delimiter='\t'):
        if r['status']=='PASS': members[r['protein_id']] = r

results=[]
for pid,m in members.items():
    if pid not in anchors: continue
    a=anchors[pid]
    exact=[]; posmatch=[]
    with tarfile.open(m['tar_path'],'r') as tar:
        raw=tar.extractfile(tar.getmember(m['cis_member']))
        with gzip.GzipFile(fileobj=raw,mode='rb') as gz:
            header=gz.readline().decode(errors='replace').strip().split()
            H={x:i for i,x in enumerate(header)}
            for b in gz:
                p=b.decode(errors='replace').strip().split()
                if len(p)<len(header): continue
                vid=p[H['ID']]; t=vid.split(':')
                if len(t)<2: continue
                row={'ID':vid,'CHROM':ii(p[H['CHROM']]),'GENPOS':ii(p[H['GENPOS']]),'ID_CHR':ii(t[0]),'ID_POS':ii(t[1]),
                     'ALLELE0':p[H['ALLELE0']].upper(),'ALLELE1':p[H['ALLELE1']].upper(),
                     'A1FREQ':ff(p[H['A1FREQ']]),'BETA':ff(p[H['BETA']])}
                if vid==a['variant_id_hg19']: exact.append(row)
                if row['ID_CHR']==a['chr_hg19'] and row['ID_POS']==a['pos_hg19']: posmatch.append(row)
    s=(exact or posmatch or [None])[0]
    if s is None:
        results.append({'protein_id':pid,'gene_symbol':a['gene_symbol'],'status':'FAIL_ANCHOR_NOT_FOUND'})
        continue
    direct=s['ALLELE1']==a['effect_allele_st9'] and s['ALLELE0']==a['other_allele_st9']
    swapped=s['ALLELE0']==a['effect_allele_st9'] and s['ALLELE1']==a['other_allele_st9']
    aligned_beta=s['BETA'] if direct else (-s['BETA'] if swapped and s['BETA'] is not None else None)
    aligned_eaf=s['A1FREQ'] if direct else (1-s['A1FREQ'] if swapped and s['A1FREQ'] is not None else None)
    status='PASS' if (s['ID_CHR']==a['chr_hg19'] and s['ID_POS']==a['pos_hg19'] and (direct or swapped)) else 'FAIL'
    results.append({
        'protein_id':pid,'gene_symbol':a['gene_symbol'],'st9_variant_id_hg19':a['variant_id_hg19'],
        'st9_chr_hg19':a['chr_hg19'],'st9_pos_hg19':a['pos_hg19'],'pgwas_ID':s['ID'],
        'pgwas_ID_chr':s['ID_CHR'],'pgwas_ID_pos':s['ID_POS'],'pgwas_GENPOS':s['GENPOS'],
        'GENPOS_minus_IDPOS':(s['GENPOS']-s['ID_POS']) if s['GENPOS'] and s['ID_POS'] else '',
        'exact_ID_match':int(s['ID']==a['variant_id_hg19']),
        'ID_pos_matches_ST9_hg19':int(s['ID_CHR']==a['chr_hg19'] and s['ID_POS']==a['pos_hg19']),
        'GENPOS_matches_ST9_hg19':int(s['CHROM']==a['chr_hg19'] and s['GENPOS']==a['pos_hg19']),
        'allele_orientation_match':'DIRECT' if direct else ('SWAPPED' if swapped else 'ALLELE_MISMATCH'),
        'beta_st9':a['beta_st9'],'beta_pgwas':aligned_beta,
        'beta_abs_diff':abs(a['beta_st9']-aligned_beta) if a['beta_st9'] is not None and aligned_beta is not None else '',
        'eaf_st9':a['eaf_st9'],'eaf_pgwas':aligned_eaf,
        'eaf_abs_diff':abs(a['eaf_st9']-aligned_eaf) if a['eaf_st9'] is not None and aligned_eaf is not None else '',
        'status':status,
    })

fields=[]
for r in results:
    for k in r:
        if k not in fields: fields.append(k)
with OUT.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(results)
summary={
    'targets':len(results),
    'pass':sum(r.get('status')=='PASS' for r in results),
    'exact_ID_match':sum(r.get('exact_ID_match')==1 for r in results),
    'ID_position_matches_ST9_hg19':sum(r.get('ID_pos_matches_ST9_hg19')==1 for r in results),
    'GENPOS_matches_ST9_hg19':sum(r.get('GENPOS_matches_ST9_hg19')==1 for r in results),
    'cis_coordinate_field':'position embedded in ID',
    'cis_coordinate_build':'GRCh37/hg19',
}
SUMMARY.write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
