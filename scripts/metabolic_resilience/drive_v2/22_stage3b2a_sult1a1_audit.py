#!/usr/bin/env python3
from pathlib import Path
import csv, tarfile, gzip, math, json

ROOT=Path('/srv/is-analysis')
WINDOWS=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B2_GENE_CIS_WINDOWS.tsv'
PREFLIGHT=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B_CIS_FILE_PREFLIGHT.tsv'
OUT=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B2A_SULT1A1_AUDIT.tsv'
SUMMARY=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B2A_SULT1A1_SUMMARY.json'
PID='SULT1A1:P50225:OID21031:v1'; LOG10P_THRESHOLD=-math.log10(1.7e-11)

window=None
with WINDOWS.open() as f:
    for r in csv.DictReader(f,delimiter='\t'):
        if r['protein_id']==PID: window=r; break
if window is None: raise SystemExit(20)
chrom=int(window['chrom']); gene_start=int(window['gene_start']); gene_end=int(window['gene_end']); cis_start=int(window['cis_start']); cis_end=int(window['cis_end']); anchor=int(window['st9_pos'])
distance=0 if gene_start<=anchor<=gene_end else min(abs(anchor-gene_start),abs(anchor-gene_end))

pf=None
with PREFLIGHT.open() as f:
    for r in csv.DictReader(f,delimiter='\t'):
        if r['protein_id']==PID: pf=r; break
if pf is None: raise SystemExit(21)

rows=[]; n_scanned=n_cis=n_sig=n_f10=0
with tarfile.open(pf['tar_path'],'r') as tar:
    raw=tar.extractfile(tar.getmember(pf['cis_member']))
    with gzip.GzipFile(fileobj=raw,mode='rb') as gz:
        header=gz.readline().decode(errors='replace').strip().split(); H={x:i for i,x in enumerate(header)}
        for b in gz:
            n_scanned+=1; p=b.decode(errors='replace').strip().split()
            if len(p)<len(header): continue
            vid=p[H['ID']]; tok=vid.split(':')
            if len(tok)<4: continue
            try: id_chr=int(tok[0]); id_pos=int(tok[1])
            except: continue
            if id_chr!=chrom or not (cis_start<=id_pos<=cis_end): continue
            n_cis+=1
            try:
                beta=float(p[H['BETA']]); se=float(p[H['SE']]); log10p=float(p[H['LOG10P']]); eaf=float(p[H['A1FREQ']]); info=float(p[H['INFO']])
            except: continue
            F=(beta/se)**2 if se>0 else None; sig=log10p>=LOG10P_THRESHOLD
            if sig: n_sig+=1
            if sig and F is not None and F>10: n_f10+=1
            rows.append({'protein_id':PID,'gene_symbol':'SULT1A1','ID':vid,'chrom_hg19':id_chr,'pos_hg19':id_pos,'ALLELE0':p[H['ALLELE0']],
                         'ALLELE1':p[H['ALLELE1']],'A1FREQ':eaf,'INFO':info,'N':p[H['N']],'BETA':beta,'SE':se,'F':F,'LOG10P':log10p,
                         'P_approx':10**(-log10p) if log10p<300 else 0,'p_lt_1p7e_11':int(sig),'F_gt_10':int(F is not None and F>10)})
rows.sort(key=lambda x:x['LOG10P'],reverse=True)
with OUT.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
s={'protein_id':PID,'gene_interval_grch37':[gene_start,gene_end],'conventional_cis_window':[cis_start,cis_end],'st9_anchor':anchor,
   'anchor_distance_from_gene_bp':distance,'anchor_outside_1Mb_by_bp':max(0,distance-1_000_000),'pgwas_rows_scanned':n_scanned,
   'variants_in_true_cis':n_cis,'variants_p_lt_1p7e_11':n_sig,'variants_p_lt_1p7e_11_and_F_gt_10':n_f10,
   'strongest_true_cis_variant':rows[0]['ID'] if rows else None,'strongest_true_cis_LOG10P':rows[0]['LOG10P'] if rows else None,'strongest_true_cis_F':rows[0]['F'] if rows else None}
SUMMARY.write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s,indent=2))
