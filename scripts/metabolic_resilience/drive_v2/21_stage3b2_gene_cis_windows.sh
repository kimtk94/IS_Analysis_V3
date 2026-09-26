#!/usr/bin/env bash
cd /srv/is-analysis
set +e
set +u
set +o pipefail 2>/dev/null || true

ROOT="/srv/is-analysis"
SHORT="$ROOT/results/metabolic_resilience/stage2_gwas/final_shortlist/STAGE2_FINAL_SCREENING_SHORTLIST.tsv"
EXP="$ROOT/results/metabolic_resilience/stage2_gwas/exposure_marginal/UKBPPP_ST9_STRONGEST_CIS_PER_PROTEIN.tsv.gz"
OUTDIR="$ROOT/results/metabolic_resilience/stage3_full_pgwas/audit"
PAYLOAD="$OUTDIR/STAGE3B2_ENSEMBL_PAYLOAD.json"
JSON="$OUTDIR/STAGE3B2_ENSEMBL_GRCH37_LOOKUP.json"
OUT="$OUTDIR/STAGE3B2_GENE_CIS_WINDOWS.tsv"
SUMMARY="$OUTDIR/STAGE3B2_GENE_CIS_WINDOWS_SUMMARY.json"
mkdir -p "$OUTDIR"

python3 - <<'PY'
import csv,json
from pathlib import Path
ROOT=Path('/srv/is-analysis')
p=ROOT/'results/metabolic_resilience/stage2_gwas/final_shortlist/STAGE2_FINAL_SCREENING_SHORTLIST.tsv'
out=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B2_ENSEMBL_PAYLOAD.json'
with p.open() as f: genes=[r['gene_symbol'] for r in csv.DictReader(f,delimiter='\t')]
out.write_text(json.dumps({'symbols':genes}))
print(out.read_text())
PY

curl -sS -L --retry 3 --retry-delay 2 --connect-timeout 20 --max-time 120 \
  -X POST -H 'Content-Type: application/json' -H 'Accept: application/json' \
  --data-binary @"$PAYLOAD" \
  "https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens" -o "$JSON"
CURL_RC=$?
echo "curl rc=$CURL_RC"

if [ "$CURL_RC" = "0" ] && [ -s "$JSON" ]; then
python3 - <<'PY'
from pathlib import Path
import csv,gzip,json
ROOT=Path('/srv/is-analysis')
SHORT=ROOT/'results/metabolic_resilience/stage2_gwas/final_shortlist/STAGE2_FINAL_SCREENING_SHORTLIST.tsv'
EXP=ROOT/'results/metabolic_resilience/stage2_gwas/exposure_marginal/UKBPPP_ST9_STRONGEST_CIS_PER_PROTEIN.tsv.gz'
JSON=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B2_ENSEMBL_GRCH37_LOOKUP.json'
OUT=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B2_GENE_CIS_WINDOWS.tsv'
SUMMARY=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B2_GENE_CIS_WINDOWS_SUMMARY.json'
with SHORT.open() as f: short=list(csv.DictReader(f,delimiter='\t'))
anchors={}
with gzip.open(EXP,'rt') as f:
    for r in csv.DictReader(f,delimiter='\t'): anchors[r['protein_id']]=r
data=json.loads(JSON.read_text())
rows=[]
for r in short:
    pid=r['protein_id'];gene=r['gene_symbol'];a=anchors.get(pid,{});e=data.get(gene)
    if not e: continue
    chrom=str(e.get('seq_region_name','')); start=int(e['start']); end=int(e['end']); cis_start=max(1,start-1_000_000); cis_end=end+1_000_000
    st9_chr=str(a.get('chr_hg19','')).replace('chr','')
    try: st9_pos=int(a.get('pos_hg19'))
    except: st9_pos=None
    same=chrom==st9_chr; inside_gene=same and st9_pos is not None and start<=st9_pos<=end; inside_cis=same and st9_pos is not None and cis_start<=st9_pos<=cis_end
    status='PASS'
    if str(e.get('assembly_name',''))!='GRCh37': status='FAIL_NOT_GRCH37'
    elif not same: status='FAIL_CHROM_MISMATCH'
    elif not inside_cis: status='FAIL_ST9_OUTSIDE_CIS'
    elif str(e.get('biotype',''))!='protein_coding': status='WARN_NON_PROTEIN_CODING'
    rows.append({'protein_id':pid,'gene_symbol':gene,'ensembl_gene_id':e.get('id',''),'assembly':e.get('assembly_name',''),'biotype':e.get('biotype',''),
                 'chrom':chrom,'gene_start':start,'gene_end':end,'strand':e.get('strand',''),'cis_start':cis_start,'cis_end':cis_end,
                 'st9_chr':st9_chr,'st9_pos':st9_pos,'st9_inside_gene':int(inside_gene),'st9_inside_cis_1mb':int(inside_cis),'status':status})
with OUT.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
s={'targets':len(rows),'pass':sum(x['status']=='PASS' for x in rows),'warnings':sum(x['status'].startswith('WARN') for x in rows),'fail':sum(x['status'].startswith('FAIL') for x in rows),
   'GRCh37':sum(x['assembly']=='GRCh37' for x in rows),'protein_coding':sum(x['biotype']=='protein_coding' for x in rows),'ST9_inside_gene':sum(x['st9_inside_gene']==1 for x in rows),
   'ST9_inside_cis_1Mb':sum(x['st9_inside_cis_1mb']==1 for x in rows),'cis_definition':'Ensembl GRCh37 gene interval extended 1 Mb upstream and downstream','pgwas_coordinate':'position embedded in UKB-PPP ID field'}
SUMMARY.write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s,indent=2))
PY
fi
