#!/usr/bin/env bash
ROOT="${IS_ANALYSIS_ROOT:-/srv/is-analysis}"
cd "$ROOT" || exit 1
set +e
set +u
set +o pipefail 2>/dev/null || true

TARGETS="$ROOT/results/metabolic_resilience/stage2_gwas/final_shortlist/STAGE3_FULL_PGWAS_TARGETS.txt"
MATCH="$ROOT/results/metabolic_resilience/stage2_gwas/exposure_marginal/UKBPPP_TARGET_SYN_MATCH.tsv"
VENV_PY="$ROOT/.venv-ukbppp/bin/python"
DATA="$ROOT/data/metabolic_resilience/stage3_full_pgwas/ukbppp_eur"
AUDIT="$ROOT/results/metabolic_resilience/stage3_full_pgwas/audit"
MAP="$AUDIT/STAGE3A_CANDIDATE_SYN_MAP.tsv"
DOWNLOAD_MANIFEST="$AUDIT/STAGE3A_DOWNLOAD_MANIFEST.tsv"
mkdir -p "$DATA" "$AUDIT"

echo "===================================================="
echo "STAGE 3-A : CANDIDATE-ONLY UKB-PPP FULL pGWAS"
echo "===================================================="

READY=1
for F in "$TARGETS" "$MATCH" "$VENV_PY"; do
  if [ -s "$F" ]; then echo "[PASS] $F"; else echo "[FAIL] $F"; READY=0; fi
done

if [ "$READY" = "1" ]; then
python3 - <<'PY'
from pathlib import Path
import csv,re
ROOT=Path('/srv/is-analysis')
TARGETS=ROOT/'results/metabolic_resilience/stage2_gwas/final_shortlist/STAGE3_FULL_PGWAS_TARGETS.txt'
MATCH=ROOT/'results/metabolic_resilience/stage2_gwas/exposure_marginal/UKBPPP_TARGET_SYN_MATCH.tsv'
OUT=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3A_CANDIDATE_SYN_MAP.tsv'
targets=[x.strip() for x in TARGETS.read_text().splitlines() if x.strip()]
with MATCH.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
out=[]
for pid in targets:
    parts=pid.split(':'); oid=parts[2] if len(parts)>2 else ''
    mm=[]
    for r in rows:
        blob='\t'.join(str(v) for v in r.values() if v is not None)
        if oid and oid in blob: mm.append(r)
    syn=set(); fn=set()
    for r in mm:
        for v in r.values():
            v=str(v or '').strip()
            if re.fullmatch(r'syn\d+',v,re.I): syn.add(v)
            if any(x in v.lower() for x in ['.tar','.tgz','.gz']): fn.add(v)
    status='PASS' if len(mm)>0 and len(syn)==1 else ('FAIL_NO_MATCH' if len(mm)==0 else 'FAIL_SYN_ID')
    out.append({'protein_id':pid,'gene_symbol':parts[0] if parts else '','uniprot':parts[1] if len(parts)>1 else '',
                'oid':oid,'n_manifest_rows':len(mm),'synapse_id':sorted(syn)[0] if len(syn)==1 else '',
                'manifest_filenames':';'.join(sorted(fn)),'status':status})
with OUT.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=out[0].keys(),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
print('targets=',len(out),'PASS=',sum(x['status']=='PASS' for x in out),'FAIL=',sum(x['status']!='PASS' for x in out))
PY
MAP_RC=$?
else
MAP_RC=99
fi

N_PASS=0
if [ "$MAP_RC" = "0" ] && [ -s "$MAP" ]; then
  N_PASS="$(awk -F'\t' 'NR>1 && $NF=="PASS" {n++} END {print n+0}' "$MAP")"
fi
echo "Mapped PASS=$N_PASS / 8"

if [ "$N_PASS" = "8" ]; then
"$VENV_PY" - <<'PY'
from pathlib import Path
import csv,synapseclient
ROOT=Path('/srv/is-analysis')
MAP=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3A_CANDIDATE_SYN_MAP.tsv'
DATA=ROOT/'data/metabolic_resilience/stage3_full_pgwas/ukbppp_eur'
OUT=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3A_DOWNLOAD_MANIFEST.tsv'
with MAP.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
syn=synapseclient.Synapse(); syn.login(silent=True)
out=[]
for r in rows:
    d=DATA/r['oid']; d.mkdir(parents=True,exist_ok=True)
    try:
        e=syn.get(r['synapse_id'],downloadLocation=str(d)); p=Path(e.path)
        st='PASS' if p.exists() and p.stat().st_size>0 else 'FAIL_EMPTY'
        out.append({'protein_id':r['protein_id'],'oid':r['oid'],'synapse_id':r['synapse_id'],'local_path':str(p),
                    'size_bytes':p.stat().st_size if p.exists() else 0,'status':st})
    except Exception as ex:
        print('[FAIL]',r['protein_id'],type(ex).__name__,str(ex))
        out.append({'protein_id':r['protein_id'],'oid':r['oid'],'synapse_id':r['synapse_id'],'local_path':'','size_bytes':0,'status':'FAIL_DOWNLOAD'})
with OUT.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=out[0].keys(),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
print('PASS=',sum(x['status']=='PASS' for x in out),'/',len(out))
PY
fi

echo "===== STAGE 3-A END ====="
