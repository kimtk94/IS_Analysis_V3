#!/usr/bin/env bash
ROOT="${IS_ANALYSIS_ROOT:-/srv/is-analysis}"
cd "$ROOT" || exit 1
set +e
set +u
set +o pipefail 2>/dev/null || true

PY="$ROOT/.venv-ukbppp/bin/python"
DATA="$ROOT/data/metabolic_resilience/stage3_full_pgwas/ukbppp_eur"
AUDIT="$ROOT/results/metabolic_resilience/stage3_full_pgwas/audit"
MANIFEST="$AUDIT/STAGE3A_DOWNLOAD_MANIFEST.tsv"
OID="OID21468"
DIR="$DATA/$OID"
TARGET="$DIR/IDUA_P35475_OID21468_v1_Oncology.tar"
EXPECTED_MD5="fc2c6a199fd3fd49422082ebe787fc0d"
mkdir -p "$DIR" "$AUDIT"

find "$DIR" -maxdepth 1 -type f \( -name '*.synapse_download_*' -o -name 'IDUA_P35475_OID21468_v1_Oncology.tar' \) -print -delete

IDUA_OK=0
for ATTEMPT in 1 2 3; do
  echo "IDUA DOWNLOAD ATTEMPT $ATTEMPT / 3"
  "$PY" - <<'PY'
from pathlib import Path
import synapseclient
DEST=Path('/srv/is-analysis/data/metabolic_resilience/stage3_full_pgwas/ukbppp_eur/OID21468')
s=synapseclient.Synapse(); s.login(silent=True)
e=s.get(SID,downloadLocation=str(DEST),ifcollision='overwrite.local')
print(e.path)
PY
  if [ -s "$TARGET" ]; then
    MD5="$(md5sum "$TARGET" | awk '{print $1}')"
    echo "expected=$EXPECTED_MD5 actual=$MD5"
    if [ "$MD5" = "$EXPECTED_MD5" ]; then IDUA_OK=1; break; fi
    rm -f "$TARGET"
  fi
done

if [ "$IDUA_OK" = "1" ]; then
  tar -tf "$TARGET" >/dev/null 2>&1
  TAR_RC=$?
  echo "tar rc=$TAR_RC"
else
  TAR_RC=99
fi

if [ "$IDUA_OK" = "1" ] && [ "$TAR_RC" = "0" ] && [ -s "$MANIFEST" ]; then
python3 - <<'PY'
from pathlib import Path
import csv
p=Path('/srv/is-analysis/results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3A_DOWNLOAD_MANIFEST.tsv')
t=Path('/srv/is-analysis/data/metabolic_resilience/stage3_full_pgwas/ukbppp_eur/OID21468/IDUA_P35475_OID21468_v1_Oncology.tar')
with p.open() as f: rd=csv.DictReader(f,delimiter='\t'); fields=rd.fieldnames; rows=list(rd)
for r in rows:
    if r['synapse_id']=='syn51470923':
        r['local_path']=str(t);r['size_bytes']=str(t.stat().st_size);r['status']='PASS'
with p.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
print('[PASS] manifest updated')
PY
fi

echo "IDUA_OK=$IDUA_OK"
