#!/usr/bin/env python3
from pathlib import Path
import gzip
import zipfile
import json
import os

ROOT=Path(os.environ.get("IS_ANALYSIS_ROOT","/srv/is-analysis"))
DATA=ROOT/"data/metabolic_resilience/stage2_gwas"

TARGETS=[
    ("lipids/eur",["*.gz","*.tsv","*.txt"]),
    ("blood_pressure/eur",["*.gz","*.tsv","*.txt"]),
    ("adiposity/eur",["*.gz","*.tsv","*.txt"]),
    ("glycemia/eur",["*.gz","*.tsv","*.txt"]),
    ("t2d/eur",["*.gz","*.zip","*.tsv","*.txt"]),
]

KEYWORDS=("hdl","tg","trig","bmi","whr","waist","sbp","dbp","gcst","mahajan","t2d")

def header(path):
    try:
        if path.suffix.lower()==".zip":
            with zipfile.ZipFile(path) as z:
                members=[x for x in z.namelist() if not x.endswith("/")]
                if not members:
                    return ""
                name=members[0]
                with z.open(name) as f:
                    for b in f:
                        s=b.decode("utf-8",errors="replace").strip()
                        if s and not s.startswith("##"):
                            return s[:500]
                return ""
        opener=gzip.open if str(path).endswith(".gz") else open
        with opener(path,"rt",encoding="utf-8",errors="replace") as f:
            for line in f:
                s=line.strip()
                if s and not s.startswith("##"):
                    return s[:500]
    except Exception as e:
        return f"<HEADER_ERROR {type(e).__name__}: {e}>"
    return ""

rows=[]
for rel,patterns in TARGETS:
    base=DATA/rel
    if not base.exists():
        rows.append({"dir":rel,"file":"<DIR_MISSING>","size":0,"header":""})
        continue
    seen={}
    for pat in patterns:
        for p in base.glob(pat):
            if p.is_file():
                seen[str(p.resolve())]=p
    files=sorted(seen.values(),key=lambda p:p.name.lower())
    for p in files:
        low=p.name.lower()
        interesting=any(k in low for k in KEYWORDS) or rel in ("adiposity/eur","lipids/eur")
        if not interesting:
            continue
        rows.append({
            "dir":rel,
            "file":p.name,
            "size":p.stat().st_size,
            "header":header(p),
        })

print("="*100)
print("STAGE 3-D0 OUTCOME INVENTORY")
print("="*100)
for r in rows:
    print()
    print(f"[{r['dir']}] {r['file']} size={r['size']}")
    print("HEADER:",r["header"])

out=ROOT/"results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3D0_OUTCOME_INVENTORY.json"
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(rows,indent=2)+"\n",encoding="utf-8")
print()
print("inventory_json =",out)
