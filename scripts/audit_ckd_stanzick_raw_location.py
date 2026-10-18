#!/usr/bin/env python3
"""Read-only search for original Stanzick 2021 eGFRcrea input."""
import os,json
from pathlib import Path
root=Path("/srv/is-analysis/data/ckd")
needle="metal_eGFR_meta_ea1.TBL.map.annot.gc.gz"
found=[]
for current,dirs,files in os.walk(root):
 depth=len(Path(current).relative_to(root).parts)
 if depth>=9:dirs.clear()
 for name in files:
  if name==needle or "stanzick" in name.lower():
   p=Path(current)/name
   found.append({"path":str(p),"size":p.stat().st_size})
print(json.dumps({"name":needle,"found":found},indent=2))
