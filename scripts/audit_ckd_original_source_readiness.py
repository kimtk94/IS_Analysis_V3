#!/usr/bin/env python3
"""Read-only original-source readiness and selected variant trace for CKD.
Only searches known research locations and files; never fetches or mutates input.
"""
import csv,gzip,json,os
from pathlib import Path
TARGET={"rs35646972","rs6663222","rs11834807","rs11860566"}
ROOT=Path("/srv/is-analysis")
SEARCH=[ROOT/"results/ckd/stage1",ROOT/"results/ckd/stage2",ROOT/"results/ckd/stage2b_coloc",ROOT/"data/ckd"]
TOKENS=("pqtl","gwas","egfr","harmon","ukb","source")
def scan():
    files=[]
    for root in SEARCH:
        if not root.is_dir():continue
        for base,dirs,names in os.walk(root):
            depth=len(Path(base).relative_to(root).parts)
            dirs[:]=[d for d in dirs if not d.startswith(".")] if depth<5 else []
            for name in names:
                low=name.lower()
                if not any(t in low for t in TOKENS):continue
                if not low.endswith((".tsv",".tsv.gz",".csv",".csv.gz",".txt",".txt.gz")):continue
                p=Path(base)/name
                try: size=p.stat().st_size
                except OSError:continue
                files.append((str(p),size))
            if len(files)>5000:return files[:5000]
    return files
def main():
    files=scan()
    bytype={"pqtl":[],"gwas":[],"harmonized":[]}
    for p,size in files:
        low=p.lower()
        kind="harmonized" if "harmon" in low or "stage2b_coloc" in low else "pqtl" if "pqtl" in low or "ukb" in low else "gwas"
        if len(bytype[kind])<30:bytype[kind].append({"path":p,"size":size})
    print(json.dumps({"audit":"CKD_TOP3_ORIGINAL_SOURCE_DISCOVERY","files_scanned":len(files),
        "candidate_sources_by_class":bytype,
        "original_consortium_allele_verified":False,
        "original_consortium_maf_verified":False,
        "note":"File discovery does not imply original unharmonized consortium data present or comparable."},indent=2))
if __name__=="__main__":main()
