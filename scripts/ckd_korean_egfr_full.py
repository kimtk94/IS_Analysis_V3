#!/usr/bin/env python3
"""Full publisher-verified Korean eGFR scan for frozen discovery variants."""
import argparse
import collections
import gzip
import json
import math
import re
from pathlib import Path
from ckd_korean_validation import anchors, digest, orient, summarize, write

EXPECTED_SIZE=720566637
EXPECTED_MD5="aa0056c580cc28ac77140c62dd87e079"
EXPECTED_SHA_META="7bc909e64138a0abfdb529e74f56a9cef1218ba60ee6066d23fe158eaf7caf64"
def getcol(header,*names):
    h=[re.sub(r"[^a-z0-9]","",s.lower()) for s in header]
    inds=[i for i,key in enumerate(h) if key in names]
    if len(inds)!=1:
        raise ValueError("Unexpected source column identity "+str(names)+" header "+repr(header))
    return inds[0]

def run(root,out):
    panel=anchors(root)
    meta=json.loads((root/"zenodo.json").read_text())
    if str(meta["id"])!="16976585" or meta["metadata"]["doi"]!="10.5281/zenodo.16976585" or meta["metadata"]["access_right"]!="open":
        raise ValueError("Zenodo identity/access changed")
    fs=[x for x in meta["files"] if x["key"]=="GWAS_summary_statistics.txt.gz" and x["size"]==EXPECTED_SIZE and x["checksum"]=="md5:"+EXPECTED_MD5]
    if len(fs)!=1:
        raise ValueError("Publisher file metadata changed")
    desc=meta["metadata"]["description"]
    if not all(k in desc for k in ("POS_hg19","Effect_Allele","Obeta","Effect size for the effect allele")):
        raise ValueError("Allele/build documentation changed")
    path=root/"koges.gz"
    if path.stat().st_size!=EXPECTED_SIZE or digest(path,"md5")!=EXPECTED_MD5:
        raise ValueError("Publisher original size/MD5 mismatch")
    targets={(v[0]["chr37"],v[0]["pos37"]):k for k,v in panel.items()}
    matches=collections.defaultdict(list)
    n=bad=0
    with gzip.open(path,"rt",encoding="utf8") as f:
        h=f.readline().strip().split()
        ix={"chr":getcol(h,"chrom"),"pos":getcol(h,"poshg19"),"id":getcol(h,"id"),
            "ea":getcol(h,"effectallele"),"oa":getcol(h,"otherallele"),
            "eaf":getcol(h,"effectallelefreq"),"beta":getcol(h,"obeta"),
            "se":getcol(h,"se"),"p":getcol(h,"p")}
        # Full decompression to EOF verifies gzip CRC; all nonempty rows have the source schema.
        for line in f:
            if not line.strip():continue
            n+=1
            fields=line.split()
            if len(fields)!=len(h):
                bad+=1
                continue
            try:
                coord=(int(fields[ix["chr"]].removeprefix("chr")),int(fields[ix["pos"]]))
            except ValueError:
                bad+=1
                continue
            rs=targets.get(coord)
            if rs:
                q={k:fields[j] for k,j in ix.items()}
                if q["id"].split("_")[0]==rs:
                    matches[rs].append(q)
    if n<1000000 or bad:
        raise ValueError("Unexpected full-file rows or schema: "+repr((n,bad,h)))
    result=[]
    for rs,rr in panel.items():
        a=rr[0]
        valid=[q for q in matches[rs] if orient(a["ea37"],a["oa37"],q["ea"].upper(),q["oa"].upper()) is not None]
        if len(valid)!=1:
            result.append({"rsid":rs,"status":"HOLD_UNMATCHED_OR_AMBIGUOUS","source_rows":len(matches[rs])})
            continue
        q=valid[0]
        by,sy,p,eaf=map(float,(q["beta"],q["se"],q["p"],q["eaf"]))
        if not all(math.isfinite(x) for x in (by,sy,p,eaf)) or sy<=0 or not 0<=p<=1 or not 0<eaf<1:
            result.append({"rsid":rs,"status":"HOLD_SOURCE_NUMERIC"})
            continue
        s=orient(a["ea37"],a["oa37"],q["ea"].upper(),q["oa"].upper())
        aligned_eaf=eaf if s==1 else 1-eaf
        # Fixed outcome association primary P is provider P; normal approximation is a separate diagnostic.
        r=summarize(rs,rr,by,sy,q["ea"].upper(),q["oa"].upper(),"KoGES 2026 eGFRcrea")
        r["status"]="KOREAN_OUTCOME_WALD"
        r["source_p"]=p
        r["source_bonferroni_p_7"]=min(1,p*7)
        r["effect_allele_frequency"]=eaf
        r["aligned_eaf_minus_ckb_eaf"]=aligned_eaf-a["eaf"]
        r["passes_fixed_direction_and_source_bonferroni"]=r["egfr_direction_agrees"] and p*7<.05
        result.append(r)
    write(out,{"status":"PUBLIC_KOGES_2026_FULL_OUTCOME_SCAN","source_md5":EXPECTED_MD5,"source_sha256":digest(path),
               "source_size":EXPECTED_SIZE,"source_rows":n,"source_header":h,"malformed":bad,"family":sorted(panel),
               "matched":sum(r["status"]=="KOREAN_OUTCOME_WALD" for r in result),"results":result,
               "limits":["Same CKB exposure, separate KoGES outcome cohort; individual overlap not verified.",
                         "Korean normalized eGFR phenotype; effect magnitudes cannot be interpreted as raw eGFR units.",
                         "No CKB regional outcome colocalization, no protein intervention or independent pQTL validation.",
                         "Seven preselected variant-level P tests, not ten independent assay tests."]})
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.root,a.out)
