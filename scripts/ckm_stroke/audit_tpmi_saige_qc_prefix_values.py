#!/usr/bin/env python3
"""Source bounded diagnostic of SAIGE TPMI 433.21 QC and INFO conventions.

First 1MiB source-only; record counts of QC labels and imputationInfo bounds.
No per-variant records or effects exported, no full-file extrapolation.
"""
import argparse,collections,json,math,urllib.request,zlib
from pathlib import Path
from scripts.ckm_stroke.sample_tpmi_43321_saige_casecontrol import URL,SOURCE_SIZE,PARTIAL_BYTES,EXPECTED
def audit():
 req=urllib.request.Request(URL,headers={"Range":f"bytes=0-{PARTIAL_BYTES-1}","User-Agent":"IS-AnalysisV3-TPMI43321-QC-convention/1.0"})
 with urllib.request.urlopen(req,timeout=75) as r:
  if r.status!=206 or r.geturl()!=URL or r.headers.get("Content-Range")!=f"bytes 0-{PARTIAL_BYTES-1}/{SOURCE_SIZE}":raise ValueError("not exact source range")
  compressed=r.read(PARTIAL_BYTES+1)
 if len(compressed)!=PARTIAL_BYTES:raise ValueError("unexpected source length")
 unpack=zlib.decompressobj(zlib.MAX_WBITS|16).decompress(compressed,20000000)
 if len(unpack)>=20000000:raise ValueError("source expanded beyond bounded cap")
 split=unpack.split(b"\n")
 if tuple(split[0].decode("utf-8-sig").rstrip("\r").split("\t"))!=EXPECTED:raise ValueError("header drift")
 flags=collections.Counter();spa=collections.Counter();maf=collections.Counter();info=collections.Counter()
 for line in split[1:-1]:
  if not line.strip():continue
  row=line.decode("utf-8").rstrip("\r").split("\t")
  if len(row)!=len(EXPECTED):raise ValueError("row malformed")
  flags[row[19]]+=1
  spa[row[14]]+=1
  try:
   v=float(row[7])
   if not math.isfinite(v):info["nonfinite"]+=1
   elif v<0:info["negative"]+=1
   elif v<0.3:info["lt_0_3"]+=1
   elif v<0.8:info["0_3_to_0_8"]+=1
   elif v<=1.001:info["ge_0_8"]+=1
   else:info["gt_1_001"]+=1
  except ValueError:info["nonnumeric"]+=1
  try:
   v=float(row[6])
   if v<=0 or v>=1:maf["monomorphic"]+=1
   elif min(v,1-v)<.001:maf["MAF_lt_0_001"]+=1
   elif min(v,1-v)<.01:maf["MAF_lt_0_01"]+=1
   else:maf["MAF_ge_0_01"]+=1
  except ValueError:maf["nonnumeric"]+=1
 return {"schema":"TPMI43321-SAIGE-QC-original-prefix-v1",
         "technical_status":"PASS_BOUNDED_SOURCE_QC_CONVENTION",
         "complete_prefix_rows":len(split)-2,
         "original_QC_labels":dict(flags),
         "Is_SPA_labels":dict(spa),
         "imputation_INFO_bins":dict(info),
         "AF_Allele2_bins":dict(maf),
         "full_source_QC_distribution_verified":False,
         "classify_nontrivial_QC_as_all_fail":False,
         "LDSC_MR_coloc_allowed":False,
         "scientific_status":"HOLD"}
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
 try:r=audit()
 except Exception as e:r={"schema":"TPMI43321-SAIGE-QC-original-prefix-v1","technical_status":"HOLD_ORIGINAL_SOURCE","error":str(e)[:250],"scientific_status":"HOLD"}
 a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n")
 print(json.dumps(r,ensure_ascii=False))
 if not r["technical_status"].startswith("PASS"):raise SystemExit(1)
if __name__=="__main__":main()
