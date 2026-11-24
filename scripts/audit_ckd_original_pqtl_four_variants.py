#!/usr/bin/env python3
"""Targeted original UKB-PPP per-variant pQTL lookup, read-only archived gzip."""
import csv,gzip,io,json,tarfile
from pathlib import Path
BASE=Path("/srv/is-analysis/data/ckd/stage2_pqtl/EUR")
SPEC=[
 ("SDCCAG8","SDCCAG8_Q86SQ7_OID31423_v1_Oncology_II.tar","1",{"rs35646972":"242690143","rs6663222":"242990635"}),
 ("INHBC","INHBC_P55103_OID21093_v1_Neurology.tar","12",{"rs11834807":"58433974"}),
 ("UMOD","UMOD_P07911_OID20237_v1_Cardiometabolic.tar","16",{"rs11860566":"20208930"}),
]
results=[]
for gene,archive,chrom,targets in SPEC:
 path=BASE/gene/archive
 result={"gene":gene,"archive_exists":path.is_file(),"targets":list(targets),"records":[]}
 if not path.is_file():
  results.append(result);continue
 with tarfile.open(path,"r:") as tf:
  members=[m for m in tf.getmembers() if m.isfile() and ("discovery_chr"+chrom+"_") in m.name and m.name.endswith(".gz")]
  result["matching_members"]=[m.name for m in members]
  for m in members:
   handle=tf.extractfile(m)
   if handle is None:continue
   with gzip.GzipFile(fileobj=handle,mode="rb") as gz, io.TextIOWrapper(gz,encoding="utf-8",errors="replace") as txt:
    names=txt.readline().split()
    result.setdefault("diagnostic_columns",names[:20])
    result.setdefault("sample_format",[])
    sample_count=0
    result["nearest"]={k:{"distance":10**15} for k in targets}
    result["coordinate_diagnostic"]={"rows_checked":0,"genpos_idpos_different":0}
    for line in txt:
     values=line.split()
     if len(values)<len(names):continue
     row=dict(zip(names,values))
     if sample_count<2:
      result["sample_format"].append({"field_count":len(values),"first_chrom":row.get("CHROM"),"first_genpos":row.get("GENPOS"),"first_id_pattern":str(row.get("ID",""))[:12]})
      sample_count+=1
     pos=row.get("GENPOS","")
     rid=row.get("ID","")
     parts=rid.split(":")
     idpos=parts[1] if len(parts)>1 else ""
     if pos.isdigit() and idpos.isdigit():
      result["coordinate_diagnostic"]["rows_checked"]+=1
      result["coordinate_diagnostic"]["genpos_idpos_different"]+=int(pos!=idpos)
      for target,tp in targets.items():
       delta=min(abs(int(pos)-int(tp)),abs(int(idpos)-int(tp)))
       if delta<result["nearest"][target]["distance"]:
        result["nearest"][target]={"distance":delta,"GENPOS":pos,"ID":rid,"A1FREQ":row.get("A1FREQ"),"INFO":row.get("INFO")}
     if idpos in targets.values():
      pass
     if not (pos in targets.values() or rid in targets or idpos in targets.values()):continue
     result["records"].append({k:row.get(k) for k in ("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","BETA","SE","N","TEST")})
 results.append(result)
print(json.dumps({"audit":"ORIGINAL_UKBPPP_ARCHIVE_TARGETED_VARIANT_LOOKUP","results":results},indent=2))
