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
    reader=csv.DictReader(txt,delimiter=" ",skipinitialspace=True)
    for row in reader:
     pos=row.get("GENPOS","")
     rid=row.get("ID","")
     if not (pos in targets.values() or rid in targets):continue
     result["records"].append({k:row.get(k) for k in ("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","BETA","SE","N","TEST")})
 results.append(result)
print(json.dumps({"audit":"ORIGINAL_UKBPPP_ARCHIVE_TARGETED_VARIANT_LOOKUP","results":results},indent=2))
