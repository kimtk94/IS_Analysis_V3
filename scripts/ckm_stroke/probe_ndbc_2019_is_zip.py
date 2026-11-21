#!/usr/bin/env python3
"""Inspect official NBDC NHA000096 ZIP central directory via HTTP range requests.
Do not download the 2.3 GB archive to identify its members. Fail closed on bad ranges.
"""
import argparse,json,re,struct,urllib.request
from pathlib import Path

URL="https://humandbs-production.ddbj.nig.ac.jp/files/hum0014/hum0014.v17.CI.v1.zip"
def request(method="GET",start=None,end=None):
 headers={"User-Agent":"IS_Analysis_V3-NBDC-BBJ-2019-ZIP-probe/1.0"}
 if start is not None:headers["Range"]=f"bytes={start}-{end}"
 req=urllib.request.Request(URL,method=method,headers=headers)
 with urllib.request.urlopen(req,timeout=80) as r:
  if r.geturl()!=URL: raise RuntimeError("Unapproved redirect")
  data=r.read((end-start+1)+1 if start is not None else 0)
  if start is not None and (r.status!=206 or len(data)!=end-start+1):raise RuntimeError("Range request not honored")
  return {"http":r.status,"content_length":r.headers.get("Content-Length"),
          "content_range":r.headers.get("Content-Range"),
          "accept_ranges":r.headers.get("Accept-Ranges")},data
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
 out={"schema":"ndbc-BBJ-CI-original-zip-members-v1","url":URL,"dataset":"NHA000096",
      "target":"GRCh37:19:15573379:T:C","scientific_status":"HOLD"}
 try:
  h,_=request("HEAD");size=int(h["content_length"])
  if not 1000000<size<8000000000:raise ValueError("Unexpected source size")
  start=max(0,size-262144);r,blob=request(start=start,end=size-1)
  signature=b"PK\\x05\\x06"
  # Construct byte sequence from escaped text to keep script readable
  signature=bytes([80,75,5,6])
  ix=blob.rfind(signature)
  if ix<0:raise ValueError("ZIP end-of-central-directory not found in last 256 KB")
  if len(blob)<ix+22:raise ValueError("truncated ZIP EOCD")
  disk,cd_disk,n_disk,total,cd_size,cd_offset,comment=struct.unpack_from("<4H2LH",blob,ix+4)
  if disk or cd_disk or n_disk!=total or 0xffffffff in (cd_size,cd_offset):
   raise ValueError("ZIP64/spanned archive requires separate verified parser")
  if cd_size>4000000:raise ValueError("central directory exceeds 4 MB cap")
  cdr,cd=request(start=cd_offset,end=cd_offset+cd_size-1)
  items=[];cur=0
  for _ in range(total):
   if cd[cur:cur+4]!=bytes([80,75,1,2]):raise ValueError("invalid central directory entry")
   vals=struct.unpack_from("<4s6H3L5H2L",cd,cur)
   # compressed bytes, original bytes, name length, extra, comment
   compressed,uncompressed,name_len,extra_len,comment_len=vals[8],vals[9],vals[10],vals[11],vals[12]
   name=cd[cur+46:cur+46+name_len].decode("utf-8","replace")
   items.append({"filename":name[:200],"compressed_bytes":compressed,"uncompressed_bytes":uncompressed})
   cur+=46+name_len+extra_len+comment_len
  if cur!=len(cd):raise ValueError("central directory trailing bytes")
  out.update({"technical_status":"PASS_OFFICIAL_ZIP_CENTRAL_DIRECTORY",
             "zip_bytes":size,"members":items,"zip_member_count":len(items),
             "range_requests_verified":True,
             "source_association_observed":False})
 except Exception as exc:
  out.update({"technical_status":"HOLD_SOURCE_ZIP_STRUCTURE",
              "error_type":type(exc).__name__,"error":str(exc)[:240]})
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
 print(json.dumps(out,ensure_ascii=False))
if __name__=="__main__":main()
