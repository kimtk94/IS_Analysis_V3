#!/usr/bin/env python3
"""Selective KCPS2 2025 Zenodo large-ZIP member by verified HTTP206 byte ranges.

13.18 GB archive is NOT downloaded. Reads official central directory from the
last 128 KiB and ZIP local header to locate members, selects only requested
ALCO_AMOUNT, SBP or DBP, downloads compressed member ranges in parallel,
streams raw-DEFLATE, checks original ZIP CRC32/uncompressed byte length, then
atomic-promotes inner .txt.gz. Full original ZIP MD5 is NOT checked, and audit
records that distinction. Default PLAN only. Public Zenodo file is CC-BY-4.0.
"""
import argparse,concurrent.futures,hashlib,json,struct,zlib
from pathlib import Path
import requests
from download_g0022_japanese_alcohol_ranges import fetch_range

ROOT=Path("/srv/is-analysis/data/is/gwas_exposure/kcps2_2025_selective")
META=ROOT/"ZENODO_15132424_METADATA.json"
TAIL=ROOT/"KCPS2_SAIGE_sumstats.zip.tail131072.bin"
FILES={"alcohol":"SAIGE_ALCO_AMOUNT_INFO.txt.gz",
       "sbp":"SAIGE_SBP_INFO.txt.gz",
       "dbp":"SAIGE_DBP_INFO.txt.gz"}
URL="https://zenodo.org/records/15132424/files/KCPS2_SAIGE_sumstats.zip?download=1"
CHUNK=4*1024*1024
def info(root):
    r=json.loads((root/"ZENODO_15132424_METADATA.json").read_text())
    if str(r["id"])!="15132424":raise ValueError("Wrong KCPS2 provenance")
    f=r["files"]
    if len(f)!=1 or f[0]["key"]!="KCPS2_SAIGE_sumstats.zip":
        raise ValueError("KCPS2 ZIP original source metadata mismatch")
    if f[0]["checksum"]!="md5:0c6e3f69ca73af97417ce8c4ea8d307c":
        raise ValueError("Unexpected original ZIP checksum")
    return int(f[0]["size"])
def central(root):
    s=(root/"KCPS2_SAIGE_sumstats.zip.tail131072.bin").read_bytes()
    if len(s)!=131072 or s.rfind(b"PK\x05\x06")<0:
        raise ValueError("Missing ZIP central directory/EOCD")
    entries={};i=s.find(b"PK\x01\x02")
    while i>=0 and i+46<=len(s) and s[i:i+4]==b"PK\x01\x02":
        v=struct.unpack_from("<4s6H3I5H2I",s,i)
        fnlen,extralen,comlen=v[10],v[11],v[12]
        name=s[i+46:i+46+fnlen].decode()
        extra=s[i+46+fnlen:i+46+fnlen+extralen]
        compressed,uncompressed,offset=v[8],v[9],v[16]
        e=0
        while e+4<=len(extra):
            tag,count=struct.unpack_from("<HH",extra,e);e+=4
            block=extra[e:e+count];e+=count
            if tag==1:
                k=0
                if uncompressed==0xffffffff:uncompressed=struct.unpack_from("<Q",block,k)[0];k+=8
                if compressed==0xffffffff:compressed=struct.unpack_from("<Q",block,k)[0];k+=8
                if offset==0xffffffff:offset=struct.unpack_from("<Q",block,k)[0];k+=8
        if name in entries:raise ValueError("Duplicate central-directory member")
        entries[name]={"name":name,"crc32":v[7],"compressed_bytes":compressed,
                       "uncompressed_bytes":uncompressed,"local_header_offset":offset,
                       "compression_method":v[4]}
        i+=46+fnlen+extralen+comlen
        if s[i:i+4]!=b"PK\x01\x02":break
    if len(entries)!=36:raise ValueError("KCPS2 expected 36 phenotype members")
    if not set(FILES.values()).issubset(entries):raise ValueError("Required KCPS2 traits absent")
    return entries
def get_exact(url,start,end,total):
    with requests.get(url,headers={"Range":f"bytes={start}-{end}"},timeout=(10,45)) as response:
        val=response.headers.get("Content-Range","")
        if response.status_code!=206 or val!=f"bytes {start}-{end}/{total}":
            raise RuntimeError("KCPS2 remote exact source-range response invalid")
        data=response.content
    if len(data)!=end-start+1:raise RuntimeError("KCPS2 source header truncated")
    return data
def locate_member(name,entries,total):
    d=entries[name];start=d["local_header_offset"]
    header=get_exact(URL,start,start+511,total)
    a=struct.unpack_from("<4s5H3I2H",header,0)
    if a[0]!=b"PK\x03\x04" or a[3]!=8:
        raise RuntimeError("ZIP member local header method/signature mismatch")
    fnlen,exlen=a[9],a[10]
    found=header[30:30+fnlen].decode()
    if found!=name or d["compression_method"]!=8:raise ValueError("ZIP member name/method disagreement")
    offset=start+30+fnlen+exlen
    if offset+d["compressed_bytes"]>total:raise ValueError("Outside original ZIP length")
    return offset
def one(job):
    source,start,end,size,path=job
    s=requests.Session()
    try:return start,end,fetch_range(s,URL,start,end,size,path,max_attempts=4)
    finally:s.close()
def fetch_trait(root,trait,workers=8,chunk_mib=4):
    total=info(root);entries=central(root);name=FILES[trait];item=entries[name]
    localname=root/name
    verified=root/(name+".member_crc_verified.json")
    if verified.is_file() and localname.is_file():
        v=json.loads(verified.read_text())
        h=hashlib.sha256()
        with localname.open("rb") as f:
            for block in iter(lambda:f.read(2**20),b""):h.update(block)
        if v.get("inner_gzip_sha256")==h.hexdigest() and v.get("ZIP_member_crc32")==item["crc32"]:
            print("REUSE_MEMBER_SHA256_AND_ZIP_CRC_VALID",name,flush=True);return v
    datastart=locate_member(name,entries,total)
    chunk=chunk_mib*2**20
    stage=root/(name+".zipdeflate.ranges");stage.mkdir(exist_ok=True,parents=True)
    jobs=[];parts=[]
    csize=item["compressed_bytes"]
    for begin in range(0,csize,chunk):
        end=min(csize,begin+chunk)
        startbyte=datastart+begin;lastbyte=datastart+end-1
        filepath=stage/f"{startbyte:012d}-{lastbyte:012d}.bin"
        parts.append((begin,end,filepath))
        if not (filepath.exists() and filepath.stat().st_size==end-begin):
            jobs.append((name,startbyte,lastbyte,total,filepath))
    print("KCPS2_MEMBER_START",trait,"zip_bytes",total,
          "member_compressed_bytes",csize,"uncompressed_member_bytes",item["uncompressed_bytes"],
          "chunks",len(parts),"remaining",len(jobs),flush=True)
    errors=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        pending=[pool.submit(one,job) for job in jobs]
        complete=0
        for f in concurrent.futures.as_completed(pending):
            try:
                f.result();complete+=1
                if complete==1 or complete%10==0 or complete==len(pending):
                    print("KCPS2_RANGE_DONE",trait,complete,len(pending),flush=True)
            except Exception as e:
                errors.append(str(e));print("KCPS2_RANGE_FAIL",str(e)[:300],flush=True)
    if errors:raise RuntimeError("KCPS2 incomplete source, not promoted "+repr(errors[:3]))
    temp=root/(name+".assembly.part")
    infl=zlib.decompressobj(-15);crc=0;produced=0;h=hashlib.sha256()
    with temp.open("wb") as out:
        for begin,end,path in parts:
            if path.stat().st_size!=end-begin:raise RuntimeError("KCPS2 source segment length mismatch")
            with path.open("rb") as f:
                for b in iter(lambda:f.read(1024*1024),b""):
                    dec=infl.decompress(b)
                    if dec:
                        crc=zlib.crc32(dec,crc);h.update(dec);produced+=len(dec);out.write(dec)
        dec=infl.flush()
        if dec:
            crc=zlib.crc32(dec,crc);h.update(dec);produced+=len(dec);out.write(dec)
    if not infl.eof or infl.unused_data:
        raise RuntimeError("KCPS2 central member deflate stream invalid")
    if crc!=item["crc32"] or produced!=item["uncompressed_bytes"]:
        raise RuntimeError("KCPS2 central CRC32 or uncompressed byte count mismatch")
    with temp.open("rb") as f:
        if f.read(2)!=b"\x1f\x8b":
            raise ValueError("Inner KCPS2 source is not gzip")
    temp.replace(localname)
    record={"trait":trait,"source":"zenodo_15132424_original_archive",
       "original_archive_size_bytes":total,
       "original_archive_md5_verified":False,
       "archive_member_name":name,
       "ZIP_member_crc32":crc,
       "ZIP_member_size_bytes":produced,
       "inner_gzip_sha256":h.hexdigest(),
       "source_zip_offset":datastart,
       "selected_member_only":True,
       "complete_member_retrieved":True,
       "status":"SELECTIVE_ZIP_MEMBER_CRC_VERIFIED_ARCHIVE_MD5_NOT_COMPUTED"}
    verified.write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2),flush=True)
    return record
if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--root",type=Path,default=ROOT)
    ap.add_argument("--trait",choices=(*FILES,"all"),default="alcohol")
    ap.add_argument("--workers",type=int,default=8);ap.add_argument("--chunk-mib",type=int,default=4)
    ap.add_argument("--execute",action="store_true")
    x=ap.parse_args()
    if x.workers<1 or x.workers>8 or x.chunk_mib<1 or x.chunk_mib>8:
        ap.error("1–8 workers, 1–8 MiB")
    targets=list(FILES) if x.trait=="all" else [x.trait]
    if not x.execute:
        vals=central(x.root);print(json.dumps({k:vals[FILES[k]] for k in targets},indent=2))
    else:
        for kind in targets:fetch_trait(x.root,kind,x.workers,x.chunk_mib)
