#!/usr/bin/env python3
"""Sample authoritative per-SNP N_case and N_ctrl from TPMI SAIGE 433.21.

Reads only FIRST 1 MiB of official 816,992,833-byte SAIGE gzip through
HTTP byte Range. Therefore results must be described as a bounded prefix:
not full-source QC, not a complete inventory of all per-variant Ns, and
never PLINK2 sample counts. Does not publish effect values or variant IDs.
"""
from __future__ import annotations
import argparse,collections,json,math,urllib.request,zlib
from pathlib import Path

URL="https://pheweb.ibms.sinica.edu.tw/download/saige_gwas/433.21"
SOURCE_SIZE=816_992_833
PARTIAL_BYTES=1_048_576
EXPECTED=("CHR","POS","MarkerID","Allele1","Allele2","AC_Allele2",
          "AF_Allele2","imputationInfo","BETA","SE","Tstat","var",
          "p.value","p.value.NA","Is.SPA","AF_case","AF_ctrl",
          "N_case","N_ctrl","QC")
def sample():
    req=urllib.request.Request(URL,headers={"Range":f"bytes=0-{PARTIAL_BYTES-1}",
                                             "User-Agent":"IS-AnalysisV3-SAIGE-N_case-N_ctrl-source/1.0"})
    with urllib.request.urlopen(req,timeout=70) as r:
        if r.status!=206 or r.geturl()!=URL:
            raise ValueError("source HTTP byte-range not confirmed")
        if r.headers.get("Content-Range")!=f"bytes 0-{PARTIAL_BYTES-1}/{SOURCE_SIZE}":
            raise ValueError("source Content-Range/size drift")
        data=r.read(PARTIAL_BYTES+1)
    if len(data)!=PARTIAL_BYTES or data[:2]!=b"\x1f\x8b":
        raise ValueError("source incomplete or not gzip")
    zd=zlib.decompressobj(zlib.MAX_WBITS|16)
    plain=zd.decompress(data,20_000_000)
    if len(plain)>=20_000_000:raise ValueError("prefix expansion cap")
    parts=plain.split(b"\n")
    if not parts:raise ValueError("missing source header")
    header=parts[0].decode("utf-8-sig").strip().split("\t")
    if tuple(header)!=EXPECTED:
        raise ValueError("SAIGE source header or column order drift")
    # Last line may be incomplete since the gzip stream is partial.
    complete=parts[1:-1]
    count=collections.Counter();qc=collections.Counter()
    rows=0
    for line in complete:
        if not line.strip():continue
        row=line.decode("utf-8").rstrip("\r").split("\t")
        if len(row)!=len(EXPECTED):raise ValueError("partial SAIGE line schema drift")
        try:
            ncase=int(float(row[17]));ncontrol=int(float(row[18]))
        except (TypeError,ValueError,OverflowError):
            qc["N_case_or_N_ctrl_invalid"]+=1;continue
        if ncase<=0 or ncontrol<=0:
            qc["N_case_or_N_ctrl_nonpositive"]+=1;continue
        count[(ncase,ncontrol)]+=1;rows+=1
    if not rows:raise ValueError("no complete original SAIGE source rows")
    top=[{"N_case":a,"N_ctrl":b,"prefix_rows":n,"sum":a+b,
          "unmatched_to_phenotype_page":(a,b)!=(9249,304660),
          "not_the_same_as_PLINK2_OBS_CT":a+b!=204099}
         for (a,b),n in count.most_common(12)]
    return {
       "schema":"TPMI43321-SAIGE-per-variant-casecontrol-prefix-v1",
       "technical_status":"PASS_SOURCE_PARTIAL_SAIGE_N_CASE_N_CTRL",
       "official_source":URL,"official_source_size":SOURCE_SIZE,
       "range_compressed_bytes":len(data),
       "full_source_integrity_verified":False,
       "complete_prefix_variant_rows":len(complete),
       "valid_case_control_prefix_rows":rows,
       "distinct_sample_size_pairs_in_prefix":len(count),
       "top_case_control_pairs_in_prefix":top,
       "prefix_invalid_N_QC":dict(qc),
       "PLINK2_OBS_CT":204099,
       "PheWeb_cases":9249,"PheWeb_controls":304660,
       "prefix_N_not_proof_of_all_variant_N":True,
       "SAIGE_case_control_N_cannot_be_transferred_to_PLINK2":True,
       "LDSC_Neff_approved":False,
       "scientific_status":"HOLD"}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
    try:r=sample()
    except Exception as e:
        r={"schema":"TPMI43321-SAIGE-per-variant-casecontrol-prefix-v1",
           "technical_status":"HOLD_SOURCE_OR_CASE_CONTROL_PARSE",
           "error_type":type(e).__name__,"error":str(e)[:250],"scientific_status":"HOLD"}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(r,ensure_ascii=False))
if __name__=="__main__":main()
