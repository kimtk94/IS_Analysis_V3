#!/usr/bin/env python3
"""Verify TPMI SAIGE source CHR/POS/MarkerID/Allele1/Allele2 consistency.

Authentic HTTP range gzip PREFIX ONLY (1 MiB). No chromosome-build inference,
GIGASTROKE outcome effects, or genomewide harmonization claim.
"""
from __future__ import annotations
import argparse,collections,json,re,urllib.request,zlib
from pathlib import Path
from scripts.ckm_stroke.sample_tpmi_43321_saige_casecontrol import (
 URL, SOURCE_SIZE,PARTIAL_BYTES,EXPECTED
)
def check_row(row):
    if len(row)!=len(EXPECTED):return "BAD_WIDTH"
    chrom,pos,marker,al1,al2=row[:5]
    if not pos.isdecimal() or not chrom.isdecimal():return "BAD_CHROM_POS"
    if len(al1)!=1 or len(al2)!=1 or al1==al2 or any(c not in "ACGT" for c in al1+al2):
        return "NON_SNVS"
    # Marker strings can use ':' or '_' and may encode alleles/reference on
    # either strand. First require chromosome and position, then optionally
    # count whether marker actually specifies allele strings.
    cleaned=marker.removeprefix("chr").replace("_",":")
    parts=cleaned.split(":")
    if len(parts)<2 or parts[0]!=chrom or parts[1]!=pos:
        return "MARKER_COORD_DISCORDANT"
    if len(parts)>=4 and all(len(a)==1 and a in "ACGT" for a in parts[2:4]):
        if {parts[2],parts[3]}!={al1,al2}:
            return "MARKER_ALLELES_DISCORDANT"
        return "PASS_COORD_AND_ALLELE_SET"
    return "PASS_COORD_NO_MARKER_ALLELE"
def run():
    req=urllib.request.Request(URL,headers={
       "Range":f"bytes=0-{PARTIAL_BYTES-1}",
       "User-Agent":"IS-Analysis-V3-TPMI-SAIGE-Marker-ID-convention/1.0"})
    with urllib.request.urlopen(req,timeout=75) as r:
        if r.status!=206 or r.geturl()!=URL:raise ValueError("PheWeb range unavailable")
        if r.headers.get("Content-Range")!=f"bytes 0-{PARTIAL_BYTES-1}/{SOURCE_SIZE}":
            raise ValueError("source size/version mismatch")
        b=r.read(PARTIAL_BYTES+1)
    if len(b)!=PARTIAL_BYTES:raise ValueError("partial file incomplete")
    z=zlib.decompressobj(zlib.MAX_WBITS|16).decompress(b,20_000_000)
    if len(z)>=20_000_000:raise ValueError("decompressed source prefix exceeds cap")
    lines=z.split(b"\n")
    if tuple(lines[0].decode("utf-8-sig").rstrip("\r").split("\t"))!=EXPECTED:
        raise ValueError("SAIGE official header altered")
    results=collections.Counter()
    first_coord=None
    for raw in lines[1:-1]:
        if not raw.strip():continue
        row=raw.decode("utf-8").rstrip("\r").split("\t")
        results[check_row(row)]+=1
        if first_coord is None:first_coord={"chromosome":row[0],"position":row[1]}
    return {
      "schema":"TPMI43321-SAIGE-MarkerID-prefix-coordinate-qc-v1",
      "technical_status":"PASS_PREFIX_MARKER_CHECK" if
        results["MARKER_COORD_DISCORDANT"]==results["MARKER_ALLELES_DISCORDANT"]==results["BAD_WIDTH"]==0 else
        "HOLD_MARKER_ORIENTATION_INCONSISTENT",
      "official_source_url":URL,
      "source_prefix_bytes":PARTIAL_BYTES,
      "complete_prefix_rows":sum(results.values()),
      "MarkerID_CHR_POS_Allele_consistency":dict(results),
      "first_source_chrom_pos":first_coord,
      "genome_build_GRCh38_independently_verified":False,
      "allele2_effect_convention_source_method_expected":True,
      "allele2_effect_convention_per_variant_validated":False,
      "full_source_verified":False,"full_SAIGE_EAS_AIS_harmonization_executed":False,
      "scientific_status":"HOLD"}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);args=p.parse_args()
    try:result=run()
    except Exception as exc:result={"schema":"TPMI43321-SAIGE-MarkerID-prefix-coordinate-qc-v1",
     "technical_status":"HOLD_ORIGINAL_SOURCE_OR_PARSER",
     "error_type":type(exc).__name__,"error":str(exc)[:300],
     "scientific_status":"HOLD"}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(result,ensure_ascii=False))
    if not result["technical_status"].startswith("PASS"):raise SystemExit(1)
if __name__=="__main__":main()
