#!/usr/bin/env python3
"""Model-specific TPMI SAIGE x GIGASTROKE EAS AIS header/assembly preflight.

Inspects official gzip prefixes, pinned source schema and known liftover
requirements. Does NOT merge SNPs or import PLINK2 model-specific N/effects.
Do not infer SAIGE's effect-allele orientation without primary validation.
"""
from __future__ import annotations
import argparse,json,urllib.request,urllib.parse,zlib
from pathlib import Path
from scripts.ckm_stroke.qc_tpmi_43321_saige_full_casecontrol import URL as SAIGE_URL,EXPECTED_HEADER as SAIGE_COLUMNS,SOURCE_BYTES as SAIGE_BYTES
from scripts.ckm_stroke.qc_gigastroke_eas_ais_original import URL as EAS_URL,SCHEMA as EAS_COLUMNS

def inspect(url,expected_host,max_size,expected_fields):
 req=urllib.request.Request(url,headers={"User-Agent":"IS_Analysis_V3-SAIGE-EAS-harmonization-source-preflight/1.0",
                                          "Range":"bytes=0-1048575"})
 with urllib.request.urlopen(req,timeout=90) as r:
  final=urllib.parse.urlsplit(r.geturl())
  if r.status!=206 or final.hostname!=expected_host or final.scheme!="https":
   raise ValueError("official host or range response changed")
  cr=r.headers.get("Content-Range","")
  if not cr.startswith("bytes 0-1048575/"):raise ValueError("invalid content range")
  declared=int(cr.split("/")[-1])
  if declared<1000000 or declared>max_size:raise ValueError("source size outside cap")
  b=r.read(1048577)
 if len(b)!=1048576 or b[:2]!=b"\x1f\x8b":raise ValueError("not complete official gzip prefix")
 decomp=zlib.decompressobj(zlib.MAX_WBITS|16)
 initial=decomp.decompress(b,8_000_000)
 if len(initial)>=8_000_000:raise ValueError("inflated preview exceeds cap")
 header=tuple(initial.split(b"\n",1)[0].decode("utf-8-sig").rstrip("\r").split("\t"))
 if header!=expected_fields:raise ValueError("official source schema unexpectedly changed")
 return {"source_url":url,"source_compressed_bytes":declared,
         "verified_header":list(header),
         "prefix_test_only":True,
         "full_gzip_hash_checked_this_run":False}
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
 out={"schema":"TPMI43321-SAIGE-EAS-AIS-method-independent-preflight-v1",
       "assume_SAIGE_equals_PLINK2":False,
       "reuse_PLINK2_1507261_SNP_intersection_as_SAIGE":False,
       "harmonized_pairs_computed":False,"analysis_allowed":False,"scientific_status":"HOLD"}
 try:
  s=inspect(SAIGE_URL,"pheweb.ibms.sinica.edu.tw",900000000,SAIGE_COLUMNS)
  if s["source_compressed_bytes"]!=SAIGE_BYTES:raise ValueError("SAIGE length drift")
  o=inspect(EAS_URL,"ftp.ebi.ac.uk",1500000000,EAS_COLUMNS)
  out.update({"technical_status":"PASS_BOTH_ORIGINAL_METHOD_SCHEMAS",
              "original_SAIGE":s,"original_EAS_AIS":o,
              "SAIGE_has_explicit_case_control_N":all(c in s["verified_header"] for c in ("N_case","N_ctrl")),
              "EAS_AIS_has_beta_se_p":all(c in o["verified_header"] for c in ("beta","standard_error","p_value")),
              "source_build_SAIGE_verified_in_this_run":False,
              "source_SAIGE_allele2_effect_orientation_verified_in_this_run":False,
              "full_SAIGE_N_snp_distribution_verified_in_this_run":False})
 except Exception as e:
  out.update({"technical_status":"HOLD_SOURCE_SCHEMA_OR_RANGE","error":str(e)[:300]})
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
 print(json.dumps({"status":out["technical_status"],"SAIGE_Ncase_Nctrl_schema":out.get("SAIGE_has_explicit_case_control_N"),
                   "harmonization_executed":False,"error":out.get("error"),
                   "scientific_status":"HOLD"}))
if __name__=="__main__":main()
