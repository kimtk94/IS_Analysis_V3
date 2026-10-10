#!/usr/bin/env python3
"""GRCh37 VEP annotation for FASTA-verified exploratory CS SNPs.
Requires explicit --ref-verified acknowledgement; no causal promotion.
"""
import argparse,csv,json,time,urllib.request,urllib.error
from pathlib import Path

ENDPOINT="https://grch37.rest.ensembl.org/vep/human/region"
def fetch(chrom,pos,ref,alt,retries=3):
    # VEP region input: chr start end allele strand (forward)
    variant=f"{chrom} {pos} {pos+len(ref)-1} {ref}/{alt} +"
    body=json.dumps({"variants":[variant]}).encode()
    for attempt in range(retries):
        req=urllib.request.Request(ENDPOINT,data=body,
          headers={"Content-Type":"application/json","Accept":"application/json","User-Agent":"IS-Analysis-V3-CS-annotation/1.0"},method="POST")
        try:
            with urllib.request.urlopen(req,timeout=35) as resp:
                obj=json.load(resp)
            if not isinstance(obj,list) or len(obj)!=1: raise ValueError("Unexpected VEP response shape")
            return obj[0],None
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,ValueError) as e:
            error=f"{type(e).__name__}: {str(e)[:160]}"
            if attempt<retries-1: time.sleep(2**attempt)
    return None,error

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("input_tsv",type=Path)
    ap.add_argument("--ref-verified",action="store_true",help="11/11 REF matched GRCh37 FASTA")
    ap.add_argument("--out-dir",type=Path)
    args=ap.parse_args()
    if not args.ref_verified: ap.error("FASTA verification required: --ref-verified")
    out=args.out_dir or args.input_tsv.parent
    out.mkdir(parents=True,exist_ok=True)
    with args.input_tsv.open(newline="") as f: rows=list(csv.DictReader(f,delimiter="\t"))
    if not rows: raise ValueError("Empty input")
    for row in rows:
        if row["genome_build"]!="UNRESOLVED" or row["annotation_status"]!="BLOCKED_PENDING_BUILD":
            raise ValueError("Unexpected upstream state; review input before annotation")
    result=[]
    raw=[]
    for row in rows:
        v,err=fetch(row["chrom"],int(row["pos"]),row["ref"],row["alt"])
        if v is None:
            result.append({**row,"rsid":"UNRESOLVED","vep_consequence":"UNRESOLVED","vep_status":"API_ERROR","error":err})
            continue
        raw.append({"variant_id":row["variant_id"],"response":v})
        ids=[x for x in v.get("colocated_variants",[]) if isinstance(x,dict) and str(x.get("id","")).startswith("rs")]
        rsids=sorted(set(str(x["id"]) for x in ids))
        terms=v.get("most_severe_consequence")
        result.append({**row,"rsid":";".join(rsids) if rsids else "UNRESOLVED",
          "vep_consequence":terms or "UNRESOLVED","vep_status":"ANNOTATED_UNVERIFIED_ALLELE_MAPPING",
          "error":""})
    dest=out/"ALCOHOL_CS_VEP_GRCH37_EXPLORATORY.tsv"
    cols=list(result[0].keys())
    with dest.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter="\t",extrasaction="ignore");w.writeheader();w.writerows(result)
    (out/"ALCOHOL_CS_VEP_GRCH37_RAW.json").write_text(json.dumps(raw,indent=2)+"\n")
    manifest={"status":"EXPLORATORY_NO_CAUSAL_PROMOTION","reference":"GRCh37",
      "gwas_source_build":"UNRESOLVED","ref_check":"USER_REPORTED_11_OF_11_MATCH",
      "rsid_mapping":"COLOCATED_CANDIDATES_NOT_ALLELE_VERIFIED",
      "total":len(rows),"api_errors":sum(x["vep_status"]=="API_ERROR" for x in result),
      "output":str(dest)}
    (out/"ALCOHOL_CS_VEP_GRCH37_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print("VEP_GRCH37_ANNOTATION_COMPLETE",len(rows),"api_errors",manifest["api_errors"],dest)
if __name__=="__main__":main()
