#!/usr/bin/env python3
"""Read-only allele-aware VEP colocated rsID verification, GRCh37."""
import argparse,csv,json
from pathlib import Path

def alleles_from_colocated(x):
    vals=set()
    a=x.get("allele_string")
    if isinstance(a,str): vals.update(t.upper() for t in a.split("/") if t)
    for f in x.get("frequencies",{}).keys():
        vals.add(str(f).upper())
    return vals

def main():
    p=argparse.ArgumentParser()
    p.add_argument("annotation_tsv",type=Path)
    p.add_argument("raw_json",type=Path)
    p.add_argument("--out-dir",type=Path)
    a=p.parse_args()
    out=a.out_dir or a.annotation_tsv.parent
    with a.annotation_tsv.open(newline="") as f: rows=list(csv.DictReader(f,delimiter="\t"))
    raw=json.loads(a.raw_json.read_text())
    lookup={x["variant_id"]:x["response"] for x in raw}
    results=[]
    for r in rows:
        v=lookup.get(r["variant_id"])
        status="UNRESOLVED"; matched=[]
        if v:
            for x in v.get("colocated_variants",[]):
                rs=str(x.get("id",""))
                if not rs.startswith("rs"): continue
                pos=x.get("start")
                end=x.get("end")
                if pos!=int(r["pos"]) or end!=int(r["pos"])+len(r["ref"])-1: continue
                alleles=alleles_from_colocated(x)
                # VEP colocated allele_string can be absent: no false positive.
                if r["ref"].upper() in alleles and r["alt"].upper() in alleles:
                    matched.append(rs)
            if matched: status="VERIFIED_ALLELE_SET"
            elif r["vep_status"]!="API_ERROR": status="UNRESOLVED_ALLELE_EVIDENCE"
        results.append({**r,"allele_mapping_status":status,
          "allele_verified_rsids":";".join(sorted(set(matched))) if matched else "UNRESOLVED"})
    out.mkdir(parents=True,exist_ok=True)
    dest=out/"ALCOHOL_CS_VEP_RSID_ALLELE_AUDIT.tsv"
    fields=list(results[0])
    with dest.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t");w.writeheader();w.writerows(results)
    summary={"status":"DIAGNOSTIC_ONLY","total":len(results),
      "verified_allele_set":sum(x["allele_mapping_status"]=="VERIFIED_ALLELE_SET" for x in results),
      "unresolved":sum(x["allele_mapping_status"]!="VERIFIED_ALLELE_SET" for x in results),
      "note":"Allele set and coordinate match; not independent GWAS-build or functional validation."}
    (out/"ALCOHOL_CS_VEP_RSID_ALLELE_AUDIT.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("RSID_ALLELE_AUDIT_COMPLETE",json.dumps(summary))
if __name__=="__main__":main()
