#!/usr/bin/env python3
"""Read-only source consistency check for original 646 ABF input/result rows."""
import argparse
import csv
import hashlib
import json
import statistics
from pathlib import Path

KEY=("locus","dataset_key","gene_base")
POSTERIOR=("PP.H0","PP.H1","PP.H2","PP.H3","PP.H4")


def load(file):
    with file.open(newline="",encoding="utf8") as h:
        x=list(csv.DictReader(h,delimiter="\t"))
    if len(x)!=646:
        raise ValueError("SOURCE_EXPECTED_646_ROWS:"+str(file))
    if any(not all(r.get(k) for k in KEY) for r in x):
        raise ValueError("MISSING_ABF_COMPOSITE_KEY")
    z={tuple(r[k] for k in KEY):r for r in x}
    if len(z)!=646:
        raise ValueError("NONUNIQUE_ABF_COMPOSITE_KEYS")
    return z


def run(index,master,annotated,out):
    if out.exists():
        raise FileExistsError("REFUSE_OVERWRITE")
    idx=load(index)
    raw=load(master)
    ann=load(annotated)
    if set(idx)!=set(raw) or set(raw)!=set(ann):
        raise ValueError("SOURCE_COMPOSITE_KEYS_MISMATCH")
    spreads=[]
    for key in sorted(idx):
        x,y,z=idx[key],raw[key],ann[key]
        if x["status"]!="READY" or y["status"]!="PASS" or z["status"]!="PASS":
            raise ValueError("SOURCE_STATUS_MISMATCH")
        if int(x["nsnps"])!=int(y["nsnps"]) or int(y["nsnps"])!=int(z["nsnps"]):
            raise ValueError("SOURCE_SNP_COUNT_MISMATCH")
        if abs(float(x["qtl_n_first"])-float(y["qtl_n"]))>1e-9 or abs(float(y["qtl_n"])-float(z["qtl_n"]))>1e-9:
            raise ValueError("SOURCE_QTL_N_MISMATCH")
        qmin,qmax=float(x["qtl_n_min"]),float(x["qtl_n_max"])
        if qmin<=0 or qmax<qmin or float(x["qtl_n_first"])<qmin or float(x["qtl_n_first"])>qmax:
            raise ValueError("SOURCE_QTL_N_RANGE_INVALID")
        for h in POSTERIOR:
            a,b=float(y[h]),float(z[h])
            if abs(a-b)>1e-10:
                raise ValueError("ANNOTATED_POSTERIOR_MODIFIED")
        spreads.append((qmax-qmin,key))
    deltas=[x[0] for x in spreads]
    def qtile(percentile):
        # Use inclusive linear interpolation, matching pandas default.
        v=sorted(deltas);h=(len(v)-1)*percentile/100
        a=int(h);f=h-a
        return v[a]*(1-f)+v[min(a+1,len(v)-1)]*f
    return dict(
        schema="IS_ABF_SOURCE_646_INDEX_AUDIT_V1",
        status="INDEX_MASTER_ANNOTATED_CONSISTENT_ALL_646",
        rows=646,
        unique_loci=len({k[0] for k in idx}),
        unique_datasets=len({k[1] for k in idx}),
        unique_gene_ids=len({k[2] for k in idx}),
        indexed_ready=646,
        raw_pass=646,
        annotated_pass=646,
        matching_snp_count=646,
        matching_first_scalar_QTL_N=646,
        matching_H0_to_H4_annotation=646,
        pairs_with_variable_snp_specific_QTL_N=sum(d>0 for d in deltas),
        N_range_median=statistics.median(deltas),
        N_range_q25=qtile(25),
        N_range_q75=qtile(75),
        N_range_max=max(deltas),
        source_sha256={
            "index":hashlib.sha256(index.read_bytes()).hexdigest(),
            "raw_master":hashlib.sha256(master.read_bytes()).hexdigest(),
            "annotated_master":hashlib.sha256(annotated.read_bytes()).hexdigest()
        },
        source_paths={"index":str(index),"raw_master":str(master),"annotated_master":str(annotated)},
        claim_gate="PROVENANCE_ONLY_NOT_INDEPENDENT_646_SNP_REPLAY"
    )


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for field in ("index","master","annotated","out"):
        parser.add_argument("--"+field,type=Path,required=True)
    a=parser.parse_args()
    audit=run(a.index,a.master,a.annotated,a.out)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(audit,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:audit[k] for k in ("status","rows","pairs_with_variable_snp_specific_QTL_N","N_range_median","N_range_q25","N_range_q75","N_range_max")},indent=2))


if __name__=="__main__":
    main()
