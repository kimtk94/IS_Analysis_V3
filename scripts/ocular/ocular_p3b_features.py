#!/usr/bin/env python3
"""P3b: compare decompressed feature row identity/order across all 18 matrices."""
import argparse, gzip, hashlib, json, tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from ocular_p1_integrity import download, file_url

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument("--output",required=True)
    args=pa.parse_args()
    result={"run_utc":datetime.now(timezone.utc).isoformat(),
            "dataset":"GSE313736","comparison":"uncompressed 10x features.tsv.gz",
            "features":[], "status":"HOLD","errors":[]}
    baseline=None
    with tempfile.TemporaryDirectory(prefix="ocular-feature-") as temp:
        for i in range(9):
            pool=chr(ord("A")+i)
            for treatment,gsm,prefix in (
              ("vehicle",f"GSM{9373276+i}",f"pool_{pool}"),
              ("dex",f"GSM{9373285+i}",f"pool_dex{pool}")):
                row={"pool":pool,"condition":treatment,"gsm":gsm}
                try:
                    p=Path(temp)/(gsm+"features.tsv.gz")
                    row["transfer"]=download(file_url(gsm,prefix,"features.tsv.gz"),p,5*1024*1024)
                    with gzip.open(p,"rt",encoding="utf-8") as f:
                        values=[line.rstrip("\r\n").split("\t") for line in f]
                    row["n_features"]=len(values)
                    if not values or any(len(x)<2 for x in values):
                        raise ValueError("INVALID_FEATURE_COLUMNS")
                    ids=[x[0] for x in values]
                    row["distinct_gene_ids"]=len(set(ids))
                    row["uncompressed_gene_order_sha256"]=hashlib.sha256(
                        ("\n".join(ids)+"\n").encode()).hexdigest()
                    row["uncompressed_complete_row_sha256"]=hashlib.sha256(
                        ("\n".join("\t".join(x) for x in values)+"\n").encode()).hexdigest()
                    if len(set(ids))!=len(ids):raise ValueError("DUPLICATE_GENE_IDS")
                    if len(ids)!=36601:raise ValueError("UNEXPECTED_GENE_COUNT")
                    if baseline is None:
                        baseline=ids
                        result["reference_gsm"]=gsm
                    row["order_mismatches"]=sum(a!=b for a,b in zip(ids,baseline))
                    row["gene_id_set_symmetric_difference"]=len(set(ids)^set(baseline))
                    row["status"]="PASS" if not row["order_mismatches"] and not row["gene_id_set_symmetric_difference"] else "ORDER_DIFF"
                    if row["status"]!="PASS":result["errors"].append(gsm+"_FEATURE_ORDER_DIFFERENCE")
                except Exception as e:
                    row["status"]="ERROR"
                    row["reason"]=type(e).__name__+": "+str(e)[:150]
                    result["errors"].append(gsm+"_FEATURES_NOT_VERIFIED")
                result["features"].append(row)
    result["status"]="PASS_IDENTICAL_GENE_ORDER_ALL_18" if len(result["features"])==18 and not result["errors"] else "HOLD_FEATURE_ORDER"
    Path(args.output).write_text(json.dumps(result,indent=2)+"\n")
    print(result["status"],"verified",sum(x["status"]=="PASS" for x in result["features"]),"errors",result["errors"])
    if result["status"]!="PASS_IDENTICAL_GENE_ORDER_ALL_18":raise SystemExit(2)
if __name__=="__main__":main()
