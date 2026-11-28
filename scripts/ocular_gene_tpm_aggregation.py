#!/usr/bin/env python3
"""Convert verified kallisto transcript abundance to gene TPM for fixed rank scores."""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

TX_PIN="cb7305d9d6c0d3f336c388554a801513c4004d163cc7035ea448a3476fdc111e"

def require(ok,msg):
    if not ok: raise ValueError(msg)

def read_map(path):
    require(hashlib.sha256(path.read_bytes()).hexdigest()==TX_PIN,"tx2gene SHA256 mismatch")
    with path.open(newline="",encoding="utf8") as f:
        rd=csv.DictReader(f,delimiter="\t")
        require(rd.fieldnames==["transcript_id","gene_id_versioned","gene_id"],"mapping schema")
        mapping={}
        for row in rd:
            t=row["transcript_id"]
            require(t not in mapping,"duplicate mapping transcript")
            mapping[t]=(row["gene_id_versioned"],row["gene_id"])
    require(len(mapping)==252835 and len({g for _,g in mapping.values()})==62700,"mapping cardinality")
    return mapping

def aggregate(rows,mapping):
    seen=set()
    tpm=defaultdict(list);counts=defaultdict(list)
    for r in rows:
        # kallisto may return the original pipe-delimited FASTA header.
        parts=r["target_id"].split("|")
        t=parts[0]
        require(t in mapping and t not in seen,"unmapped or duplicate transcript")
        version,g=mapping[t]
        if len(parts)>1: require(parts[1]==version,"FASTA header gene identity mismatch")
        vals=[float(r[k]) for k in ("length","eff_length","est_counts","tpm")]
        require(all(math.isfinite(x) and x>=0 for x in vals),"invalid abundance")
        require(vals[0]>0,"invalid transcript length")
        seen.add(t);counts[g].append(vals[2]);tpm[g].append(vals[3])
    require(seen==set(mapping),"incomplete transcript universe")
    result={g:(math.fsum(tpm[g]),math.fsum(counts[g])) for g in tpm}
    total=math.fsum(v[0] for v in result.values())
    require(abs(total-1e6)<=100,"TPM sum outside rounding tolerance")
    return result

def selftest():
    m={"t1":("g1.1","g1"),"t2":("g1.1","g1"),"t3":("g2.1","g2")}
    def row(t,p,c=1):
        return dict(target_id=t,length="100",eff_length="90",est_counts=str(c),tpm=str(p))
    good=[row("t1|g1.1|x",200000,3),row("t2",300000,4),row("t3",500000,5)]
    assert aggregate(good,m)=={"g1":(500000.,7.),"g2":(500000.,5.)}
    bad=[good[:-1],good+[good[0]], [row("t1|g2.1|x",200000),*good[1:]], [row("t1",float("nan")),*good[1:]], [row("t1",0),*good[1:]], [row("unknown",200000),*good[1:]]]
    for rows in bad:
        try: aggregate(rows,m)
        except ValueError: pass
        else: raise AssertionError("expected rejection")
    print("GENE_TPM_SYNTHETIC_PASS: 7 cases; no real expression inputs")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--selftest",action="store_true")
    p.add_argument("--map",type=Path)
    p.add_argument("--abundance",type=Path)
    p.add_argument("--expected-abundance-sha")
    p.add_argument("--out",type=Path)
    a=p.parse_args()
    if a.selftest: selftest();return
    if not all((a.map,a.abundance,a.expected_abundance_sha,a.out)):p.error("all input paths and exact abundance SHA required")
    h=hashlib.sha256(a.abundance.read_bytes()).hexdigest()
    require(h==a.expected_abundance_sha,"abundance hash mismatch")
    m=read_map(a.map)
    with a.abundance.open(newline="",encoding="utf8") as f:
        rd=csv.DictReader(f,delimiter="\t")
        require(rd.fieldnames==["target_id","length","eff_length","est_counts","tpm"],"abundance schema")
        data=aggregate(rd,m)
    a.out.mkdir(parents=True,exist_ok=True)
    path=a.out/"gene_tpm.tsv"
    with path.open("w",newline="",encoding="utf8") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(["gene_id","tpm","estimated_counts"])
        for g,(t,c) in sorted(data.items()):w.writerow([g,format(t,".17g"),format(c,".17g")])
    result={"status":"GENE_TPM_TECHNICAL_AGGREGATION_ONLY","scientific_validation":False,"input_sha256":h,"tx2gene_sha256":TX_PIN,"gene_table_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"genes":len(data),"total_tpm":math.fsum(v[0] for v in data.values()),"total_estimated_counts":math.fsum(v[1] for v in data.values()),"nonzero_tpm_genes":sum(v[0]>0 for v in data.values()),"note":"Sum transcript TPM per gene, not a count-based differential expression test; no treatment contrast."}
    (a.out/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
if __name__=="__main__":main()
