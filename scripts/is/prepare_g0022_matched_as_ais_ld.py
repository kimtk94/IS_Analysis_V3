#!/usr/bin/env python3
"""Match AS/AIS GWAS z scores to exactly same SNP order and signed EAS LD.

This isolates SNP-set effects in LD mismatch comparison. Related cohorts NOT
independent replicates, and source effective N is an approximation.
"""
import csv,gzip,json,argparse
from pathlib import Path
import numpy as np

ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/susie_g0022_pilot")
GWAS=ROOT.parent/"IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz"
AS="GCST90104544"
AIS="GCST90104545"
def shared(gwas):
    idx={}
    with gzip.open(gwas,"rt") as f:
        for row in csv.DictReader(f,delimiter="\t"):
            if row["group_id"]!="IS_XDATA_G0022" or row["dataset"]!=AIS:
                continue
            if row["qc_status"] not in ("MATCH_ALT_EFFECT","MATCH_REF_EFFECT"):
                continue
            variant=row["variant_id"]
            try:
                z=float(row["alt_effect_beta"])/float(row["se"])
                eaf=float(row["alt_effect_eaf"])
            except (TypeError,ValueError,ZeroDivisionError):
                continue
            if np.isfinite(z) and np.isfinite(eaf):
                idx[variant]=(z,eaf)
    return idx

def run(root,gwas):
    ais=shared(gwas)
    fpaths=sorted(root.glob("GCST90104544_AS_signal*/variants.tsv"))
    if len(fpaths)!=2:raise ValueError("Need two original AS LD-mismatch windows")
    dest=root/"matched_as_ais_ld_sensitivity"
    dest.mkdir(parents=True,exist_ok=True)
    records=[]
    for path in fpaths:
        directory=path.parent
        with path.open() as f:
            vals=list(csv.DictReader(f,delimiter="\t"))
        as_z=np.asarray([float(r["gwas_z"]) for r in vals])
        ld=np.loadtxt(directory/"signed_ld.tsv.gz",delimiter="\t",ndmin=2)
        if len(vals)!=ld.shape[0] or ld.shape[0]!=ld.shape[1]:
            raise ValueError("Original LD dimension mismatch")
        take=[];zs=[]
        for i,row in enumerate(vals):
            record=ais.get(row["variant_id"])
            if record is None:continue
            z,eaf=record
            ref=float(row["reference_alt_af"])
            if abs(eaf-ref)>0.15:continue
            take.append(i);zs.append(z)
        if len(take)<100:raise ValueError("Insufficient exact matched AS/AIS SNPs")
        newld=ld[np.ix_(take,take)]
        ordered=[vals[i]["variant_id"] for i in take]
        out=dest/directory.name
        out.mkdir(exist_ok=True)
        with (out/"variants.tsv").open("w") as f:
            f.write("variant_id\n"+"\n".join(ordered)+"\n")
        np.savetxt(out/"ld.tsv.gz",newld,delimiter="\t",fmt="%.17g")
        np.savetxt(out/"as_z.tsv",as_z[take],fmt="%.12g")
        np.savetxt(out/"ais_z.tsv",np.asarray(zs),fmt="%.12g")
        both=np.asarray(zs)
        r=float(np.corrcoef(as_z[take],both)[0,1])
        if not np.isfinite(r):raise ValueError("Undefined study Z-score correlation")
        rec=dict(window=directory.name,original_as_variants=len(vals),
          shared_study_variant_pairs=len(take),
          study_z_correlation=round(r,6),
          minimum_shared_variants=100,
          proof_level="EXACT_VARIANT_MATCHED_NOT_INDEPENDENT_REPLICATION")
        records.append(rec)
        print(json.dumps(rec))
    (dest/"MATCHED_AS_AIS_SHARED_VARIANTS_SUMMARY.json").write_text(json.dumps(records,indent=2))
    return records

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--gwas",type=Path,default=GWAS)
    a=p.parse_args()
    run(a.root,a.gwas)
