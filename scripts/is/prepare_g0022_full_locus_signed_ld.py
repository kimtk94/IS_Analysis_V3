#!/usr/bin/env python3
"""Full 4.05-Mb G0022 EAS AIS reference LD input, explicit --execute only.

Read-only canonical GWAS and reference PGEN. Creates a row-major binary signed
SNP-correlation matrix with blockwise bounded RAM. NOT fine-mapping validity.
"""
import argparse,csv,gzip,hashlib,json,math,os,subprocess
from pathlib import Path
from collections import Counter
import numpy as np

ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2")
REFERENCE=Path("/srv/is-analysis/data/is/ld_reference/broad_v1/1kg_eas/IS_XDATA_G0022.stable")
OUT=ROOT/"g0022_full_locus_ais_v1"
LEADS=["12:112241766:G:A","12:113031474:G:A","12:110675363:C:T"]
GWAS="GCST90104545"
def source(path):
    selected={}
    with gzip.open(path,"rt") as f:
        for x in csv.DictReader(f,delimiter="\t"):
            if x["group_id"]!="IS_XDATA_G0022" or x["dataset"]!=GWAS:
                continue
            if x["qc_status"] not in ("MATCH_ALT_EFFECT","MATCH_REF_EFFECT"):
                continue
            key=x["variant_id"]
            if key in selected:
                raise ValueError("Duplicate variant ID in GWAS "+key)
            if {x["ref"],x["alt"]} in ({"A","T"},{"C","G"}):continue
            z=float(x["alt_effect_beta"])/float(x["se"])
            if not math.isfinite(z):continue
            selected[key]=(z,x)
    if len(selected)<3000:raise ValueError("Insufficient genome-wide matched G0022 AIS variants")
    if not set(LEADS).issubset(set(selected)):
        raise ValueError("An AIS clump lead is unavailable from GWAS candidate set")
    return selected
def mem_available_bytes():
    with open("/proc/meminfo") as f:
        for l in f:
            if l.startswith("MemAvailable:"):
                return int(l.split()[1])*1024
    raise RuntimeError("Cannot inspect available RAM")
def build(args):
    chosen=source(args.root/"IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz")
    if len(chosen)>args.max_variants:
        raise RuntimeError(f"{len(chosen)} variants exceed budget {args.max_variants}")
    avail=mem_available_bytes()
    plan={"source":GWAS,"group_id":"IS_XDATA_G0022",
        "source_harmonized_variants":len(chosen),
        "LD_reference":"1000Genomes_EAS_504","expected_memory_available_GiB":round(avail/2**30,2),
        "matrix_byte_estimate":8*len(chosen)**2,
        "status":"PLAN_ONLY","claim":"NO_VALIDATED_FINE_MAPPING"}
    if not args.execute:
        print(json.dumps(plan,indent=2));return plan
    if avail<3*2**30:raise RuntimeError("Available RAM under 3GiB: refuse full-locus LD")
    args.out.mkdir(parents=True,exist_ok=True)
    extract=args.out/"requested_allele_matched_ids.txt"
    extract.write_text("\n".join(sorted(chosen))+"\n")
    traw=args.out/"reference_export.traw"
    if not traw.is_file():
        cmd=["plink2","--pfile",str(args.reference),"--extract",str(extract),
            "--export","Av","--out",str(args.out/"reference_export")]
        r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,timeout=300)
        if r.returncode:
            raise RuntimeError("PLINK export failed "+r.stdout[-700:]+r.stderr[-700:])
    allrows=[];X=[];qc=Counter()
    with traw.open() as f:
        h=next(f).rstrip("\n").split("\t")
        if h[:6]!=["CHR","SNP","(C)M","POS","COUNTED","ALT"]:
            raise ValueError("Unexpected PLINK dosage schema")
        n=len(h)-6
        if n!=504:raise ValueError("Expected exactly 504 reference samples")
        for line in f:
            a=line.rstrip("\n").split("\t",6)
            if len(a)<7:raise ValueError("Truncated genotype record")
            vid=a[1];match=chosen.get(vid)
            if match is None:raise ValueError("Unlisted variant exported: "+vid)
            z,gw=match
            geno=np.fromstring(a[6].replace("NA","nan"),sep="\t")
            if len(geno)!=n:raise ValueError("Wrong number of dosages")
            if a[4]==gw["ref"] and a[5]==gw["alt"]:
                alt=2-geno
            elif a[4]==gw["alt"] and a[5]==gw["ref"]:
                alt=geno
            else:raise RuntimeError("Conflicting counted allele "+vid)
            mask=np.isfinite(alt)
            if int(mask.sum())<int(math.ceil(.95*n)):
                qc["too_many_missing"]+=1;continue
            nonmiss=alt[mask]
            freq=float(nonmiss.mean()/2)
            if not .01<=freq<=.99:qc["maf_filter"]+=1;continue
            eaf=float(gw["alt_effect_eaf"]) if gw["alt_effect_eaf"] else None
            if eaf is not None and abs(eaf-freq)>.15:
                qc["eaf_discordant"]+=1;continue
            if not mask.all():alt[~mask]=nonmiss.mean()
            if np.std(alt,ddof=1)<1e-7:
                qc["monomorphic"]+=1;continue
            X.append(alt)
            allrows.append({"variant_id":vid,"chr":gw["chr"],"pos":int(gw["pos"]),
                "ref":gw["ref"],"alt":gw["alt"],"z":z,"p":gw["p"],
                "reference_alt_af":freq,"gwas_alt_eaf":gw["alt_effect_eaf"],
                "status":"ALT_SIGNED_CONFIRMED"})
    if len(allrows)<2500:raise ValueError("Too few variants after real genotype QC")
    if not set(LEADS).issubset({r["variant_id"] for r in allrows}):
        raise ValueError("A clump index was removed by genotype QC")
    order=sorted(range(len(allrows)),key=lambda i:(int(allrows[i]["pos"]),allrows[i]["variant_id"]))
    X=np.asarray([X[i] for i in order],dtype=np.float64)
    allrows=[allrows[i] for i in order]
    X-=X.mean(axis=1,keepdims=True)
    X/=X.std(axis=1,ddof=1,keepdims=True)
    p=X.shape[0]
    temp=args.out/"ld.rowmajor.f64.part"
    if temp.exists():raise FileExistsError("Partial LD matrix already exists; examine first")
    md=np.memmap(temp,mode="w+",dtype="<f8",shape=(p,p))
    block=256
    for i in range(0,p,block):
        md[i:i+block,:]=X[i:i+block,:]@X.T/(n-1)
    md.flush()
    diff=0.
    for i in range(0,p,block):
        diff=max(diff,float(np.max(np.abs(md[i:i+block,:]-md[:,i:i+block].T))))
    if diff>1e-7:raise RuntimeError("LD symmetry failed")
    if np.max(np.abs(np.diag(md)-1))>1e-6:raise RuntimeError("LD diagonal failed")
    del md
    ld=args.out/"ld.rowmajor.f64"
    temp.replace(ld)
    with (args.out/"variants.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter="\t")
        w.writeheader();w.writerows(allrows)
    report={**plan,"status":"FULL_G0022_SIGNED_LD_INPUT_QC_PASS",
       "genotype_snp_count":p,"genotype_n":n,
       "matrix_dimension":p,"LD_matrix_bytes":ld.stat().st_size,
       "matrix_max_asymmetry":diff,"excluded":dict(qc),
       "lead_ids_present":LEADS,
       "LD_matrix_sha256":hashlib.sha256(ld.read_bytes()).hexdigest(),
       "science":"FULL_LOCUS_INPUT_ONLY_NO_FINE_MAPPING_OR_CAUSAL_GENE_VALIDATION"}
    (args.out/"FULL_LOCUS_INPUT_QC.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    return report
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,default=ROOT)
    ap.add_argument("--reference",type=Path,default=REFERENCE)
    ap.add_argument("--out",type=Path,default=OUT)
    ap.add_argument("--max-variants",type=int,default=6500)
    ap.add_argument("--execute",action="store_true")
    build(ap.parse_args())
