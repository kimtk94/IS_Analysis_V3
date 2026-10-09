#!/usr/bin/env python3
"""Prepare isolated G0022 AIS/AS LD+z matrices for *exploratory* SuSiE-RSS.

No causal claims: 1KG-EAS reference is small (504), effective case-control N
unknown, GIGASTROKE AS/AIS overlap, windows truncated at +/-250kb.
"""
import argparse,csv,gzip,json,subprocess
from pathlib import Path
from collections import defaultdict
import numpy as np
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2")
REF=Path("/srv/is-analysis/data/is/ld_reference/broad_v1/1kg_eas/IS_XDATA_G0022.stable")
CLUMPS=ROOT/"clump_screen_v1"
OUT=ROOT/"susie_g0022_pilot"
DATASETS={"GCST90104545":"AIS","GCST90104544":"AS"}
def read_clumps(path):
    with path.open() as f:
        return list(csv.DictReader(f,delimiter="\t"))
def source_variants(path):
    g=defaultdict(dict)
    with gzip.open(path,"rt") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            if (r["group_id"]!="IS_XDATA_G0022" or
                r["dataset"] not in DATASETS or
                r["qc_status"] not in ("MATCH_REF_EFFECT","MATCH_ALT_EFFECT")):continue
            if {r["ref"],r["alt"]} in ({"A","T"},{"C","G"}):continue
            if r["alt_effect_beta"] in ("","NA"):continue
            key=(r["dataset"],r["variant_id"])
            try:z=float(r["alt_effect_beta"])/float(r["se"])
            except (ValueError,TypeError,ZeroDivisionError):continue
            if not np.isfinite(z):continue
            g[r["dataset"]][r["variant_id"]]=(z,r)
    return g
def traw_rows(path):
    with path.open() as f:
        head=next(f).strip().split("\t")
        if head[:6]!=["CHR","SNP","(C)M","POS","COUNTED","ALT"]:
            raise ValueError("Unexpected PLINK dosage export schema")
        n=len(head)-6
        for l in f:
            a=l.rstrip("\n").split("\t",6)
            if len(a)!=7:continue
            v=np.fromstring(a[6].replace("NA","nan"),sep="\t")
            if len(v)!=n:raise ValueError("Incomplete dosage row")
            yield a[:6],v
def dosage_alt_orientation(dose, counted, other, ref, alt):
    if counted==ref and other==alt:
        return 2-dose
    if counted==alt and other==ref:
        return dose.copy()
    raise ValueError("Dosage counted allele is incompatible with GWAS REF/ALT")
def build(args):
    if args.radius_bp<50000 or args.radius_bp>500000:
        raise ValueError("Pilot radius must be 50-500kb")
    src=source_variants(args.root/"IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz")
    out=args.out;out.mkdir(parents=True,exist_ok=True)
    reports=[]
    for ds,trait in DATASETS.items():
        cl=read_clumps(args.root/"clump_screen_v1"/(ds+"__IS_XDATA_G0022.clumps"))
        for i,entry in enumerate(cl,1):
            center=int(entry["POS"])
            if entry["ID"] not in src[ds]:
                raise ValueError("GWS clump index absent from allele-oriented GWAS "+entry["ID"])
            name=f"{ds}_{trait}_signal{i}_{center}"
            folder=out/name;folder.mkdir(exist_ok=True)
            prefix=folder/"reference_1kg_eas"
            traw=Path(str(prefix)+".traw")
            if not traw.exists():
                cmd=["plink2","--pfile",str(args.reference),
                    "--chr","12","--from-bp",str(center-args.radius_bp),
                    "--to-bp",str(center+args.radius_bp),
                    "--export","Av","--out",str(prefix)]
                p=subprocess.run(cmd,capture_output=True,text=True,timeout=180)
                if p.returncode:
                    raise RuntimeError("PLINK dosage export failed: "+p.stderr+p.stdout[-900:])
            dosage=[];zs=[];idx=[];stats=defaultdict(int)
            for fields,vec in traw_rows(traw):
                chrom,vid,_,pos,counted,other=fields
                match=src[ds].get(vid)
                if match is None:continue
                z,r=match
                ref,alt=r["ref"],r["alt"]
                try:
                    alt_dose=dosage_alt_orientation(vec,counted,other,ref,alt)
                except ValueError:
                    stats["dosage_orientation_mismatch"]+=1
                    continue
                missing=~np.isfinite(alt_dose)
                nonmiss=alt_dose[~missing]
                if len(nonmiss)<0.95*len(alt_dose):
                    stats["high_missing"]+=1;continue
                if len(nonmiss)<100:continue
                freq=float(nonmiss.mean()/2)
                if not (0.01<=freq<=0.99):stats["low_reference_maf"]+=1;continue
                gwas_freq=float(r["alt_effect_eaf"]) if r["alt_effect_eaf"] else None
                if gwas_freq is not None and abs(gwas_freq-freq)>0.15:
                    stats["frequency_discordance"]+=1;continue
                if missing.any():alt_dose[missing]=nonmiss.mean()
                if float(alt_dose.std())<1e-7:continue
                dosage.append(alt_dose)
                zs.append(z)
                idx.append({"variant_id":vid,"pos":pos,"ref":ref,"alt":alt,
                    "gwas_z":z,"gwas_p":r["p"],
                    "reference_alt_af":freq,
                    "gwas_alt_eaf":gwas_freq if gwas_freq is not None else ""})
            if not dosage:raise ValueError("No matching dosages in "+name)
            mat=np.asarray(dosage,dtype=np.float64)
            arr=(mat-mat.mean(axis=1,keepdims=True))/mat.std(axis=1,ddof=1,keepdims=True)
            ld=arr@arr.T/(arr.shape[1]-1)
            ld=(ld+ld.T)*0.5
            if not np.isfinite(ld).all():raise ValueError("NaN LD matrix")
            if np.max(np.abs(np.diag(ld)-1))>1e-6:raise ValueError("Bad LD diagonal")
            if entry["ID"] not in {x["variant_id"] for x in idx}:
                stats["lead_failed_qc"]+=1
                report={"study":ds,"trait":trait,"signal":name,"status":"BLOCKED_LEAD_QC",
                    "variant_count":len(idx),"reason":"Clump index absent after allele/MAF/FREQ QC"}
                reports.append(report)
                continue
            # Write reproducible artifacts. Rows/columns are index-order identical.
            with (folder/"variants.tsv").open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=list(idx[0]),delimiter="\t");w.writeheader();w.writerows(idx)
            # Keep full double precision: roundoff can introduce negative
            # eigenvalues in rank-deficient 504-person reference LD.
            np.savetxt(folder/"signed_ld.tsv.gz",ld,delimiter="\t",fmt="%.17g")
            with (folder/"z.tsv").open("w") as f:
                f.write("z\n")
                for z in zs:f.write(format(z,".12g")+"\n")
            report={"study":ds,"trait":trait,"signal":name,"clump_index":entry["ID"],
                "center":center,"radius_bp":args.radius_bp,"status":"PILOT_INPUT_QC_PASS",
                "variant_count":len(idx),"reference_samples":mat.shape[1],
                "mean_gwas_abs_z":round(float(np.mean(np.abs(zs))),4),
                "max_gwas_abs_z":round(float(np.max(np.abs(zs))),4),
                "dosage_allele_orientation":"ALT_SIGNED","matrix_diagonal_max_deviation":float(np.max(np.abs(np.diag(ld)-1))),
                "filters":dict(stats),"limitation":"Input only; 250kb truncation, uncertain case/control effective N, no biological causal inference"}
            (folder/"input_qc.json").write_text(json.dumps(report,indent=2))
            reports.append(report)
            print(json.dumps(report),flush=True)
    (out/"G0022_SUSIE_PILOT_INPUT_MANIFEST.json").write_text(json.dumps(reports,indent=2))
    return reports
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--reference",type=Path,default=REF)
    p.add_argument("--out",type=Path,default=OUT)
    p.add_argument("--radius-bp",type=int,default=250000)
    build(p.parse_args())
