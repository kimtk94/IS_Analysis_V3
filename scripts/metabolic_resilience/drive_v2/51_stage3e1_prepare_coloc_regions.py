#!/usr/bin/env python3
from pathlib import Path
import csv
import gzip
import zipfile
import json
import math

ROOT = Path("/srv/is-analysis")
WINDOWS = ROOT / "results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B2_GENE_CIS_WINDOWS.tsv"
REGISTRY = ROOT / "results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D0_OUTCOME_REGISTRY.tsv"
CIS = ROOT / "results/metabolic_resilience/stage3_full_pgwas/cis_marginal"
REFVCF = ROOT / "data/metabolic_resilience/stage2_gwas/ld_reference_1kg_eur/stage3_regional_vcf"
OUT = ROOT / "results/metabolic_resilience/stage3_coloc/inputs"
AUDIT = ROOT / "results/metabolic_resilience/stage3_full_pgwas/audit"
OUT.mkdir(parents=True, exist_ok=True)

# Load candidate windows.
windows=[]
with WINDOWS.open("r",encoding="utf-8",newline="") as f:
    for r in csv.DictReader(f,delimiter="\t"):
        if r["assembly"]=="GRCh37" and r["chrom"]:
            windows.append({
                "gene":r["gene_symbol"],
                "chrom":str(r["chrom"]).replace("chr",""),
                "start":int(r["cis_start"]),
                "end":int(r["cis_end"]),
            })

# Regional rsID -> GRCh37 location mapping from the exact EUR reference VCFs.
rsmap={}
for w in windows:
    p=REFVCF/f"{w['gene']}.1000G_EUR.b37.vcf.gz"
    if not p.exists():
        continue
    with gzip.open(p,"rt",encoding="utf-8",errors="replace") as f:
        for line in f:
            if line.startswith("#"): continue
            q=line.rstrip("\n").split("\t")
            if len(q)<5: continue
            chrom=q[0].replace("chr","")
            pos=int(q[1]); rsid=q[2]
            if rsid not in ("","."):
                for one in rsid.split(";"):
                    rsmap[one]=(chrom,pos)

# Exposure maps, exact position + allele set.
exposure={}
for p in CIS.glob("*_cis1Mb_marginal.tsv.gz"):
    with gzip.open(p,"rt",encoding="utf-8",newline="") as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    if not rows: continue
    gene=rows[0]["gene_symbol"]
    m={}
    for r in rows:
        key=(str(r["chrom_hg19"]),int(r["pos_hg19"]),frozenset([r["allele0"].upper(),r["allele1"].upper()]))
        m[key]=r
    exposure[gene]=m

registry=list(csv.DictReader(REGISTRY.open("r",encoding="utf-8"),delimiter="\t"))

def open_text(spec):
    path=Path(spec["path"])
    if path.suffix.lower()==".zip":
        z=zipfile.ZipFile(path)
        import io
        return io.TextIOWrapper(z.open(spec["member"]),encoding="utf-8",errors="replace"),z
    if str(path).endswith(".gz"):
        return gzip.open(path,"rt",encoding="utf-8",errors="replace"),None
    return open(path,"r",encoding="utf-8",errors="replace"),None

summary=[]

for spec in registry:
    trait=spec["trait"]
    if spec["status"]!="PASS":
        summary.append({"trait":trait,"gene":"*","matched":0,"status":"SKIP_REGISTRY"})
        continue

    # T2D input can be prepared, but case-control coloc execution will remain HOLD until s is known.
    capture={w["gene"]:[] for w in windows}

    fh,z=open_text(spec)
    header=None;H=None;tab=True
    try:
        for line in fh:
            s=line.strip()
            if not s or s.startswith("##"): continue
            if header is None:
                tab="\t" in line
                header=line.lstrip("#").rstrip("\n").split("\t") if tab else line.lstrip("#").split()
                H={x:i for i,x in enumerate(header)}
                continue
            p=line.rstrip("\n").split("\t") if tab else line.split()
            if len(p)<len(header): continue

            def val(c):
                return p[H[c]] if c and c in H and H[c]<len(p) else ""

            rsid=val(spec["col_rsid"])
            chrom="";pos=None

            # Prefer regional 1000G rsID mapping. This safely resolves source build differences.
            if rsid and rsid in rsmap:
                chrom,pos=rsmap[rsid]
            elif spec["build"]=="GRCh37" and spec["col_chr"] and spec["col_pos"]:
                chrom=val(spec["col_chr"]).replace("chr","")
                try: pos=int(float(val(spec["col_pos"])))
                except Exception: pos=None
            else:
                continue

            gene=None
            for w in windows:
                if chrom==w["chrom"] and pos is not None and w["start"]<=pos<=w["end"]:
                    gene=w["gene"];break
            if gene is None: continue

            ea=val(spec["col_ea"]).upper()
            oa=val(spec["col_oa"]).upper()
            if not ea or not oa: continue

            key=(chrom,pos,frozenset([ea,oa]))
            ex=exposure.get(gene,{}).get(key)
            if ex is None: continue

            try:
                beta=float(val(spec["col_beta"])) if spec["col_beta"] else math.log(float(val(spec["col_or"])))
                se=float(val(spec["col_se"]))
            except Exception:
                continue

            def fnum(c):
                if not c:return None
                try:return float(val(c))
                except Exception:return None

            capture[gene].append({
                "gene":gene,"trait":trait,"domain":spec["domain"],"outcome_type":spec["type"],
                "rsid":rsid,"chrom":chrom,"pos":pos,
                "exp_a0":ex["allele0"],"exp_a1":ex["allele1"],
                "exp_eaf":ex["eaf"],"beta_exp":ex["beta"],"se_exp":ex["se"],"n_exp":ex["n"],
                "out_ea":ea,"out_oa":oa,"out_eaf":fnum(spec["col_eaf"]),
                "beta_out":beta,"se_out":se,"n_out":fnum(spec["col_n"]),
            })
    finally:
        fh.close()
        if z is not None:z.close()

    for gene,rows in capture.items():
        dst=OUT/f"{gene}__{trait}.coloc.tsv.gz"
        fields=[
            "gene","trait","domain","outcome_type","rsid","chrom","pos",
            "exp_a0","exp_a1","exp_eaf","beta_exp","se_exp","n_exp",
            "out_ea","out_oa","out_eaf","beta_out","se_out","n_out"
        ]
        with gzip.open(dst,"wt",encoding="utf-8",newline="") as f:
            wr=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
            wr.writeheader();wr.writerows(rows)
        summary.append({
            "trait":trait,"gene":gene,"matched":len(rows),
            "status":"PASS" if len(rows)>=50 else "WARN_LOW_NSNP",
            "file":str(dst)
        })

sf=AUDIT/"STAGE3E1_COLOC_INPUT_SUMMARY.tsv"
with sf.open("w",encoding="utf-8",newline="") as f:
    fields=list(summary[0].keys()) if summary else ["status"]
    wr=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
    wr.writeheader();wr.writerows(summary)

print(json.dumps({
    "pairs":len(summary),
    "pass":sum(r["status"]=="PASS" for r in summary),
    "warn":sum(r["status"].startswith("WARN") for r in summary),
},indent=2))
