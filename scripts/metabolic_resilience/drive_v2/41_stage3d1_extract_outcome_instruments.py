#!/usr/bin/env python3
from pathlib import Path
import csv,gzip,zipfile,math,json,os,sys,io

ROOT=Path(os.environ.get("IS_ANALYSIS_ROOT","/srv/is-analysis"))
INST=ROOT/"results/metabolic_resilience/stage3_full_pgwas/instruments/STAGE3C4_FROZEN_INSTRUMENTS_ANNOTATED.tsv"
REG=ROOT/"results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D0_OUTCOME_REGISTRY.tsv"
OUT=ROOT/"results/metabolic_resilience/stage3_confirmatory_mr/standardized_outcomes"
AUDIT=ROOT/"results/metabolic_resilience/stage3_full_pgwas/audit"
OUT.mkdir(parents=True,exist_ok=True)

def split_line(s,tab):
    return s.rstrip("\n").split("\t") if tab else s.split()

def open_text(path,member=""):
    path=Path(path)
    if path.suffix.lower()==".zip":
        z=zipfile.ZipFile(path)
        return io.TextIOWrapper(z.open(member),encoding="utf-8",errors="replace"),z
    if str(path).endswith(".gz"):
        return gzip.open(path,"rt",encoding="utf-8",errors="replace"),None
    return open(path,"r",encoding="utf-8",errors="replace"),None

if not INST.exists() or not REG.exists():
    print("[FAIL] annotated instrument or outcome registry missing")
    sys.exit(2)

with INST.open("r",encoding="utf-8",newline="") as f:
    inst_rows=list(csv.DictReader(f,delimiter="\t"))

# Unique canonical IVs; primary/stringent membership is restored in D2 from INST.
uniq_inst={}
for r in inst_rows:
    uniq_inst.setdefault(r["variant_id"],r)

by_rsid={}
by_pos={}
for vid,r in uniq_inst.items():
    rsid=r.get("rsid","")
    if rsid not in ("",".","NA"):
        by_rsid.setdefault(rsid,[]).append(vid)
    key=(str(r["chrom_hg19"]).replace("chr",""),int(r["pos_hg19"]))
    by_pos.setdefault(key,[]).append(vid)

with REG.open("r",encoding="utf-8",newline="") as f:
    registry=list(csv.DictReader(f,delimiter="\t"))

summary=[]
fail=0

for spec in registry:
    trait=spec["trait"]
    if spec["status"]!="PASS":
        print(f"[FAIL] {trait}: registry status={spec['status']}")
        summary.append({"trait":trait,"matched_unique_variants":0,"status":"FAIL_REGISTRY_NOT_PASS","file":""})
        fail+=1
        continue

    path=spec["path"]; member=spec.get("member","")
    fh,z=open_text(path,member)
    hits={}
    duplicate_conflicts=0
    header=None; H=None; tab=True

    try:
        for line in fh:
            s=line.strip()
            if not s or s.startswith("##"):
                continue
            if header is None:
                tab="\t" in line
                header=split_line(line.lstrip("#"),tab)
                H={x:i for i,x in enumerate(header)}
                continue

            p=split_line(line,tab)
            if len(p)<len(header):
                continue

            def val(col):
                return p[H[col]] if col and col in H and H[col]<len(p) else ""

            rsid=val(spec["col_rsid"]).strip()
            chrom=val(spec["col_chr"]).replace("chr","") if spec["col_chr"] else ""
            pos=None
            if spec["col_pos"]:
                try: pos=int(float(val(spec["col_pos"])))
                except Exception: pos=None

            matched_vids=[]
            match_method=""

            if rsid and rsid in by_rsid:
                matched_vids=by_rsid[rsid]
                match_method="RSID"
            elif spec["build"]!="RSID_PREFERRED" and chrom and pos is not None:
                matched_vids=by_pos.get((chrom,pos),[])
                if matched_vids:
                    match_method="GRCH37_POSITION"

            if not matched_vids:
                continue

            ea=val(spec["col_ea"]).upper()
            oa=val(spec["col_oa"]).upper()
            if not ea or not oa:
                continue

            beta=None
            if spec["col_beta"]:
                try: beta=float(val(spec["col_beta"]))
                except Exception: beta=None
            elif spec["col_or"]:
                try:
                    orv=float(val(spec["col_or"]))
                    beta=math.log(orv) if orv>0 else None
                except Exception: beta=None

            try: se=float(val(spec["col_se"]))
            except Exception: se=None
            if beta is None or se is None or se<=0:
                continue

            def fnum(col):
                if not col: return None
                try: return float(val(col))
                except Exception: return None

            for vid in matched_vids:
                row={
                    "trait":trait,"domain":spec["domain"],"outcome_type":spec["type"],
                    "source_build":spec["build"],"matched_variant_id":vid,
                    "matched_instrument_rsid":uniq_inst[vid].get("rsid",""),
                    "match_method":match_method,"rsid":rsid,
                    "chrom_source":chrom,"pos_source":pos if pos is not None else "",
                    "effect_allele":ea,"other_allele":oa,"beta":beta,"se":se,
                    "p":fnum(spec["col_p"]),"eaf":fnum(spec["col_eaf"]),
                    "n":fnum(spec["col_n"]),"source_path":path,
                }
                if vid in hits:
                    old=hits[vid]
                    same=(old["effect_allele"],old["other_allele"],old["beta"],old["se"])==(ea,oa,beta,se)
                    if not same:
                        duplicate_conflicts+=1
                        continue
                hits[vid]=row
    finally:
        fh.close()
        if z is not None: z.close()

    matched=list(hits.values())
    dst=OUT/f"{trait}.instrument_outcomes.tsv.gz"
    fields=[
        "trait","domain","outcome_type","source_build","matched_variant_id",
        "matched_instrument_rsid","match_method","rsid","chrom_source","pos_source",
        "effect_allele","other_allele","beta","se","p","eaf","n","source_path"
    ]
    with gzip.open(dst,"wt",encoding="utf-8",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
        wr.writeheader(); wr.writerows(matched)

    status="PASS" if matched and duplicate_conflicts==0 else ("FAIL_DUPLICATE_CONFLICT" if duplicate_conflicts else "FAIL_NO_MATCH")
    if status!="PASS": fail+=1
    summary.append({
        "trait":trait,"matched_unique_variants":len(matched),
        "instrument_universe_unique":len(uniq_inst),
        "coverage_pct":round(100*len(matched)/len(uniq_inst),3) if uniq_inst else 0,
        "duplicate_conflicts":duplicate_conflicts,"status":status,"file":str(dst)
    })
    print(f"[{status}] {trait}: matched={len(matched)}/{len(uniq_inst)} duplicate_conflicts={duplicate_conflicts}")

sumfile=AUDIT/"STAGE3D1_OUTCOME_INSTRUMENT_EXTRACTION.tsv"
with sumfile.open("w",encoding="utf-8",newline="") as f:
    fields=list(summary[0].keys()) if summary else ["status"]
    wr=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
    wr.writeheader(); wr.writerows(summary)

(AUDIT/"STAGE3D1_OUTCOME_INSTRUMENT_EXTRACTION.json").write_text(
    json.dumps({"traits":len(summary),"pass":sum(r["status"]=="PASS" for r in summary),
                "fail":sum(r["status"]!="PASS" for r in summary)},indent=2)+"\n",encoding="utf-8"
)
sys.exit(0 if fail==0 else 2)
