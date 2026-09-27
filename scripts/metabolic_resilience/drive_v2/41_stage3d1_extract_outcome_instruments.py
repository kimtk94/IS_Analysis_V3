#!/usr/bin/env python3
from pathlib import Path
import csv
import gzip
import zipfile
import math
import json
import os
import sys
import io

ROOT=Path(os.environ.get("IS_ANALYSIS_ROOT","/srv/is-analysis"))
INST=ROOT/"results/metabolic_resilience/stage3_full_pgwas/instruments/STAGE3C4_FROZEN_INSTRUMENTS_ANNOTATED.tsv"
REG=ROOT/"results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D0_OUTCOME_REGISTRY.tsv"
OUT=ROOT/"results/metabolic_resilience/stage3_confirmatory_mr/standardized_outcomes"
AUDIT=ROOT/"results/metabolic_resilience/stage3_full_pgwas/audit"
DUPAUDIT=AUDIT/"STAGE3D1_DUPLICATE_RESOLUTION.tsv"
OUT.mkdir(parents=True,exist_ok=True)

COMP=str.maketrans("ACGT","TGCA")

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

def snp_comp(a):
    return a.translate(COMP) if len(a)==1 and set(a)<=set("ACGT") else None

def allele_orientation(inst_ea,inst_oa,out_ea,out_oa):
    ie=inst_ea.upper(); ioa=inst_oa.upper()
    oe=out_ea.upper(); oo=out_oa.upper()

    if oe==ie and oo==ioa:
        return "DIRECT"
    if oe==ioa and oo==ie:
        return "SWAPPED"

    # Strand flips are valid only for SNPs. Never complement indels.
    if all(len(x)==1 for x in (ie,ioa,oe,oo)):
        ce=snp_comp(oe); co=snp_comp(oo)
        if ce==ie and co==ioa:
            return "STRAND_DIRECT"
        if ce==ioa and co==ie:
            return "STRAND_SWAPPED"

    return ""

def normalized_beta(row,orientation):
    b=row["beta"]
    if orientation in ("DIRECT","STRAND_DIRECT"):
        return b
    if orientation in ("SWAPPED","STRAND_SWAPPED"):
        return -b
    return None

def almost_equal(a,b,rtol=1e-10,atol=1e-12):
    return abs(a-b) <= max(atol,rtol*max(abs(a),abs(b),1.0))

if not INST.exists() or not REG.exists():
    print("[FAIL] annotated instrument or outcome registry missing")
    raise SystemExit(2)

with INST.open("r",encoding="utf-8",newline="") as f:
    inst_rows=list(csv.DictReader(f,delimiter="\t"))

# Unique canonical IVs; primary/stringent membership is restored in D2.
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
duplicate_audit=[]
fail=0

for spec in registry:
    trait=spec["trait"]

    if spec["status"]!="PASS":
        print(f"[FAIL] {trait}: registry status={spec['status']}")
        summary.append({
            "trait":trait,"matched_unique_variants":0,
            "instrument_universe_unique":len(uniq_inst),"coverage_pct":0,
            "raw_candidate_rows":0,"allele_incompatible_rows":0,
            "duplicate_variant_keys":0,"duplicate_rows_resolved":0,
            "unresolved_duplicate_conflicts":0,
            "status":"FAIL_REGISTRY_NOT_PASS","file":""
        })
        fail+=1
        continue

    path=spec["path"]
    member=spec.get("member","")
    fh,z=open_text(path,member)

    # Collect all allele-compatible source rows by frozen variant ID first,
    # then resolve duplicates deterministically.
    candidates_by_vid={}
    raw_candidate_rows=0
    allele_incompatible_rows=0
    header=None
    H=None
    tab=True

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
                try:
                    pos=int(float(val(spec["col_pos"])))
                except Exception:
                    pos=None

            # Build-aware identity rule:
            # - GRCh37 sources: exact GRCh37 position first; rsID only fallback.
            # - RSID_PREFERRED sources: rsID only, never cross-build position fallback.
            matched_vids=[]
            match_method=""

            if spec["build"]=="RSID_PREFERRED":
                if rsid and rsid in by_rsid:
                    matched_vids=by_rsid[rsid]
                    match_method="RSID"
            else:
                if chrom and pos is not None:
                    matched_vids=by_pos.get((chrom,pos),[])
                    if matched_vids:
                        match_method="GRCH37_POSITION"
                if not matched_vids and rsid and rsid in by_rsid:
                    matched_vids=by_rsid[rsid]
                    match_method="RSID_FALLBACK"

            if not matched_vids:
                continue

            ea=val(spec["col_ea"]).upper()
            oa=val(spec["col_oa"]).upper()
            if not ea or not oa:
                continue

            beta=None
            if spec["col_beta"]:
                try:
                    beta=float(val(spec["col_beta"]))
                except Exception:
                    beta=None
            elif spec["col_or"]:
                try:
                    orv=float(val(spec["col_or"]))
                    beta=math.log(orv) if orv>0 else None
                except Exception:
                    beta=None

            try:
                se=float(val(spec["col_se"]))
            except Exception:
                se=None

            if beta is None or not math.isfinite(beta) or se is None or not math.isfinite(se) or se<=0:
                continue

            def fnum(col):
                if not col:
                    return None
                try:
                    x=float(val(col))
                    return x if math.isfinite(x) else None
                except Exception:
                    return None

            for vid in matched_vids:
                raw_candidate_rows+=1
                inst=uniq_inst[vid]
                inst_ea=inst["exposure_effect_allele"].upper()
                inst_oa=inst["exposure_other_allele"].upper()

                orientation=allele_orientation(inst_ea,inst_oa,ea,oa)
                if not orientation:
                    allele_incompatible_rows+=1
                    continue

                row={
                    "trait":trait,
                    "domain":spec["domain"],
                    "outcome_type":spec["type"],
                    "source_build":spec["build"],
                    "matched_variant_id":vid,
                    "matched_instrument_rsid":inst.get("rsid",""),
                    "match_method":match_method,
                    "rsid":rsid,
                    "chrom_source":chrom,
                    "pos_source":pos if pos is not None else "",
                    "effect_allele":ea,
                    "other_allele":oa,
                    "beta":beta,
                    "se":se,
                    "p":fnum(spec["col_p"]),
                    "eaf":fnum(spec["col_eaf"]),
                    "n":fnum(spec["col_n"]),
                    "source_path":path,
                    "d1_allele_orientation":orientation,
                }

                candidates_by_vid.setdefault(vid,[]).append(row)

    finally:
        fh.close()
        if z is not None:
            z.close()

    hits={}
    duplicate_variant_keys=0
    duplicate_rows_resolved=0
    unresolved_duplicate_conflicts=0

    for vid,arr in candidates_by_vid.items():
        if len(arr)==1:
            hits[vid]=arr[0]
            continue

        duplicate_variant_keys+=1

        # Collapse byte/logically identical association rows.
        uniq={}
        for x in arr:
            key=(
                x["effect_allele"],x["other_allele"],
                x["beta"],x["se"],x["p"],x["eaf"],x["n"],
                x["chrom_source"],x["pos_source"],x["rsid"],
            )
            uniq.setdefault(key,x)

        dedup=list(uniq.values())

        if len(dedup)==1:
            hits[vid]=dedup[0]
            duplicate_rows_resolved+=len(arr)-1
            duplicate_audit.append({
                "trait":trait,"variant_id":vid,
                "raw_rows":len(arr),"compatible_unique_rows":1,
                "resolution":"IDENTICAL_DUPLICATES_COLLAPSED",
                "selected_beta":dedup[0]["beta"],
                "selected_se":dedup[0]["se"],
            })
            continue

        # Normalize effect direction to the frozen exposure EA and check whether
        # rows are equivalent representations of the same association.
        normalized=[]
        for x in dedup:
            nb=normalized_beta(x,x["d1_allele_orientation"])
            normalized.append((x,nb))

        ref_beta=normalized[0][1]
        ref_se=normalized[0][0]["se"]

        equivalent=(
            ref_beta is not None
            and all(nb is not None and almost_equal(nb,ref_beta) for _,nb in normalized)
            and all(almost_equal(x["se"],ref_se) for x,_ in normalized)
        )

        if equivalent:
            # Prefer exact GRCh37-position match, then direct orientation,
            # then deterministic lexical ordering.
            inst=uniq_inst[vid]
            ipos=int(inst["pos_hg19"])

            def rank(x):
                exact_pos=(
                    str(x["pos_source"]) not in ("","None")
                    and int(float(x["pos_source"]))==ipos
                ) if str(x["pos_source"]) not in ("","None") else False
                direct=x["d1_allele_orientation"] in ("DIRECT","STRAND_DIRECT")
                return (
                    0 if exact_pos else 1,
                    0 if direct else 1,
                    x["effect_allele"],
                    x["other_allele"],
                )

            chosen=sorted(dedup,key=rank)[0]
            hits[vid]=chosen
            duplicate_rows_resolved+=len(arr)-1
            duplicate_audit.append({
                "trait":trait,"variant_id":vid,
                "raw_rows":len(arr),"compatible_unique_rows":len(dedup),
                "resolution":"EQUIVALENT_AFTER_ALLELE_NORMALIZATION",
                "selected_beta":chosen["beta"],
                "selected_se":chosen["se"],
            })
            continue

        unresolved_duplicate_conflicts+=1
        duplicate_audit.append({
            "trait":trait,"variant_id":vid,
            "raw_rows":len(arr),"compatible_unique_rows":len(dedup),
            "resolution":"UNRESOLVED_CONFLICT",
            "selected_beta":"",
            "selected_se":"",
        })

    matched=list(hits.values())

    dst=OUT/f"{trait}.instrument_outcomes.tsv.gz"
    fields=[
        "trait","domain","outcome_type","source_build","matched_variant_id",
        "matched_instrument_rsid","match_method","rsid","chrom_source","pos_source",
        "effect_allele","other_allele","beta","se","p","eaf","n","source_path",
        "d1_allele_orientation"
    ]

    with gzip.open(dst,"wt",encoding="utf-8",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
        wr.writeheader()
        wr.writerows(matched)

    if not matched:
        status="FAIL_NO_MATCH"
    elif unresolved_duplicate_conflicts>0:
        status="FAIL_DUPLICATE_CONFLICT"
    else:
        status="PASS"

    if status!="PASS":
        fail+=1

    summary.append({
        "trait":trait,
        "matched_unique_variants":len(matched),
        "instrument_universe_unique":len(uniq_inst),
        "coverage_pct":round(100*len(matched)/len(uniq_inst),3) if uniq_inst else 0,
        "raw_candidate_rows":raw_candidate_rows,
        "allele_incompatible_rows":allele_incompatible_rows,
        "duplicate_variant_keys":duplicate_variant_keys,
        "duplicate_rows_resolved":duplicate_rows_resolved,
        "unresolved_duplicate_conflicts":unresolved_duplicate_conflicts,
        "status":status,
        "file":str(dst),
    })

    print(
        f"[{status}] {trait}: matched={len(matched)}/{len(uniq_inst)} "
        f"raw={raw_candidate_rows} incompatible={allele_incompatible_rows} "
        f"duplicate_keys={duplicate_variant_keys} resolved={duplicate_rows_resolved} "
        f"unresolved={unresolved_duplicate_conflicts}"
    )

sumfile=AUDIT/"STAGE3D1_OUTCOME_INSTRUMENT_EXTRACTION.tsv"
with sumfile.open("w",encoding="utf-8",newline="") as f:
    fields=list(summary[0].keys()) if summary else ["status"]
    wr=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
    wr.writeheader()
    wr.writerows(summary)

dup_fields=[
    "trait","variant_id","raw_rows","compatible_unique_rows",
    "resolution","selected_beta","selected_se"
]
with DUPAUDIT.open("w",encoding="utf-8",newline="") as f:
    wr=csv.DictWriter(f,fieldnames=dup_fields,delimiter="\t",lineterminator="\n")
    wr.writeheader()
    wr.writerows(duplicate_audit)

payload={
    "traits":len(summary),
    "pass":sum(r["status"]=="PASS" for r in summary),
    "fail":sum(r["status"]!="PASS" for r in summary),
    "duplicate_variant_keys":sum(r["duplicate_variant_keys"] for r in summary),
    "duplicate_rows_resolved":sum(r["duplicate_rows_resolved"] for r in summary),
    "unresolved_duplicate_conflicts":sum(r["unresolved_duplicate_conflicts"] for r in summary),
    "rule":"Duplicates are resolved only when frozen-IV allele compatibility identifies a single/equivalent association representation. Genuine discordant associations remain FAIL.",
}
(AUDIT/"STAGE3D1_OUTCOME_INSTRUMENT_EXTRACTION.json").write_text(
    json.dumps(payload,indent=2)+"\n",encoding="utf-8"
)
print(json.dumps(payload,indent=2))

raise SystemExit(0 if fail==0 else 2)
