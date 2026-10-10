#!/usr/bin/env python3
"""Direct variant-level AIS GWAS lookup of exploratory alcohol SuSiE CS SNPs.

Works on *canonical GRCh37* BBJ-IS Japanese and GIGASTROKE EAS AIS files.
Allele orientation is aligned to the INPUT ALT (not alcohol exposure effect).
The GWAS input source build is NOT considered resolved by this check.
Missing SNPs are MISSING_FROM_DATASET, not negative associations.
This is diagnostic association overlap, never evidence of causal mediation.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

DEFAULT_BBJ="/srv/is-analysis/data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz"
DEFAULT_GIGA="/srv/is-analysis/data/is/processed/gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz"
SOURCES={"BBJ_JAPAN_IS":DEFAULT_BBJ, "GIGASTROKE_EAS_AIS":DEFAULT_GIGA}

def number(s):
    try:
        value=float(s)
        return value if math.isfinite(value) else None
    except (TypeError,ValueError,OverflowError):
        return None

def harmonize_alt(row, alt, ref):
    ea=str(row.get("effect_allele","")).upper()
    oa=str(row.get("other_allele","")).upper()
    b=number(row.get("beta"))
    se=number(row.get("se"))
    af=number(row.get("eaf"))
    if ea==alt and oa==ref: sign=1;orientation="EFFECT_ALT"
    elif ea==ref and oa==alt: sign=-1;orientation="EFFECT_REF_FLIPPED"
    else: return {"orientation":"INVALID_EFFECT_ALLELES","beta_alt":"","se":"","eaf_alt":"","or_alt":""}
    beta=b*sign if b is not None else None
    eaf=af if sign==1 else (1-af if af is not None else None)
    return {"orientation":orientation,"beta_alt":beta if beta is not None else "",
      "se":se if se is not None else "","eaf_alt":eaf if eaf is not None else "",
      "or_alt":math.exp(beta) if beta is not None and -700<beta<700 else ""}

def scan(path, targets):
    """Fast gzip+grep chromosome-position prefilter, then allele-level validation.

    Returned count is number of target-position rows, not full GWAS row count.
    """
    import subprocess
    import tempfile

    position_map=defaultdict(list)
    for x in targets: position_map[(x["chrom"],int(x["pos"]))].append(x)
    matched=defaultdict(list)
    colocated=defaultdict(list)
    opener=gzip.open if str(path).endswith(".gz") else open
    with opener(path,"rt",newline="") as f:
        header=f.readline().rstrip("\n\r").split("\t")
    required={"chr","pos","effect_allele","other_allele","beta","se","p","build","variant_id"}
    missing=required-set(header)
    if missing: raise ValueError(f"{path}: missing columns {sorted(missing)}")

    with tempfile.NamedTemporaryFile("w",prefix="is_ais_target_positions_",delete=True) as pats:
        for chrom,pos in position_map:
            pats.write(f"\t{chrom}\t{pos}\t\n")
        pats.flush()
        source=(["gzip","-dc",str(path)] if str(path).endswith(".gz")
                else ["cat",str(path)])
        decomp=subprocess.Popen(source,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        filt=subprocess.Popen(["grep","-F","-f",pats.name],stdin=decomp.stdout,
                              stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        decomp.stdout.close()
        count=0
        try:
            for raw in filt.stdout:
                count+=1
                values=raw.rstrip("\n\r").split("\t")
                if len(values)!=len(header):
                    raise ValueError(f"GWAS row has {len(values)} fields, expected {len(header)}")
                line=dict(zip(header,values))
                chrom=str(line["chr"]).removeprefix("chr").upper()
                pos=int(line["pos"])
                opts=position_map.get((chrom,pos),[])
                if not opts: continue
                if line["build"]!="GRCh37":
                    raise ValueError(f"{path}: targeted SNP build is {line['build']}, not GRCh37")
                if not ("ref" in line and "alt" in line):
                    parts=line["variant_id"].split(":")
                    if len(parts)!=4 or parts[0].removeprefix("chr")!=chrom or int(parts[1])!=pos:
                        raise ValueError(f"{path}: cannot derive REF/ALT from variant_id")
                    line["ref"],line["alt"]=parts[2],parts[3]
                for target in opts:
                    vid=target["variant_id"]
                    if line["ref"].upper()==target["ref"].upper() and line["alt"].upper()==target["alt"].upper():
                        matched[vid].append(line)
                    else:
                        colocated[vid].append(line)
        except Exception:
            filt.kill()
            decomp.kill()
            raise
        finally:
            filt.stdout.close()
            filter_stderr=filt.stderr.read()
            filt.wait()
            source_err=decomp.stderr.read().decode(errors="replace")
            decomp.wait()
            filt.stderr.close()
            decomp.stderr.close()
        if filt.returncode not in (0,1) or decomp.returncode!=0:
            raise RuntimeError(f"scan subprocess failed: grep={filt.returncode}, decompress={decomp.returncode}, stderr={filter_stderr[:200]} {source_err[:200]}")
    return matched,colocated,count

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--audit-tsv",type=Path,required=True)
    p.add_argument("--bbj",type=Path,default=Path(DEFAULT_BBJ))
    p.add_argument("--giga",type=Path,default=Path(DEFAULT_GIGA))
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    with a.audit_tsv.open(newline="") as f: targets=list(csv.DictReader(f,delimiter="\t"))
    if len(targets)!=11 or len({x["variant_id"] for x in targets})!=11:
        raise ValueError("Expected exactly 11 distinct curated SNPs; fail closed")
    for t in targets:
        if t["allele_mapping_status"]!="VERIFIED_ALLELE_SET":
            raise ValueError("Not all rsID-allele verifications passed")
        if t["ref"].upper() not in ("A","C","G","T") or t["alt"].upper() not in ("A","C","G","T"):
            raise ValueError("Non-SNP or ambiguous input allele")
    sources={"BBJ_JAPAN_IS":a.bbj,"GIGASTROKE_EAS_AIS":a.giga}
    for src in sources.values():
        if not src.is_file(): raise FileNotFoundError(src)
    if a.out_dir.exists() and any(a.out_dir.iterdir()):
        raise FileExistsError(f"Output dir not empty (do not overwrite): {a.out_dir}")
    a.out_dir.mkdir(parents=True,exist_ok=True)
    outrows=[]; meta={}
    for dataset,path in sources.items():
        found,other,count=scan(path,targets)
        meta[dataset]={"source":str(path),"target_position_rows":count}
        for t in targets:
            vid=t["variant_id"]
            hits=found.get(vid,[])
            alt_hits=other.get(vid,[])
            reversed_pairs=[x for x in alt_hits
                if x["ref"].upper()==t["alt"].upper() and x["alt"].upper()==t["ref"].upper()]
            status=("MATCH_EXACT" if len(hits)==1 else
                    "DUPLICATE_MATCH_REVIEW" if len(hits)>1 else
                    "REF_ALT_SWAP_REVIEW" if len(reversed_pairs)==1 else
                    "DUPLICATE_SWAP_REVIEW" if len(reversed_pairs)>1 else
                    "ALLELE_MISMATCH_SAME_POSITION" if alt_hits else
                    "MISSING_FROM_DATASET")
            # Swapped REF/ALT is NOT a canonical GRCh37 REF match; flag separately.
            hit=(hits[0] if len(hits)==1 else
                 reversed_pairs[0] if not hits and len(reversed_pairs)==1 else {})
            h=harmonize_alt(hit,t["alt"].upper(),t["ref"].upper()) if hit else {}
            if len(hits)==1 and h.get("orientation")=="INVALID_EFFECT_ALLELES":
                status="EFFECT_ALLELE_INVALID"
            row={"locus":t["locus"],"variant_id":vid,"rsid_verified":t["allele_verified_rsids"],
              "consequence":t["vep_consequence"],"dataset":dataset,
              "build":"GRCh37","status":status,"n_exact_rows":len(hits),
              "n_other_allele_rows":len(alt_hits),"n_ref_alt_swaps":len(reversed_pairs),
              "effect_allele":hit.get("effect_allele",""),
              "other_allele":hit.get("other_allele",""),
              "orientation":h.get("orientation",""),
              "beta_alt":h.get("beta_alt",""),"se":h.get("se",""),
              "p":hit.get("p",""),"eaf_alt":h.get("eaf_alt",""),
              "or_alt":h.get("or_alt",""),
              "sample_n":hit.get("n",""),
              "gwas_rsid":hit.get("rsid",""),
              "source_variant_id":hit.get("variant_id","")}
            # RSID field may be stale/ambiguous in GWAS; record but do not require exact equality
            if hit and hit.get("variant_id","")!=vid and status=="MATCH_EXACT":
                row["status"]="INPUT_VARIANT_ID_MISMATCH"
            outrows.append(row)
        print(f"SCANNED {dataset} target_position_rows={count} exact={sum(r['status']=='MATCH_EXACT' for r in outrows if r['dataset']==dataset)}",flush=True)
    dest=a.out_dir/"ALCOHOL_CS_AIS_GWAS_DIRECT_OVERLAP.tsv"
    with dest.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(outrows[0]),delimiter="\t");w.writeheader();w.writerows(outrows)
    manifest={"analysis":"ALCOHOL_CS_AIS_DIRECT_OVERLAP","timestamp_utc":datetime.now(timezone.utc).isoformat(),
      "interpretation":"OBSERVATIONAL_ASSOCIATION_LOOKUP_ONLY_NOT_MR_NOT_COLOCALIZATION",
      "source_alcohol_gwas_build":"UNRESOLVED_OFFICIAL_SOURCE_BUILD",
      "coordinates":"GRCh37_FASTA_11_OF_11_REF_MATCH_USER_REPORTED",
      "rsid":"VEP_GRCh37_ALLELE_SET_11_OF_11_MATCH",
      "effect":"AIS GWAS effect harmonized to FASTA-verified input ALT; not aligned to alcohol GWAS effect",
      "ref_alt_swaps":"REF_ALT_SWAP_REVIEW rows do not have canonical GRCh37 REF consistency; beta and p shown provisionally",
      "sources":meta,"total_pairs":len(outrows),
      "by_dataset":{ds:{s:sum(r["dataset"]==ds and r["status"]==s for r in outrows)
                          for s in sorted({r["status"] for r in outrows if r["dataset"]==ds})}
                    for ds in sources}}
    (a.out_dir/"ALCOHOL_CS_AIS_GWAS_DIRECT_OVERLAP_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print("AIS_DIRECT_OVERLAP_COMPLETE",json.dumps(manifest["by_dataset"]),flush=True)
    for r in outrows:
        print("RESULT",r["dataset"],r["locus"],r["rsid_verified"],r["status"],
              f"p={r['p']}",f"beta_alt={r['beta_alt']}",flush=True)

if __name__=="__main__": main()
