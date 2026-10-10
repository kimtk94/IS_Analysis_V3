#!/usr/bin/env python3
"""Prepare a fail-closed, build-unresolved variant manifest from SuSiE CS TSV."""
import argparse,csv,json,re
from collections import defaultdict
from pathlib import Path

PAT=re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|MT):([1-9][0-9]*):([ACGT]+):([ACGT]+)$",re.I)
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("input_tsv",type=Path)
    ap.add_argument("--out-dir",type=Path)
    args=ap.parse_args()
    out=args.out_dir or args.input_tsv.parent
    out.mkdir(parents=True,exist_ok=True)
    seen={}
    groups=defaultdict(set)
    with args.input_tsv.open(newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        needed={"locus","model","cs","variant_id","pip"}
        if not needed.issubset(reader.fieldnames or []): raise ValueError("Missing required columns")
        for row in reader:
            raw=row["variant_id"].strip()
            m=PAT.fullmatch(raw)
            if not m: raise ValueError("Unsupported SNP ID/alleles: "+raw)
            chrom,pos,ref,alt=m.groups()
            chrom=chrom.upper();ref=ref.upper();alt=alt.upper()
            if ref==alt: raise ValueError("Identical REF/ALT: "+raw)
            canonical=f"{chrom}:{int(pos)}:{ref}:{alt}"
            pip=float(row["pip"])
            if not 0<=pip<=1: raise ValueError("Invalid PIP: "+raw)
            locus=row["locus"]
            if locus not in {"ADH1B","ALDH2"}: raise ValueError("Unexpected locus")
            if canonical in seen and seen[canonical]["locus"]!=locus:
                raise ValueError("Cross-locus ID collision")
            seen.setdefault(canonical,{"locus":locus,"chrom":chrom,"pos":int(pos),
              "ref":ref,"alt":alt,"variant_id":canonical,"genome_build":"UNRESOLVED",
              "rsid":"UNRESOLVED","vep_consequence":"NOT_RUN",
              "annotation_status":"BLOCKED_PENDING_BUILD"})
            groups[canonical].add((row["model"],row["cs"]))
    dest=out/"ALCOHOL_CS_VARIANTS_BUILD_UNRESOLVED.tsv"
    columns=["locus","variant_id","chrom","pos","ref","alt","genome_build",
             "rsid","vep_consequence","annotation_status","cs_memberships"]
    with dest.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter="\t")
        w.writeheader()
        for key,item in sorted(seen.items(),key=lambda x:(x[1]["chrom"],x[1]["pos"])):
            w.writerow({**item,"cs_memberships":";".join(sorted(a+":"+b for a,b in groups[key]))})
    manifest={"status":"BUILD_UNRESOLVED_NO_VEP_REQUEST","unique_variants":len(seen),
      "source":str(args.input_tsv),"output":str(dest),
      "prohibited_inferences":["rsid","genome_build","consequence","AIS_causal_variant"]}
    (out/"ALCOHOL_CS_ANNOTATION_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print("ANNOTATION_MANIFEST_OK",len(seen),dest)
if __name__=="__main__":main()
