#!/usr/bin/env python3
"""Prepare locus-matched UKB-PPP pQTL and CKDGen eGFR data for colocalization.

Key design points:
- UKB-PPP full per-protein summaries are GRCh38.
- CKDGen Stanzick 2021 eGFR is GRCh37/hg19.
- Extract each protein cis region in GRCh38, lift SNP coordinates to GRCh37,
  then match the outcome by GRCh37 chromosome/position + allele pair.
- Only biallelic SNVs are used in this initial coloc stage. This avoids
  ambiguous indel normalization across builds.
- UKB-PPP ALLELE1 is treated as the effect allele for BETA; A1FREQ is its
  frequency. Outcome beta is aligned to the pQTL ALLELE1 where possible.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import re
import tarfile
from collections import defaultdict
from pathlib import Path

BASES = set("ACGT")
COMP = {"A":"T","T":"A","C":"G","G":"C"}


def fnum(x):
    if x is None:
        return None
    s=str(x).strip()
    if not s or s.lower() in {"na","nan",".","none"}:
        return None
    try:
        v=float(s)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def norm_chr(x):
    s=str(x).strip()
    if s.lower().startswith("chr"):
        s=s[3:]
    if s=="23":
        return "X"
    return s.upper()


def norm_a(x):
    return str(x).strip().upper()


def complement(a):
    return "".join(COMP[b] for b in a) if a and all(b in COMP for b in a) else None


def is_snv(a0,a1):
    return len(a0)==1 and len(a1)==1 and a0 in BASES and a1 in BASES and a0!=a1


def p_from_log10(logp):
    if logp is None:
        return None
    if logp >= 323:
        return 0.0
    return 10.0 ** (-logp)


def read_tsv(path: Path):
    opener=gzip.open if path.suffix==".gz" else open
    with opener(path,"rt",encoding="utf-8",newline="") as fh:
        yield from csv.DictReader(fh,delimiter="\t")


def write_tsv_gz(path: Path, rows, fields):
    path.parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(path,"wt",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def load_candidates(path):
    rows=list(read_tsv(path))
    if not rows:
        raise SystemExit(f"no candidates: {path}")
    return rows


def load_gene_coords(path):
    out={}
    for r in read_tsv(path):
        gene=(r.get("gene_symbol") or r.get("gene") or "").upper()
        if not gene:
            continue
        build=(r.get("genome_build") or "").upper()
        if build and build not in {"GRCH38","HG38"}:
            continue
        out[gene]={
            "chr":norm_chr(r.get("chr","")),
            "start":int(float(r["start"])),
            "end":int(float(r["end"])),
        }
    return out


def archive_for_gene(root: Path, gene: str):
    hits=sorted((root/"EUR"/gene).glob("*.tar"))
    if len(hits)!=1:
        raise SystemExit(f"{gene}: expected exactly one EUR .tar under {root/'EUR'/gene}, got {len(hits)}")
    return hits[0]


def chromosome_member(tar, chrom):
    pat=re.compile(rf"(^|/)discovery_chr{re.escape(chrom)}_",re.I)
    hits=[m for m in tar.getmembers() if m.isfile() and pat.search(m.name)]
    if len(hits)!=1:
        raise SystemExit(f"expected one chr{chrom} member, got {len(hits)}: {[m.name for m in hits[:10]]}")
    return hits[0]


def text_member(tar, member):
    raw=tar.extractfile(member)
    if raw is None:
        raise RuntimeError(f"cannot extract {member.name}")
    binary=raw
    if member.name.lower().endswith(".gz"):
        binary=gzip.GzipFile(fileobj=raw,mode="rb")
    return io.TextIOWrapper(binary,encoding="utf-8",errors="replace",newline="")


def split_ws(line, n=14):
    vals=line.strip().split()
    if len(vals)<n:
        return None
    if len(vals)>n:
        vals=vals[:n-1]+[" ".join(vals[n-1:])]
    return vals


def extract_locus(candidate, coord, archive, lo, window_bp, min_info):
    chrom=coord["chr"]
    lo_start=max(1,coord["start"]-window_bp)
    lo_end=coord["end"]+window_bp
    rows=[]
    counts=defaultdict(int)

    with tarfile.open(archive,"r:*") as tar:
        member=chromosome_member(tar,chrom)
        with text_member(tar,member) as fh:
            header=fh.readline().strip().split()
            required=["CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA"]
            missing=[x for x in required if x not in header]
            if missing:
                raise SystemExit(f"{candidate['gene_symbol']}: missing pQTL columns {missing}; header={header}")
            idx={x:i for i,x in enumerate(header)}
            for line in fh:
                if not line.strip():
                    continue
                v=split_ws(line,len(header))
                if v is None:
                    counts["malformed"]+=1
                    continue
                if norm_chr(v[idx["CHROM"]]) != chrom:
                    continue
                pos38=fnum(v[idx["GENPOS"]])
                if pos38 is None:
                    counts["bad_pos"]+=1
                    continue
                pos38=int(pos38)
                if pos38 < lo_start or pos38 > lo_end:
                    continue
                counts["in_region"]+=1
                if v[idx["TEST"]].upper()!="ADD":
                    counts["non_add"]+=1
                    continue
                a0,a1=norm_a(v[idx["ALLELE0"]]),norm_a(v[idx["ALLELE1"]])
                if not is_snv(a0,a1):
                    counts["non_snv"]+=1
                    continue
                af=fnum(v[idx["A1FREQ"]]); info=fnum(v[idx["INFO"]])
                beta=fnum(v[idx["BETA"]]); se=fnum(v[idx["SE"]]); n=fnum(v[idx["N"]]); logp=fnum(v[idx["LOG10P"]])
                if af is None or beta is None or se is None or se<=0 or n is None or not (0<af<1):
                    counts["bad_stats"]+=1
                    continue
                if info is not None and info < min_info:
                    counts["low_info"]+=1
                    continue

                lifts=lo.convert_coordinate("chr"+chrom,pos38-1)
                lifts=[x for x in lifts if norm_chr(x[0])==chrom]
                # Require a unique same-chromosome mapping.
                uniq={(norm_chr(x[0]),int(round(x[1]))+1,x[2]) for x in lifts}
                if len(uniq)!=1:
                    counts["liftover_nonunique_or_missing"]+=1
                    continue
                chr37,pos37,strand=next(iter(uniq))
                if strand=="-":
                    a0=complement(a0); a1=complement(a1)
                    if a0 is None or a1 is None:
                        counts["bad_reverse_allele"]+=1
                        continue

                rows.append({
                    "protein_id":candidate["protein_id"],
                    "gene_symbol":candidate["gene_symbol"],
                    "variant_id_grch38":v[idx["ID"]],
                    "chr38":chrom,
                    "pos38":pos38,
                    "chr37":chr37,
                    "pos37":pos37,
                    "allele0_pqtl_hg37":a0,
                    "allele1_pqtl_hg37":a1,
                    "a1freq_pqtl":af,
                    "maf_pqtl":min(af,1-af),
                    "info_pqtl":"" if info is None else info,
                    "n_pqtl":int(n),
                    "beta_pqtl":beta,
                    "se_pqtl":se,
                    "p_pqtl":p_from_log10(logp),
                    "liftover_strand":strand,
                    "source_member":member.name,
                })
                counts["kept"]+=1
    return rows,dict(counts),lo_start,lo_end


def outcome_header(path):
    with gzip.open(path,"rt",encoding="utf-8",errors="replace",newline="") as fh:
        reader=csv.DictReader(fh,delimiter="\t")
        return reader.fieldnames or []


def load_outcome_matches(path, wanted_positions):
    out=defaultdict(list)
    parsed=0
    with gzip.open(path,"rt",encoding="utf-8",errors="replace",newline="") as fh:
        reader=csv.DictReader(fh,delimiter="\t")
        required={"RSID","Allele1","Allele2","Effect","StdErr.GC","P.value.GC","Freq1","n","chr","pos"}
        missing=required-set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Stanzick eGFR missing columns: {sorted(missing)}; header={reader.fieldnames}")
        for r in reader:
            parsed+=1
            chrom=norm_chr(r.get("chr",""))
            pos=fnum(r.get("pos"))
            if pos is None:
                continue
            key=(chrom,int(pos))
            if key not in wanted_positions:
                continue
            ea,oa=norm_a(r.get("Allele1","")),norm_a(r.get("Allele2",""))
            if not is_snv(ea,oa):
                continue
            beta=fnum(r.get("Effect")); se=fnum(r.get("StdErr.GC")); eaf=fnum(r.get("Freq1")); n=fnum(r.get("n"))
            if beta is None or se is None or se<=0 or eaf is None or not (0<eaf<1) or n is None:
                continue
            out[key].append({
                "rsid":r.get("RSID","").strip(),
                "effect_allele_outcome":ea,
                "other_allele_outcome":oa,
                "beta_outcome":beta,
                "se_outcome":se,
                "p_outcome":fnum(r.get("P.value.GC")),
                "eaf_outcome":eaf,
                "maf_outcome":min(eaf,1-eaf),
                "n_outcome":int(n),
            })
    return out,parsed


def match_outcome(p, choices):
    a0,a1=p["allele0_pqtl_hg37"],p["allele1_pqtl_hg37"]
    same=[]
    for o in choices:
        ea,oa=o["effect_allele_outcome"],o["other_allele_outcome"]
        if {ea,oa}!={a0,a1}:
            continue
        orient="direct" if (ea,oa)==(a1,a0) else "swapped" if (ea,oa)==(a0,a1) else "pair_match"
        same.append((0 if orient=="direct" else 1,o,orient))
    if not same:
        return None
    same.sort(key=lambda x:(x[0], fnum(x[1].get("p_outcome")) if fnum(x[1].get("p_outcome")) is not None else 1.0))
    _,o,orient=same[0]
    sign=1 if orient=="direct" else -1
    return o,orient,sign


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--pqtl-root",type=Path,required=True)
    ap.add_argument("--gene-coordinates",type=Path,required=True)
    ap.add_argument("--chain",type=Path,required=True)
    ap.add_argument("--outcome-egfr",type=Path,required=True)
    ap.add_argument("--output-root",type=Path,required=True)
    ap.add_argument("--window-bp",type=int,default=1_000_000)
    ap.add_argument("--min-info",type=float,default=0.3)
    ap.add_argument("--min-shared-snps",type=int,default=50)
    args=ap.parse_args()

    try:
        from pyliftover import LiftOver
    except ImportError as exc:
        raise SystemExit("pyliftover required in CKD venv") from exc

    candidates=load_candidates(args.candidates)
    coords=load_gene_coords(args.gene_coordinates)
    lo=LiftOver(str(args.chain))
    extracted={}
    extraction_meta={}
    wanted=set()

    locus_fields=[
        "protein_id","gene_symbol","variant_id_grch38","chr38","pos38","chr37","pos37",
        "allele0_pqtl_hg37","allele1_pqtl_hg37","a1freq_pqtl","maf_pqtl","info_pqtl",
        "n_pqtl","beta_pqtl","se_pqtl","p_pqtl","liftover_strand","source_member",
    ]
    for c in candidates:
        gene=c["gene_symbol"].upper()
        if gene not in coords:
            raise SystemExit(f"{gene}: missing GRCh38 coordinate")
        coord=coords[gene]
        if norm_chr(c["anchor_chr_hg19"]) != coord["chr"]:
            raise SystemExit(f"{gene}: chromosome mismatch candidate={c['anchor_chr_hg19']} coord={coord['chr']}")
        archive=archive_for_gene(args.pqtl_root,gene)
        rows,counts,start38,end38=extract_locus(c,coord,archive,lo,args.window_bp,args.min_info)
        if not rows:
            raise SystemExit(f"{gene}: no pQTL SNVs survived locus extraction")
        extracted[gene]=rows
        wanted.update((r["chr37"],r["pos37"]) for r in rows)
        write_tsv_gz(args.output_root/"pqtl_locus"/f"{gene}.tsv.gz",rows,locus_fields)
        extraction_meta[gene]={
            "archive":str(archive),"chr38":coord["chr"],"gene_start38":coord["start"],"gene_end38":coord["end"],
            "region_start38":start38,"region_end38":end38,"counts":counts,
        }

    outcome,source_rows=load_outcome_matches(args.outcome_egfr,wanted)

    coloc_fields=[
        "protein_id","gene_symbol","snp","chr37","pos37",
        "allele0_pqtl","allele1_pqtl","a1freq_pqtl","maf_pqtl","n_pqtl","beta_pqtl","se_pqtl","p_pqtl",
        "effect_allele_outcome_original","other_allele_outcome_original","eaf_outcome","maf_outcome",
        "n_outcome","beta_outcome_original","beta_outcome_aligned","se_outcome","p_outcome","harmonization",
        "variant_id_grch38","chr38","pos38",
    ]
    qc=[]
    for c in candidates:
        gene=c["gene_symbol"].upper()
        matched=[]
        reasons=defaultdict(int)
        for p in extracted[gene]:
            m=match_outcome(p,outcome.get((p["chr37"],p["pos37"]),[]))
            if m is None:
                reasons["no_position_allele_match"]+=1
                continue
            o,orient,sign=m
            snp=o["rsid"] if o["rsid"].lower().startswith("rs") else f"{p['chr37']}:{p['pos37']}:{p['allele0_pqtl_hg37']}:{p['allele1_pqtl_hg37']}"
            matched.append({
                "protein_id":p["protein_id"],"gene_symbol":gene,"snp":snp,
                "chr37":p["chr37"],"pos37":p["pos37"],
                "allele0_pqtl":p["allele0_pqtl_hg37"],"allele1_pqtl":p["allele1_pqtl_hg37"],
                "a1freq_pqtl":p["a1freq_pqtl"],"maf_pqtl":p["maf_pqtl"],"n_pqtl":p["n_pqtl"],
                "beta_pqtl":p["beta_pqtl"],"se_pqtl":p["se_pqtl"],"p_pqtl":p["p_pqtl"],
                "effect_allele_outcome_original":o["effect_allele_outcome"],
                "other_allele_outcome_original":o["other_allele_outcome"],
                "eaf_outcome":o["eaf_outcome"],"maf_outcome":o["maf_outcome"],"n_outcome":o["n_outcome"],
                "beta_outcome_original":o["beta_outcome"],"beta_outcome_aligned":o["beta_outcome"]*sign,
                "se_outcome":o["se_outcome"],"p_outcome":o["p_outcome"],"harmonization":orient,
                "variant_id_grch38":p["variant_id_grch38"],"chr38":p["chr38"],"pos38":p["pos38"],
            })
            reasons["matched"]+=1

        # coloc requires unique SNP identifiers; keep the row with lowest pQTL p when duplicates occur.
        best={}
        for r in matched:
            key=r["snp"]
            score=fnum(r["p_pqtl"])
            if key not in best or (score is not None and score < (fnum(best[key]["p_pqtl"]) or 1.0)):
                best[key]=r
        matched=list(best.values())
        matched.sort(key=lambda r:(int(r["pos37"]),r["snp"]))
        write_tsv_gz(args.output_root/"coloc_input"/f"{gene}.tsv.gz",matched,coloc_fields)

        qc.append({
            "gene_symbol":gene,
            "protein_id":c["protein_id"],
            "pqtl_locus_snvs":len(extracted[gene]),
            "shared_unique_snps":len(matched),
            "match_fraction":len(matched)/len(extracted[gene]) if extracted[gene] else 0,
            "min_shared_snps_pass":int(len(matched)>=args.min_shared_snps),
            "harmonization_direct":sum(r["harmonization"]=="direct" for r in matched),
            "harmonization_swapped":sum(r["harmonization"]=="swapped" for r in matched),
            "anchor_rsid":c["anchor_rsid"],
        })
        if len(matched)<args.min_shared_snps:
            raise SystemExit(f"{gene}: only {len(matched)} shared SNVs (<{args.min_shared_snps})")

    qc_fields=list(qc[0])
    args.output_root.mkdir(parents=True,exist_ok=True)
    with (args.output_root/"STAGE2B_MATCH_QC.tsv").open("w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=qc_fields,delimiter="\t",lineterminator="\n")
        w.writeheader(); w.writerows(qc)

    summary={
        "stage":"CKD Stage 2B locus preparation",
        "builds":{"UKB_PPP":"GRCh38","CKDGen_Stanzick_eGFR":"GRCh37"},
        "matching":"GRCh38 pQTL SNVs lifted to GRCh37; exact chromosome/position + allele-pair match",
        "effect_orientation":"UKB-PPP ALLELE1 is exposure effect allele; outcome beta aligned to exposure ALLELE1",
        "window_bp":args.window_bp,
        "min_info":args.min_info,
        "snv_only":True,
        "outcome_source_rows_scanned":source_rows,
        "n_candidates":len(candidates),
        "extraction":extraction_meta,
        "qc":qc,
        "notes":[
            "Initial coloc excludes indels to avoid cross-build normalization ambiguity.",
            "The next R stage uses coloc.abf as a single-causal-variant screen.",
            "Complex/multi-signal loci, especially HLA-E/MHC, require LD-aware follow-up before interpretation.",
        ],
    }
    (args.output_root/"STAGE2B_PREP_SUMMARY.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"n_candidates":len(candidates),"qc":qc},indent=2))
    print(f"CKD_STAGE2B_PREP_PASS output={args.output_root}")

if __name__=="__main__":
    main()
