#!/usr/bin/env python3
"""Read-only provenance audit for GIGASTROKE EAS AIS allele orientation.

Cross-check GWAS Catalog raw GRCh37 EA/OA records, canonical derived records,
the 11 FASTA-verified alcohol-CS targets, and GRCh37 FASTA itself.
The derivative TSV is a *correctly oriented* diagnostic dataset, never
overwrites either input or changes AIS causal/fine-mapping status.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

RAW_DEFAULT="/srv/is-analysis/data/is/reference/gigastroke/eas/GCST90104545_buildGRCh37.tsv.gz"
CANON_DEFAULT="/srv/is-analysis/data/is/processed/gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz"
FASTA_DEFAULT="/srv/is-analysis/reference/grch37_1000g/human_g1k_v37.fasta"
EXPECTED_SOURCE_MD5="904e12e33ec3c2eeefe4dc93bec4995c"

def indexed_fasta(path):
    idx={}
    with open(str(path)+".fai") as fh:
        for line in fh:
            v=line.rstrip("\n").split("\t")
            if len(v)<5: raise ValueError("Invalid FASTA index")
            idx[v[0]]=(int(v[1]),int(v[2]),int(v[3]),int(v[4]))
    return idx

def ref_at(handle,idx,chrom,pos):
    record=idx.get(str(chrom)) or idx.get("chr"+str(chrom))
    if not record: raise ValueError(f"No FASTA chromosome {chrom}")
    size,offset,per_line,full_line=record
    if pos<1 or pos>size: raise ValueError(f"Invalid reference position {chrom}:{pos}")
    z=pos-1
    handle.seek(offset+(z//per_line)*full_line+z%per_line)
    base=handle.read(1).decode("ascii").upper()
    if base not in "ACGT": raise ValueError(f"Non-ACGT FASTA {chrom}:{pos}")
    return base

def normalized_alleles(ea,oa,ref):
    """Return FASTA-anchored REF/ALT, ALT-effect beta sign and ALT EAF orientation.

    This helper rejects alleles not represented in GRCh37 plus one single alternate.
    """
    ea,oa,ref=ea.upper(),oa.upper(),ref.upper()
    if not all(len(x)==1 and x in "ACGT" for x in (ea,oa,ref)):
        return None
    if ea==oa or ref not in (ea,oa): return None
    alt=oa if ea==ref else ea
    sign=-1 if ea==ref else 1
    return ref,alt,sign

def md5(path):
    h=hashlib.md5()
    with path.open("rb") as f:
        for part in iter(lambda:f.read(8*1024*1024),b""): h.update(part)
    return h.hexdigest()

def awk_scan(path,chrom_col,pos_col,positions,sample_stride=0):
    """One full streaming decompression per file, with compiled awk filtering.

    Positions are validated integer chromosome/position values, not user expressions.
    Returns matched-position rows and a deterministic thinned sample.
    """
    conditions=[f'(${chrom_col}=="{chrom}" && ${pos_col}=="{pos}")'
                for chrom,pos in sorted(positions)]
    clause=" || ".join(conditions) or "0"
    sample_clause=(f" || (NR>1 && NR%{sample_stride}==0)" if sample_stride else "")
    expression=f'NR==1 || ({clause}){sample_clause}'
    # awk uses $ fields: escape of $ above is removed below.
    expression=expression.replace("\\$","$")
    dc=subprocess.Popen(["gzip","-dc",str(path)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    aw=subprocess.Popen(["awk","-F","\t",expression],stdin=dc.stdout,
                        stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    dc.stdout.close()
    result,err=aw.communicate()
    decompress_error=dc.stderr.read()
    dc_rc=dc.wait()
    dc.stderr.close()
    if aw.returncode!=0 or dc_rc!=0:
        raise RuntimeError(f"gzip/awk failed: {dc_rc}/{aw.returncode} {decompress_error[:200]!r} {err[:200]!r}")
    lines=result.decode("utf-8").splitlines()
    if not lines: raise ValueError(f"No header in {path}")
    header=lines[0].split("\t")
    rows=[dict(zip(header,line.split("\t"))) for line in lines[1:]]
    if any(len(line.split("\t"))!=len(header) for line in lines[1:]):
        raise ValueError(f"Malformed TSV row in {path}")
    return rows

def as_number(value):
    try:
        n=float(value)
        return n if math.isfinite(n) else None
    except (ValueError,TypeError):
        return None

def same_numeric(a,b):
    x,y=as_number(a),as_number(b)
    return x is not None and y is not None and math.isclose(x,y,rel_tol=1e-10,abs_tol=1e-12)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--audit-tsv",type=Path,required=True)
    p.add_argument("--raw",type=Path,default=Path(RAW_DEFAULT))
    p.add_argument("--canonical",type=Path,default=Path(CANON_DEFAULT))
    p.add_argument("--fasta",type=Path,default=Path(FASTA_DEFAULT))
    p.add_argument("--metadata-yaml",type=Path,default=Path("/srv/is-analysis/data/is/reference/gigastroke/metadata_phase8/GCST90104545-meta.yaml"))
    p.add_argument("--out-dir",type=Path,required=True)
    p.add_argument("--sample-stride",type=int,default=5000)
    p.add_argument("--expected-md5",default=EXPECTED_SOURCE_MD5)
    args=p.parse_args()
    if args.sample_stride<1: p.error("--sample-stride must be >=1")
    with args.audit_tsv.open(newline="") as f: targets=list(csv.DictReader(f,delimiter="\t"))
    if len(targets)!=11 or len({r["variant_id"] for r in targets})!=11:
        raise ValueError("Expected 11 unique upstream verified variants")
    if any(r["allele_mapping_status"]!="VERIFIED_ALLELE_SET" for r in targets):
        raise ValueError("VEP upstream allele audit not complete")
    metadata=args.metadata_yaml.read_text()
    if "genome_assembly: GRCh37" not in metadata or "sample_size: 256274" not in metadata:
        raise ValueError("GIGASTROKE metadata version or assembly unexpected")
    actual_md5=md5(args.raw)
    if actual_md5.lower()!=args.expected_md5.lower():
        raise ValueError(f"Raw GIGASTROKE checksum mismatch {actual_md5}")
    positions={(str(x["chrom"]),int(x["pos"])) for x in targets}
    raw_rows=awk_scan(args.raw,1,2,positions)
    canon_rows=awk_scan(args.canonical,5,6,positions,args.sample_stride)
    idx=indexed_fasta(args.fasta)
    with args.fasta.open("rb") as fa:
        raw_by={}
        for r in raw_rows:
            key=(str(r["chromosome"]),int(r["base_pair_location"]))
            if key in positions: raw_by.setdefault(key,[]).append(r)
        canon_by={}
        sampling=Counter()
        sample_total=0
        # Target rows might overlap sampling: exclude all targets from sample.
        for r in canon_rows:
            key=(str(r["chr"]),int(r["pos"]))
            if key in positions:
                canon_by.setdefault(key,[]).append(r)
                continue
            ea,oa=r["effect_allele"].upper(),r["other_allele"].upper()
            try: ref=ref_at(fa,idx,*key)
            except ValueError:
                sampling["REFERENCE_LOOKUP_FAILED"]+=1
                continue
            normalized=normalized_alleles(ea,oa,ref)
            if normalized is None: sampling["ALLELE_PAIR_NOT_REF_RESOLVED"]+=1;continue
            sample_total+=1
            vid=r["variant_id"].split(":")
            if len(vid)!=4: sampling["MALFORMED_VARIANT_ID"]+=1;continue
            if vid[2].upper()==ref: sampling["CANONICAL_ID_REF_EQUALS_FASTA"]+=1
            else: sampling["CANONICAL_ID_REF_NOT_FASTA"]+=1
            if vid[2].upper()==oa and vid[3].upper()==ea:
                sampling["CANONICAL_ID_EQUALS_OTHER_EFFECT"]+=1
            else:
                sampling["CANONICAL_ID_NOT_OTHER_EFFECT"]+=1
        result=[]
        for t in targets:
            chrom,pos=str(t["chrom"]),int(t["pos"])
            key=(chrom,pos)
            rr=raw_by.get(key,[])
            cc=canon_by.get(key,[])
            ref=ref_at(fa,idx,chrom,pos)
            if ref!=t["ref"].upper(): raise ValueError(f"FASTA REF drift {key}")
            common={"locus":t["locus"],"rsid":t["allele_verified_rsids"],
              "variant_id_fasta":t["variant_id"],"source":"GCST90104545_EAS_AIS",
              "raw_rows_at_position":len(rr),"canonical_rows_at_position":len(cc)}
            if not rr and not cc:
                result.append({**common,"status":"MISSING_BOTH","source_md5_pass":True})
                continue
            if len(rr)!=1 or len(cc)!=1:
                result.append({**common,"status":"MULTIPLE_OR_DISCORDANT_RECORDS","source_md5_pass":True})
                continue
            raw=rr[0]; can=cc[0]
            ea,oa=raw["effect_allele"].upper(),raw["other_allele"].upper()
            norm=normalized_alleles(ea,oa,ref)
            if norm is None or {ea,oa}!={t["ref"].upper(),t["alt"].upper()}:
                result.append({**common,"status":"RAW_ALLELE_PAIR_INVALID","source_md5_pass":True})
                continue
            _,alt,sign=norm
            if alt!=t["alt"].upper(): raise ValueError(f"ALT mismatch {t['variant_id']}")
            comparisons={
                "ea":can["effect_allele"].upper()==ea,
                "oa":can["other_allele"].upper()==oa,
                "beta":same_numeric(raw["beta"],can["beta"]),
                "se":same_numeric(raw["standard_error"],can["se"]),
                "p":same_numeric(raw["p_value"],can["p"]),
                "eaf":same_numeric(raw["effect_allele_frequency"],can["eaf"]),
            }
            stats_ok=all(comparisons.values())
            tokens=can["variant_id"].split(":")
            expected=f"{chrom}:{pos}:{ref}:{alt}"
            can_ref_ok=len(tokens)==4 and tokens[2].upper()==ref and tokens[3].upper()==alt
            canonical_id_is_oa_ea=can["variant_id"]==f"{chrom}:{pos}:{oa}:{ea}"
            if not stats_ok: status="RAW_CANONICAL_STAT_DISCORDANCE"
            elif can_ref_ok: status="PASS_FASTA_REF_ALT"
            elif canonical_id_is_oa_ea: status="CANONICAL_ID_IS_OTHER_EFFECT_NOT_REF_ALT"
            else: status="CANONICAL_ID_UNEXPECTED"
            b=as_number(raw["beta"])
            af=as_number(raw["effect_allele_frequency"])
            result.append({**common,"status":status,"source_md5_pass":True,
              "original_canonical_id":can["variant_id"],
              "validated_fasta_variant_id":expected,
              "fasta_ref":ref,"fasta_alt":alt,"raw_effect_allele":ea,
              "raw_other_allele":oa,"effect_allele_is_ref":ea==ref,
              "canonical_id_equals_other_effect":canonical_id_is_oa_ea,
              "source_canonical_stats_equal":stats_ok,
              "source_canonical_equality_checks":json.dumps(comparisons,sort_keys=True),
              "ais_beta_alt":b*sign if b is not None else "",
              "ais_se":raw["standard_error"],"ais_p":raw["p_value"],
              "ais_eaf_alt":(af if sign==1 else 1-af) if af is not None else "",
              "ais_beta_source_effect":raw["beta"],
              "allele_orientation":"EFFECT_ALT" if sign==1 else "EFFECT_REF_FLIPPED"})
    critical=[r for r in result if r["status"] not in
              {"MISSING_BOTH","PASS_FASTA_REF_ALT","CANONICAL_ID_IS_OTHER_EFFECT_NOT_REF_ALT"}]
    if critical: raise ValueError("Critical source/reference discordance; no outputs created: "+
                                  json.dumps(critical[:3]))
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise FileExistsError("Output directory already has files")
    args.out_dir.mkdir(parents=True,exist_ok=True)
    output=args.out_dir/"GIGASTROKE_AIS_11_SNP_FASTA_PROVENANCE.tsv"
    fields=list(dict.fromkeys(k for row in result for k in row.keys()))
    with output.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
        w.writeheader();w.writerows(result)
    summary={"analysis":"GIGASTROKE_EAS_GRCH37_REF_ALT_AUDIT",
      "status":"DIAGNOSTIC_READONLY_NO_CANONICAL_MUTATION",
      "timestamp_utc":datetime.now(timezone.utc).isoformat(),
      "original_source_file":str(args.raw),
      "original_source_md5":actual_md5,
      "original_source_md5_expected":args.expected_md5,
      "original_source_metadata":str(args.metadata_yaml),
      "genome_build":"GRCh37",
      "canonical":str(args.canonical),
      "fasta":str(args.fasta),
      "n_variants":len(result),"n_sampled_reference_resolved":sample_total,
      "sample_stride":args.sample_stride,
      "variant_counts":dict(Counter(r["status"] for r in result)),
      "thinned_canonical_sample_counts":dict(sampling),
      "evidence_scope":"Raw-to-derived 11 loci + deterministic thinned genome-wide FASTA audit; not all rows",
      "official_alcohol_GWAS_genome_build":"hg19_GRCh37_JCGE_OFFICIAL_SITE_VERIFIED",
      "study_overlap":"BBJ_DISCOVERY_OVERLAP_POSSIBLE_NOT_INDEPENDENT_REPLICATION",
      "fine_mapping_ALDH2":"BLOCKED","fine_mapping_ADH1B":"EXPLORATORY"}
    (args.out_dir/"GIGASTROKE_AIS_REFALT_AUDIT_MANIFEST.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("GIGASTROKE_REFALT_PROVENANCE_PASS",json.dumps(summary["variant_counts"],sort_keys=True),flush=True)
    print("CANONICAL_REFERENCE_SAMPLE",json.dumps(summary["thinned_canonical_sample_counts"],sort_keys=True),flush=True)
    for r in result:
        print("VARIANT",r["rsid"],r["status"],r.get("original_canonical_id",""),
              r.get("validated_fasta_variant_id",""),flush=True)

if __name__=="__main__":
    main()
