#!/usr/bin/env python3
"""Seven GWS-group GWAS/reference overlap and effect-allele audit (GRCh37).

Read-only GWAS inputs; output is explicitly a *reference allele audit*,
not a finemapping-ready harmonized meta-analysis. No strand flips inferred.
"""
import argparse
import csv
import gzip
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
DATA=Path("/srv/is-analysis/data/is/processed")
REF_NEW=Path("/srv/is-analysis/data/is/ld_reference/broad_v1/1kg_eas")
REF_OLD=Path("/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/pgen_gt")
ANCHORS={"IS_XDATA_G0007":"L001","IS_XDATA_G0015":"L002",
         "IS_XDATA_G0022":"L003","IS_XDATA_G0023":"L004"}
CLASSES=("MATCH_ALT_EFFECT","MATCH_REF_EFFECT","PALINDROMIC_REVIEW",
         "SOURCE_REF_ALT_CONFLICT","REFERENCE_ALLELES_DIFFER","NO_REFERENCE_POSITION",
         "UNSUPPORTED_SOURCE_ALLELE","MALFORMED_EFFECT")

def fnum(value):
    try:
        v=float(value)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None

def pvar_path(group,ref_new,ref_old):
    expanded = ref_new/f"{group}.stable.pvar"
    expanded_prefix = ref_new/f"{group}.stable"
    triplet = [Path(str(expanded_prefix)+ext) for ext in (".pgen",".pvar",".psam")]
    if all(p.is_file() for p in triplet):
        p = expanded
    elif any(p.is_file() for p in triplet):
        raise RuntimeError("Incomplete expanded PGEN triplet: "+group)
    elif group in ANCHORS:
        p = ref_old/f"BBJ_IS_{ANCHORS[group]}.1KG_EAS.GRCh37.pvar"
    else:
        p = expanded
    if not p.is_file():raise FileNotFoundError("Missing genotype PVAR: "+str(p))
    return p

def ref_index(groups,ref_new,ref_old):
    indexes={}
    for group in groups:
        ref=defaultdict(list)
        with pvar_path(group,ref_new,ref_old).open() as f:
            for line in f:
                if line.startswith("#"):continue
                a=line.rstrip("\n").split("\t")
                if len(a)<5:continue
                chrom,pos,vid,ral,aal=a[:5]
                if vid=="." or ":" not in vid:
                    raise ValueError("Reference variant IDs not stable for "+group)
                if len(ral)!=1 or len(aal)!=1 or ral not in "ACGT" or aal not in "ACGT":
                    continue
                ref[(chrom.removeprefix("chr"),int(pos))].append((ral,aal,vid))
        indexes[group]=ref
    return indexes

def datasets(data):
    files=[data/"japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz"]
    files+=sorted((data/"gigastroke/eas").glob("*_GRCh37.canonical.tsv.gz"))
    if len(files)!=6 or not all(p.is_file() for p in files):
        raise RuntimeError("Expected BBJ + 5 GIGASTROKE EAS canonical inputs; found "+str(files))
    return files

def allele_class(r,a,ea,oa,source_ref,source_alt,eaf):
    ref,alt,_=r
    if source_ref and source_alt and (source_ref!=ref or source_alt!=alt):
        return "SOURCE_REF_ALT_CONFLICT",None,None
    if set([ea,oa]) != set([ref,alt]):
        return "REFERENCE_ALLELES_DIFFER",None,None
    if set([ref,alt]) in ({"A","T"},{"C","G"}) and (eaf is None or .4<=eaf<=.6):
        return "PALINDROMIC_REVIEW",None,None
    if ea==alt:
        return "MATCH_ALT_EFFECT",1,eaf
    return "MATCH_REF_EFFECT",-1,1-eaf if eaf is not None else None

def process(root,data,ref_new,ref_old):
    with (root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv").open() as f:
        all_groups=list(csv.DictReader(f,delimiter="\t"))
    groups=[g for g in all_groups if g["has_gws"]=="1"]
    if len(groups)!=7:raise ValueError("Expected exactly seven GWS groups")
    indexes=ref_index([g["group_id"] for g in groups],ref_new,ref_old)
    bychr=defaultdict(list)
    for g in groups:
        bychr[g["chr"]].append(g)
    results=root/"IS_GWS_GWAS_REFERENCE_OVERLAP.tsv"
    records=root/"IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz"
    summary=defaultdict(Counter)
    scanned=Counter()
    temp=Path(str(records)+".part")
    fields=["dataset","phenotype","group_id","chr","pos","variant_id",
            "ref","alt","effect_allele","other_allele","beta","se","p","eaf",
            "alt_effect_beta","alt_effect_eaf","qc_status","source_variant_id",
            "source_ref","source_alt","ancestry","build"]
    try:
        with gzip.open(temp,"wt",newline="") as o:
            w=csv.DictWriter(o,fieldnames=fields,delimiter="\t");w.writeheader()
            for file in datasets(data):
                with gzip.open(file,"rt",newline="") as h:
                    src=csv.DictReader(h,delimiter="\t")
                    names=set(src.fieldnames or [])
                    if not {"dataset","phenotype","chr","pos","effect_allele",
                            "other_allele","beta","se","p","eaf","build"}.issubset(names):
                        raise ValueError("Invalid input schema "+str(file))
                    for r in src:
                        ds=r["dataset"]; scanned[ds]+=1
                        chrom=r["chr"].removeprefix("chr")
                        relevant=bychr.get(chrom,[])
                        if not relevant:continue
                        try:pos=int(r["pos"])
                        except (TypeError,ValueError):continue
                        for group in relevant:
                            if pos<int(group["window_start"]) or pos>int(group["window_end"]):
                                continue
                            tag=group["group_id"]; counter=summary[(ds,tag)]
                            counter["window_variants"]+=1
                            ref_candidates=indexes[tag].get((chrom,pos),[])
                            if not ref_candidates:
                                counter["NO_REFERENCE_POSITION"]+=1;continue
                            ea=r["effect_allele"].upper();oa=r["other_allele"].upper()
                            if not ea or not oa or len(ea)!=1 or len(oa)!=1 or ea==oa or ea not in "ACGT" or oa not in "ACGT":
                                counter["UNSUPPORTED_SOURCE_ALLELE"]+=1;continue
                            eaf=fnum(r.get("eaf"))
                            beta=fnum(r.get("beta"));se=fnum(r.get("se"));p=fnum(r.get("p"))
                            if any(x is None for x in (beta,se,p)) or se<=0 or not 0<=p<=1:
                                counter["MALFORMED_EFFECT"]+=1;continue
                            sr=r.get("ref","").upper();sa=r.get("alt","").upper()
                            matches=[v for v in ref_candidates if set([v[0],v[1]])==set([ea,oa])]
                            if not matches:
                                counter["REFERENCE_ALLELES_DIFFER"]+=1;continue
                            if len(matches)!=1:
                                counter["DUPLICATE_MATCHED_ALLELE_PAIR"]+=1;continue
                            ref,alt,vid=matches[0]
                            kind,direction,alt_eaf=allele_class(matches[0],None,ea,oa,sr,sa,eaf)
                            counter[kind]+=1
                            if kind not in ("MATCH_ALT_EFFECT","MATCH_REF_EFFECT","PALINDROMIC_REVIEW"):
                                continue
                            w.writerow({"dataset":ds,"phenotype":r["phenotype"],"group_id":tag,
                                "chr":chrom,"pos":pos,"variant_id":vid,"ref":ref,"alt":alt,
                                "effect_allele":ea,"other_allele":oa,"beta":beta,"se":se,
                                "p":p,"eaf":r.get("eaf",""),
                                "alt_effect_beta":beta*direction if direction is not None else "",
                                "alt_effect_eaf":alt_eaf if alt_eaf is not None else "",
                                "qc_status":kind,"source_variant_id":r.get("variant_id",""),
                                "source_ref":sr,"source_alt":sa,
                                "ancestry":r.get("ancestry",""),"build":r["build"]})
                            counter["exported_matching_pairs"]+=1
        temp.replace(records)
    except BaseException:
        if temp.exists():temp.unlink()
        raise
    scores=[]
    for (ds,group),c in sorted(summary.items()):
        window=c["window_variants"]
        high=c["MATCH_ALT_EFFECT"]+c["MATCH_REF_EFFECT"]
        scores.append({"dataset":ds,"group_id":group,"window_variants":window,
            "matched_oriented_nonpal":high,
            "palindromic_to_review":c["PALINDROMIC_REVIEW"],
            "no_reference_position":c["NO_REFERENCE_POSITION"],
            "mismatch_alleles":c["REFERENCE_ALLELES_DIFFER"],
            "bbj_ref_alt_conflict":c["SOURCE_REF_ALT_CONFLICT"],
            "unusable_alleles":c["UNSUPPORTED_SOURCE_ALLELE"],
            "malformed_effect":c["MALFORMED_EFFECT"],
            "exported_matching_pairs":c["exported_matching_pairs"],
            "nonpal_orientation_fraction":round(high/window,6) if window else 0,
            "readiness":"REFERENCE_MATCH_AUDITED_NOT_FINE_MAPPING"})
    with results.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(scores[0]),delimiter="\t")
        w.writeheader();w.writerows(scores)
    agg={"n_input_datasets":len(scanned),"source_scan_rows":dict(scanned),
         "total_windows_with_source_data":len(scores),
         "overall_window_variants":sum(x["window_variants"] for x in scores),
         "nonpal_oriented":sum(x["matched_oriented_nonpal"] for x in scores),
         "palindromic_review":sum(x["palindromic_to_review"] for x in scores),
         "source_ref_conflicts":sum(x["bbj_ref_alt_conflict"] for x in scores),
         "reference_unmatched_positions":sum(x["no_reference_position"] for x in scores),
         "note":"Allele matching and effect orientation only. No LD independency or joint credible sets established"}
    (root/"IS_GWS_GWAS_REFERENCE_OVERLAP_SUMMARY.json").write_text(json.dumps(agg,indent=2))
    print(json.dumps(agg,indent=2))
    return agg

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,default=ROOT)
    ap.add_argument("--data",type=Path,default=DATA)
    ap.add_argument("--ref-new",type=Path,default=REF_NEW)
    ap.add_argument("--ref-old",type=Path,default=REF_OLD)
    a=ap.parse_args()
    process(a.root,a.data,a.ref_new,a.ref_old)
