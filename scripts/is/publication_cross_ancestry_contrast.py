#!/usr/bin/env python3
"""Frozen 24 EAS-region markers: conditional EAS/EUR AIS reported-effect contrasts.

This is retrospective same-release summary-statistic comparison, not a novel
ancestry-specific mechanism, independent replication, or proof of comparable
log-odds units. Source model differences and winner's curse remain.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, math
from collections import defaultdict
from pathlib import Path
from scripts.is import audit_gws_eas_24_leads_original as old
from scripts.is import p3a_publisher_md5_gate as gate

def digest(path, algo):
    with path.open("rb") as f:
        return hashlib.file_digest(f, algo).hexdigest()

def write(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+"\n")
    print(json.dumps(obj, indent=2, allow_nan=False))

def extract(ancestry, root, legacy, out):
    spec=gate.SOURCES[ancestry]
    data=root/"source.gz"
    fields=gate.publisher_metadata((root/"source.yaml").read_text(encoding="utf-8-sig"), spec)
    if (data.stat().st_size, digest(data,"md5"), digest(data,"sha256")) != (spec["bytes"],spec["publisher_md5_pinned"],spec["full_sha256_pinned"]):
        raise ValueError("Full original publisher file integrity failed")
    panel=old.read_legacy_180(legacy)
    wanted={(x["chromosome"],x["position"]):x for x in panel["selected"]}
    hits=defaultdict(list)
    count=0
    with gzip.open(data,"rt",encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        h=reader.fieldnames or []
        if len(set(h))!=len(h) or not set(old.HEADER)<=set(h):
            raise ValueError("Original schema mismatch "+repr(h))
        for row in reader:
            count+=1
            if count>9_000_000:raise ValueError("Unbounded rows")
            chrom=str(row["chromosome"]).strip().upper().removeprefix("CHR")
            try: key=chrom,int(row["base_pair_location"])
            except (ValueError,TypeError):continue
            if key in wanted:
                hits[key].append(old.inspect_original_row(row))
    if count!=spec["original_rows_in_prior_source_qc"]:
        raise ValueError("Full source row count changed "+str(count))
    result=[]
    for a in panel["selected"]:
        h=hits[a["chromosome"],a["position"]]
        r=dict(region_id=a["region_id"],variant=a["GRCh37_variant_key"],alt=a["historical_alt_allele"],ref=a["historical_ref_allele"],source_position_rows=len(h),status="HOLD")
        if len(h)==1 and h[0] is not None:
            b=h[0]
            if {b["effect_allele"],b["other_allele"]}=={r["alt"],r["ref"]}:
                sign=1 if b["effect_allele"]==r["alt"] else -1
                r.update(status="MATCHED",beta_alt=sign*b["beta"],se=b["se"],p=b["p"],alt_eaf=b["eaf"] if sign==1 else 1-b["eaf"])
                if ancestry=="EAS" and not old.compare(a,h)["joint_effect_p_eaf_reproduced"]:
                    raise ValueError("EAS archival effect provenance differs: "+a["region_id"])
        result.append(r)
    if ancestry=="EAS" and any(r["status"]!="MATCHED" for r in result):
        raise ValueError("EAS previously verified anchor not recovered")
    write(out,dict(status="EXACT_ORIGINAL_SOURCE_EXTRACT",ancestry=ancestry,accession=spec["accession"],source_sha256=spec["full_sha256_pinned"],source_md5=fields["data_file_md5sum"],source_rows=count,source_bytes=spec["bytes"],family=len(result),rows=result,withheld=panel["withheld"],missing=panel["missing"]))

def contrast(a,b,n=24):
    delta=a["beta_alt"]-b["beta_alt"]
    sa,sb=a["se"],b["se"]
    if not all(math.isfinite(x) for x in (delta,sa,sb)) or min(sa,sb)<=0:
        raise ValueError("Invalid contrast inputs")
    se0=math.hypot(sa,sb)
    sew=sa+sb
    p0=math.erfc(abs(delta)/se0/math.sqrt(2))
    pw=math.erfc(abs(delta)/sew/math.sqrt(2))
    return dict(delta_EAS_minus_EUR=delta,se_zero_covariance=se0,se_maximum_variance=sew,
                p_zero_covariance=p0,p_maximum_variance=pw,
                bonferroni_p_zero_covariance=min(1,n*p0),bonferroni_p_maximum_variance=min(1,n*pw),
                ci95_maximum_variance=[delta-1.959964*sew,delta+1.959964*sew])

def combine(root,out):
    inputs={anc:json.loads((root/(anc+".json")).read_text()) for anc in ("EAS","EUR")}
    for anc,s in inputs.items():
        if s["ancestry"]!=anc or s["source_sha256"]!=gate.SOURCES[anc]["full_sha256_pinned"] or s["family"]!=24:
            raise ValueError("Source identity drift")
    d={anc:{r["region_id"]:r for r in s["rows"]} for anc,s in inputs.items()}
    if len(d["EAS"])!=24 or set(d["EAS"])!=set(d["EUR"]):
        raise ValueError("Fixed family changed")
    results=[]
    for region,a in sorted(d["EAS"].items()):
        b=d["EUR"][region]
        if a["variant"]!=b["variant"] or a["alt"]!=b["alt"]:raise ValueError("Anchor changed")
        r=dict(region_id=region,variant=a["variant"],EAS=a,EUR=b,status="HOLD_UNMATCHED")
        if a["status"]==b["status"]=="MATCHED":
            r.update(status="CONDITIONAL_NUMERIC_CONTRAST",same_direction=a["beta_alt"]*b["beta_alt"]>0,
                     abs_alt_eaf_difference=abs(a["alt_eaf"]-b["alt_eaf"]),**contrast(a,b))
        results.append(r)
    matched=[r for r in results if r["status"]=="CONDITIONAL_NUMERIC_CONTRAST"]
    write(out,dict(status="EXPLORATORY_REPORTED_BETA_CONTRAST_NOT_CAUSAL",family=24,matched=len(matched),
        same_direction=sum(r["same_direction"] for r in matched),
        bonferroni_zero_covariance=sum(r["bonferroni_p_zero_covariance"]<.05 for r in matched),
        bonferroni_maximum_variance=sum(r["bonferroni_p_maximum_variance"]<.05 for r in matched),
        eaf_difference_ge_0_2=sum(r["abs_alt_eaf_difference"]>=.2 for r in matched),
        source_accessions={anc:s["accession"] for anc,s in inputs.items()},results=results,
        withheld=inputs["EAS"]["withheld"],missing=inputs["EAS"]["missing"],
        caveats=["Historical EAS-selected leads may have winner's curse; 24 contrasts are exploratory and not independent loci by demonstrated LD.",
                 "Conditional on comparable effect units; GIGASTROKE combines logistic and some Cox estimates; exact file-wide effect-model equivalence not established.",
                 "Maximum variance bound uses Cauchy-Schwarz on estimator covariance, not a correction for selection, scale, population structure or differing LD.",
                 "Missing EUR variants not null effects. Same marker may tag different causal variants in different ancestries.",
                 "No direct EAS cell eQTL, signed matched QTL LD, independent replication or molecular colocalization."]))

def selftest():
    a={"beta_alt":.10,"se":.01};b={"beta_alt":.02,"se":.02}
    c=contrast(a,b)
    assert math.isclose(c["delta_EAS_minus_EUR"],.08)
    assert math.isclose(c["se_maximum_variance"],.03)
    assert c["p_maximum_variance"]>=c["p_zero_covariance"]
    assert c["bonferroni_p_maximum_variance"]==min(1,24*c["p_maximum_variance"])
    assert contrast(a,a)["p_zero_covariance"]==1
    assert old.compare({"historical_ref_allele":"A","historical_alt_allele":"G","historical_beta_ALT":.1,"historical_p":.01,"historical_ALT_EAF":.2,"region_id":"test","GRCh37_variant_key":"1:1:A:G","chromosome":"1"},
       [{"effect_allele":"A","other_allele":"G","beta":-.1,"eaf":.8,"se":.03,"p":.01}])["joint_effect_p_eaf_reproduced"]
    print("Six focused conditional-contrast and orientation assertions passed")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("mode",choices=["selftest","extract","combine"])
    p.add_argument("--ancestry",choices=["EAS","EUR"])
    p.add_argument("--root",type=Path,default=Path("."))
    p.add_argument("--legacy",type=Path)
    p.add_argument("--out",type=Path,default=Path("result.json"))
    a=p.parse_args()
    if a.mode=="selftest":selftest()
    elif a.mode=="extract":extract(a.ancestry,a.root,a.legacy,a.out)
    else:combine(a.root,a.out)
if __name__=="__main__":main()
