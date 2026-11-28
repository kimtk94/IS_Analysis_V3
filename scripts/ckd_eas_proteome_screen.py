#!/usr/bin/env python3
"""Exploratory assay-wide CKB cis-sentinel screen; shared source/QC methods adapted from #338."""
import argparse
import collections
import gzip
import hashlib
import json
import math
import re
import csv
from pathlib import Path

EXPOSURE_SOURCE = "c6a6f4cfa58a2ee886a86482d7e08913c5c2aed4e42c3aa07f17c625d4fd6764"
SOURCES = {"egfr":(42774049,111041990,"TWB2_BBJ_eGFR_hg19_METAL_FUMA_noNA.gz"),
           "bun":(43218600,132381963,"TWB2_BBJ_BUN_hg19_METAL_FUMA_noNA.gz")}
COMPLEMENT = str.maketrans("ACGT","TGCA")

def digest(path, algorithm="sha256"):
    with path.open("rb") as fh:
        return hashlib.file_digest(fh,algorithm).hexdigest()

def orient(ea,oa,out_ea,out_oa):
    if {ea,oa} in ({"A","T"},{"C","G"}):
        return None
    if (ea,oa)==(out_ea,out_oa):
        return 1
    if (ea,oa)==(out_oa,out_ea):
        return -1
    return None

def wald(bx,sx,by,sy):
    if not all(math.isfinite(x) for x in (bx,sx,by,sy)) or bx==0 or min(sx,sy)<=0:
        raise ValueError("Invalid effect estimates")
    b=by/bx
    se=math.sqrt((sy/bx)**2+(by*sx/bx**2)**2)
    p=math.erfc(abs(b)/se/math.sqrt(2))
    return dict(mr_beta=b,mr_se_delta_independence=se,mr_ci95=[b-1.959964*se,b+1.959964*se],mr_p=p,
                mr_p_first_order=math.erfc(abs(by/sy)/math.sqrt(2)),
                mr_se_worst_covariance=abs(sy/bx)+abs(by*sx/bx**2))

def checked_chain(root,name,manifest):
    path=root/name
    lines=(root/manifest).read_text().splitlines()
    md5s=[s.split()[0] for s in lines if len(s.split())==2 and s.split()[1].lstrip("*")==name]
    if len(md5s)!=1 or digest(path,"md5")!=md5s[0].lower():
        raise ValueError("UCSC chain MD5 mismatch "+name)
    return path

def lift(rows,forward,reverse):
    results=[]
    for r in rows:
        m=re.fullmatch(r"chr(\d+):(\d+):([ACGT]):([ACGT])",r["sentinel"])
        out=dict(r)
        if not m:
            results.append({**out,"status":"HOLD_NON_SNV"})
            continue
        chrom,pos,ref,alt=m.groups()
        pos=int(pos)
        ea,oa=r["effect_allele"],r["other_allele"]
        if {ea,oa}!={ref,alt}:
            raise ValueError("Exposure allele pair mismatch")
        if {ea,oa} in ({"A","T"},{"C","G"}):
            results.append({**out,"status":"HOLD_PALINDROMIC_NO_OUTCOME_EAF"})
            continue
        hits=forward.convert_coordinate("chr"+chrom,pos-1)
        if not hits or len(hits)!=1:
            results.append({**out,"status":"HOLD_NONUNIQUE_LIFT"})
            continue
        c,p,s,_=hits[0]
        back=reverse.convert_coordinate(c,p)
        if c!="chr"+chrom or not back or len(back)!=1 or back[0][:3]!=("chr"+chrom,pos-1,s):
            results.append({**out,"status":"HOLD_RECIPROCAL_LIFT"})
            continue
        if s not in ("+","-"):
            raise ValueError("Invalid chain strand")
        if s=="-":
            ea,oa,ref,alt=(a.translate(COMPLEMENT) for a in (ea,oa,ref,alt))
        results.append({**out,"chr37":int(c[3:]),"pos37":p+1,"strand":s,
                        "ea37":ea,"oa37":oa,"ref37":ref,"alt37":alt,"status":"MAPPED"})
    return results

def source_meta(root,trait):
    fid,size,name=SOURCES[trait]
    meta=json.loads((root/"figshare.json").read_text())
    if meta.get("id")!=24356587:
        raise ValueError("Unexpected Figshare article identity")
    match=[r for r in meta["files"] if r["id"]==fid and r["name"]==name and r["size"]==size]
    if len(match)!=1:
        raise ValueError("Exact source file metadata not verified")
    md5=match[0].get("computed_md5") or match[0].get("supplied_md5")
    path=root/"outcome.gz"
    if not re.fullmatch(r"[a-fA-F0-9]{32}",str(md5)) or path.stat().st_size!=size or digest(path,"md5")!=md5.lower():
        raise ValueError("Outcome file size or publisher MD5 mismatch")
    return dict(file_id=fid,name=name,bytes=size,md5=md5,sha256=digest(path),article=24356587)

def scan(path,rows,families):
    targets={(r["chr37"],r["pos37"]) for r in rows if r["status"]=="MAPPED"}
    found=collections.defaultdict(list)
    malformed=0
    total=0
    with gzip.open(path,"rt",encoding="utf-8") as f:
        h=f.readline().strip().split()
        required=("CHR","BP","SNP","Allele1","Allele2","Effect","StdErr","P")
        if not set(required)<=set(h) or len(h)!=len(set(h)):
            raise ValueError("Outcome header mismatch: "+repr(h))
        ix={k:h.index(k) for k in required}
        for line in f:
            if not line.strip():
                continue
            total+=1
            parts=line.split()
            if len(parts)!=len(h):
                malformed+=1
                continue
            try:
                coord=(int(parts[ix["CHR"]]),int(parts[ix["BP"]]))
            except ValueError:
                malformed+=1
                continue
            if coord in targets:
                found[coord].append({k:parts[i] for k,i in ix.items()})
    if total<1000000 or malformed:
        raise ValueError("Incomplete or malformed outcome: "+str((total,malformed)))
    for r in rows:
        if r["status"]!="MAPPED":
            continue
        rr=found[r["chr37"],r["pos37"]]
        matches=[]
        for q in rr:
            sign=orient(r["ea37"],r["oa37"],q["Allele1"].upper(),q["Allele2"].upper())
            if sign is not None:
                matches.append((q,sign))
        if len(matches)!=1:
            r["status"]="HOLD_MISSING_OR_AMBIGUOUS_OUTCOME"
            r["outcome_coordinate_rows"]=len(rr)
            continue
        q,sign=matches[0]
        try:
            by,sy,p=map(float,(q["Effect"],q["StdErr"],q["P"]))
            if not all(math.isfinite(x) for x in (by,sy,p)) or sy<=0 or not 0<=p<=1:
                raise ValueError("Invalid outcome")
            result = wald(float(r["beta"]),float(r["se"]),sign*by,sy)
            if set(result) & set(r):
                raise ValueError("Derived field collision with source evidence")
            r.update(result)
        except (ValueError,OverflowError):
            r["status"]="HOLD_OUTCOME_NUMERIC"
            continue
        family=families[r["platform"]]
        r.update(status="EXPLORATORY_WALD",rsid=q["SNP"],outcome_beta_aligned=sign*by,
                 outcome_se=sy,outcome_source_p=p,bonferroni_family=family,
                 bonferroni_p=min(1,r["mr_p"]*family))
    return total,malformed,h,rows


MAPPING_SHA = "b1180c9f633b378b9b9793b3093c702aa393172bcabe28e83fcc3755eaff4bd0"
PLATFORMS = ("Olink","SomaScan")
EXPECTED_ROWS = {"Olink":3552,"SomaScan":2863}
def assay_anchors(path):
    from openpyxl import load_workbook
    if digest(path)!=EXPOSURE_SOURCE:
        raise ValueError("Exposure workbook digest changed")
    w=load_workbook(path,read_only=True,data_only=True)
    selected=[]
    held=[]
    families={}
    counts={}
    required={"sentinel_variant","CHR","POS","effect_allele","non_effect_allele",
              "effect_allele_frequency","INFO","BETA","SE","P","cis_trans"}
    for platform,aid in (("Olink","olink_id"),("SomaScan","somascan_id")):
        rows=w[platform].iter_rows(values_only=True)
        h=[str(v) for v in next(rows)]
        if not required|{aid}<=set(h) or len(set(h))!=len(h):
            raise ValueError("Unexpected exposure header")
        groups=collections.defaultdict(dict)
        conflicts=set()
        total=0
        cis=0
        rejected=collections.Counter()
        for row in rows:
            total+=1
            r=dict(zip(h,row))
            if str(r["cis_trans"]).strip().lower()!="cis":
                continue
            cis+=1
            assay=str(r[aid] or "").strip()
            if not assay:
                raise ValueError("Missing assay ID in cis source")
            variants=groups[assay]
            sentinel=str(r["sentinel_variant"]).strip()
            ea,oa=str(r["effect_allele"]).upper(),str(r["non_effect_allele"]).upper()
            try:
                m=re.fullmatch(r"chr(\d+):(\d+):([ACGT]+):([ACGT]+)",sentinel)
                b,se,eaf,info,p=map(float,(r["BETA"],r["SE"],r["effect_allele_frequency"],r["INFO"],r["P"]))
                if not (m and int(m[1])==int(r["CHR"]) and int(m[2])==int(r["POS"]) and
                        1<=int(m[1])<=22 and {m[3],m[4]}=={ea,oa} and ea!=oa and
                        all(math.isfinite(x) for x in (b,se,eaf,info,p)) and se>0 and 0<eaf<1 and
                        0.3<info<=1 and 0<=p<5e-8):
                    raise ValueError("Eligibility")
                f=(b/se)**2
                if not math.isfinite(f) or f<10:
                    raise ValueError("Weak instrument")
            except (TypeError,ValueError,OverflowError):
                rejected[assay]+=1
                continue
            record=dict(platform=platform,assay_id=assay,sentinel=sentinel,effect_allele=ea,
                        other_allele=oa,beta=b,se=se,eaf=eaf,info=info,source_p=p,F=f)
            if sentinel in variants and variants[sentinel]!=record:
                conflicts.add(assay)
            variants[sentinel]=record
        if total!=EXPECTED_ROWS[platform]:
            raise ValueError("Pinned exposure row count differs")
        families[platform]=len(groups)
        for assay,variants in sorted(groups.items()):
            if assay in conflicts or not variants:
                held.append(dict(platform=platform,assay_id=assay,
                    status="HOLD_CONFLICTING_EXPOSURE" if assay in conflicts else "HOLD_NO_ELIGIBLE_CIS"))
                continue
            # Freeze strongest F before outcome access. Do NOT swap if absent/palindromic.
            record=sorted(variants.values(),key=lambda r:(-r["F"],r["sentinel"]))[0].copy()
            record["eligible_sentinels"]=len(variants)
            selected.append(record)
        counts[platform]={"source_rows":total,"cis_rows":cis,"cis_assay_family":len(groups),
                          "assays_with_invalid_rows":len(rejected),"conflicting_assays":len(conflicts)}
    w.close()
    if any(n<100 for n in families.values()):
        raise ValueError("Unexpectedly small family")
    return selected,held,families,counts

def mapping(path):
    from openpyxl import load_workbook
    if digest(path)!=MAPPING_SHA:
        raise ValueError("Published reagent mapping digest changed")
    w=load_workbook(path,read_only=True,data_only=True)
    rows=w["Supplementary_Data_1"].iter_rows(values_only=True)
    h=[str(v) for v in next(rows)]
    req={"uniprot_id","olink_id","somascan_id","rho_olink_soma_non_ANML"}
    if not req<=set(h) or len(h)!=len(set(h)):
        raise ValueError("Unexpected reagent mapping")
    edges={}
    o=collections.defaultdict(set)
    s=collections.defaultdict(set)
    u=collections.defaultdict(set)
    for row in rows:
        r=dict(zip(h,row))
        if not all(r.get(k) for k in ("uniprot_id","olink_id","somascan_id")):
            continue
        up,a,b=(str(r[k]).strip() for k in ("uniprot_id","olink_id","somascan_id"))
        key=(up,a,b)
        rho=r["rho_olink_soma_non_ANML"]
        rho=float(rho) if rho is not None else None
        if rho is not None and (not math.isfinite(rho) or abs(rho)>1):
            raise ValueError("Invalid source rho")
        if key in edges and edges[key]!=rho:
            raise ValueError("Conflicting published mapping")
        edges[key]=rho
        o[a].add(key);s[b].add(key);u[up].add(key)
    w.close()
    pairs=[{"uniprot":up,"olink_id":a,"somascan_id":b,"published_rho_non_ANML":rho}
           for (up,a,b),rho in sorted(edges.items())
           if len(o[a])==len(s[b])==len(u[up])==1]
    if len(edges)<2000 or len(pairs)<1000:
        raise ValueError("Incomplete reagent mapping")
    return pairs,{"unique_edges":len(edges),"strict_one_to_one":len(pairs)}

def bh(pvalues):
    """Include untestable assays with p=1 in the denominator."""
    ordered=sorted(enumerate(pvalues),key=lambda item:(item[1],item[0]))
    out=[1.0]*len(pvalues)
    running=1.0
    for j,(i,p) in reversed(list(enumerate(ordered,1))):
        running=min(running,p*len(ordered)/j)
        out[i]=running
    return out

def summarize(rows,families,pairs):
    by={}
    all_hits=[]
    for platform in PLATFORMS:
        a=[r for r in rows if r["platform"]==platform]
        if len(a)!=families[platform] or len({r["assay_id"] for r in a})!=len(a):
            raise ValueError("Multiplicity universe changed")
        qs=bh([r.get("mr_p",1.0) for r in a])
        for r,q in zip(a,qs):
            r["bh_q_conservative_full_family"]=q
        observed=[r for r in a if r["status"]=="EXPLORATORY_WALD"]
        hits=[r for r in observed if r["bonferroni_p"]<.05]
        fhits=[r for r in observed if r["bh_q_conservative_full_family"]<.05]
        by[platform]={"family":len(a),"status_counts":dict(collections.Counter(r["status"] for r in a)),
                      "matched":len(observed),"bonferroni_hits":len(hits),"bh_hits":len(fhits)}
        all_hits.extend(hits)
    lookup={(r["platform"],r["assay_id"]):r for r in rows}
    comparisons=[]
    for pair in pairs:
        a=lookup.get(("Olink",pair["olink_id"]),{})
        b=lookup.get(("SomaScan",pair["somascan_id"]),{})
        aa=a.get("status")=="EXPLORATORY_WALD"
        bb=b.get("status")=="EXPLORATORY_WALD"
        same=aa and bb and a["rsid"]==b["rsid"] and (a["chr37"],a["pos37"])==(b["chr37"],b["pos37"])
        comparisons.append({**pair,"olink_status":a.get("status","NO_REPORTED_CIS"),
                            "soma_status":b.get("status","NO_REPORTED_CIS"),
                            "both_matched":aa and bb,"same_anchor":same,
                            "sign_agree":a["mr_beta"]*b["mr_beta"]>0 if aa and bb else None,
                            "both_bonferroni":aa and bb and a["bonferroni_p"]<.05 and b["bonferroni_p"]<.05,
                            "olink_bonferroni":a.get("bonferroni_p"),"soma_bonferroni":b.get("bonferroni_p"),
                            "olink_mr_beta":a.get("mr_beta"),"soma_mr_beta":b.get("mr_beta"),
                            "olink_rsid":a.get("rsid"),"soma_rsid":b.get("rsid")})
    both=[r for r in comparisons if r["both_matched"]]
    shared=[r for r in both if r["same_anchor"]]
    diff=[r for r in both if not r["same_anchor"]]
    joint=[r for r in both if r["both_bonferroni"]]
    ps={"strict_pairs":len(pairs),"both_matched":len(both),"same_anchor":len(shared),
        "same_anchor_direction_agree":sum(r["sign_agree"] for r in shared),
        "different_anchor":len(diff),"different_anchor_direction_agree":sum(r["sign_agree"] for r in diff),
        "both_bonferroni":len(joint),"both_bonferroni_direction_agree":sum(r["sign_agree"] for r in joint),
        "joint_hits":joint}
    keys=("platform","assay_id","rsid","mr_beta","mr_ci95","mr_p","bonferroni_p","bh_q_conservative_full_family","F")
    all_hits.sort(key=lambda r:(r["bonferroni_p"],r["platform"],r["assay_id"]))
    top=[{k:r.get(k) for k in keys} for r in all_hits[:100]]
    return by,ps,comparisons,top

def selftest():
    assert orient("A","G","G","A")==-1 and orient("A","T","T","A") is None
    x=wald(.5,.02,-.1,.01)
    assert abs(x["mr_beta"]+.2)<1e-12 and "beta" not in x
    assert bh([.01,.02,1.0,1.0])==[.04,.04,1.0,1.0]
    assert bh([1.0])==[1.0]
    assert bh([.01,.01,.5])==[.015,.015,.5]
    print("Synthetic harmonization, Wald and missing-inclusive BH tests passed.")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("trait",choices=["egfr","bun","selftest"])
    p.add_argument("--root",type=Path,default=Path("."))
    p.add_argument("--out",type=Path,default=Path("outputs"))
    a=p.parse_args()
    if a.trait=="selftest":
        selftest();return
    from pyliftover import LiftOver
    anchors,held,families,counts=assay_anchors(a.root/"ckb.xlsx")
    pairs,map_counts=mapping(a.root/"mapping.xlsx")
    anchor_hash=hashlib.sha256(json.dumps([anchors,held],sort_keys=True,separators=(",",":")).encode()).hexdigest()
    f=checked_chain(a.root,"hg38ToHg19.over.chain.gz","hg38.md5")
    b=checked_chain(a.root,"hg19ToHg38.over.chain.gz","hg19.md5")
    mapped=lift(anchors,LiftOver(str(f)),LiftOver(str(b)))
    meta=source_meta(a.root,a.trait)
    n,bad,h,results=scan(a.root/"outcome.gz",mapped,families)
    if n!=4512904 or bad:
        raise ValueError("Previously verified outcome row count mismatch")
    results.extend(held)
    by,ps,comparisons,top=summarize(results,families,pairs)
    summary=dict(status="EXPLORATORY_ASSAY_WIDE_SCREEN_NOT_CAUSAL",trait=a.trait,
                 source=meta,source_rows=n,source_headers=h,exposure_source_sha256=EXPOSURE_SOURCE,
                 mapping_sha256=MAPPING_SHA,anchor_selection_sha256=anchor_hash,
                 source_counts=counts,mapping_counts=map_counts,platforms=by,strict_pairs=ps,
                 top_bonferroni_hits=top,
                 limits=["All published cis-assay families in Supplementary Data 4, not the full CKB proteome.",
                         "Single strongest F anchor per assay selected without outcome; missing not replaced.",
                         "BUN supportive endpoint, EAS GWAS from TWB/BBJ previously inspected.",
                         "Same outcome and overlapping CKB exposures; paired signs not independent replication.",
                         "No PAV/epitope QC, dense regional coloc, independent EAS replication or ancestry difference test."])
    a.out.mkdir(parents=True,exist_ok=True)
    for name,value in (("summary.json",summary),("assay_results.json",results),("paired_results.json",comparisons)):
        (a.out/name).write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")
    print(json.dumps(summary,indent=2,allow_nan=False))

if __name__=="__main__":
    main()
