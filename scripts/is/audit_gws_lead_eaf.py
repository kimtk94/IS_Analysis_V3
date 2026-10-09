#!/usr/bin/env python3
"""Compare lead-variant ALT EAF in EAS GWAS to 1KG EAS genotype reference.

Read-only original GWAS. Frequency discrepancies are QC flags, not genetic effect evidence.
"""
import argparse,csv,gzip,json,math,subprocess
from pathlib import Path
from collections import defaultdict,Counter
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
NEW=Path("/srv/is-analysis/data/is/ld_reference/broad_v1/1kg_eas")
OLD=Path("/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/pgen_gt")
LOCUS={"IS_XDATA_G0007":"L001","IS_XDATA_G0015":"L002",
       "IS_XDATA_G0022":"L003","IS_XDATA_G0023":"L004"}
def source_prefix(group,new,old):
    return old/f"BBJ_IS_{LOCUS[group]}.1KG_EAS.GRCh37" if group in LOCUS else new/(group+".stable")
def floats(s):
    try:
        v=float(s)
        return v if math.isfinite(v) else None
    except (TypeError,ValueError):return None
def run(root,new,old):
    with (root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv").open() as f:
        groups={r["group_id"]:r for r in csv.DictReader(f,delimiter="\t") if r["has_gws"]=="1"}
    if len(groups)!=7:raise ValueError("Seven GWS regions required")
    out=root/"reference_frequency_qc";out.mkdir(exist_ok=True)
    allele_freq={}
    for group in sorted(groups):
        prefix=source_prefix(group,new,old)
        afreq=out/f"{group}.afreq"
        if not afreq.is_file():
            p=subprocess.run(["plink2","--pfile",str(prefix),"--freq",
                "--out",str(out/group)],capture_output=True,text=True,timeout=180)
            if p.returncode:raise RuntimeError("plink frequency failed for "+group+" "+p.stderr+p.stdout[-400:])
        lead=groups[group]["lead_variant"].split(":")
        with afreq.open() as f:
            for r in csv.DictReader(f,delimiter="\t"):
                pospair=r["ID"].split(":")
                if len(pospair)>=4 and pospair[:2]==lead[:2] and set(pospair[2:4])==set(lead[2:4]):
                    freq=floats(r["ALT_FREQS"])
                    allele_freq[group]=(r["ID"],r["ALT"],freq,int(r["OBS_CT"]))
                    break
        if group not in allele_freq:raise ValueError("Lead allele pair missing from reference frequency: "+group)
    with (root/"IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz").open("rb") as f:
        if f.read(2)!=b"\x1f\x8b":raise ValueError("Incomplete matching reference file")
    data=defaultdict(list)
    with gzip.open(root/"IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz","rt") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            g=r["group_id"]
            if g not in groups or int(r["pos"])!=int(groups[g]["lead_variant"].split(":")[1]):
                continue
            data[g,r["dataset"]].append(r)
    table=[]
    for (g,ds),records in sorted(data.items()):
        refid,alt,reffreq,n=allele_freq[g]
        for record in records:
            if record["variant_id"]!=refid:continue
            gwaseaf=floats(record["alt_effect_eaf"])
            delta=abs(gwaseaf-reffreq) if gwaseaf is not None and reffreq is not None else None
            table.append({"group_id":g,"dataset":ds,"lead_variant":groups[g]["lead_variant"],
                "reference_variant_id":refid,"ref_alt":alt,"ref_alt_freq":reffreq,
                "reference_allele_observation_count":n,"gwas_alt_eaf":gwaseaf if gwaseaf is not None else "",
                "abs_eaf_delta":round(delta,6) if delta is not None else "",
                "qc":"EAF_DIFFERENCE_REVIEW" if delta is not None and delta>0.10 else
                     ("FREQUENCY_UNAVAILABLE" if delta is None else "EAF_DIFFERENCE_LE_0_10"),
                "allele_status":record["qc_status"],"ancestry":record["ancestry"],
                "caveat":"Allele frequency similarity is a QC check, not proof of ancestry matching"})
    fields=["group_id","dataset","lead_variant","reference_variant_id","ref_alt",
      "ref_alt_freq","reference_allele_observation_count","gwas_alt_eaf","abs_eaf_delta",
      "qc","allele_status","ancestry","caveat"]
    with (root/"IS_GWS_LEAD_EAF_QC.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
        w.writeheader();w.writerows(table)
    summary={"lead_groups_reference_freq":len(allele_freq),
      "lead_gwas_frequency_pairs":len(table),
      "difference_above_0_10":sum(x["qc"]=="EAF_DIFFERENCE_REVIEW" for x in table),
      "difference_at_most_0_10":sum(x["qc"]=="EAF_DIFFERENCE_LE_0_10" for x in table),
      "missing_frequency":sum(x["qc"]=="FREQUENCY_UNAVAILABLE" for x in table),
      "interpretation":"Study EAF vs 1KG EAS reference ALT AF, descriptive check only"}
    (root/"IS_GWS_LEAD_EAF_QC_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    return summary
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--ref-new",type=Path,default=NEW)
    p.add_argument("--ref-old",type=Path,default=OLD)
    a=p.parse_args()
    run(a.root,a.ref_new,a.ref_old)
