#!/usr/bin/env python3
"""P5f GEO donor/metadata attrition audit using original SOFT and published TSV.

Only aggregated counts; do not publish original donor IDs, ages, PRS, genotypes.
"""
import argparse,csv,gzip,json,re,tempfile
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from ocular_p0_audit import fetch,root,soft_samples
from ocular_p1_integrity import BASE,download
def norm(x):
    m=re.search(r"\bHTMC[- ]?(\d+)\b",x,re.I)
    return m.group(1) if m else None
def parse_soft(raw):
    rows=soft_samples(raw)
    across=defaultdict(set);by_condition=defaultdict(set)
    row=[]
    for x in rows:
        txt=" ".join(x["desc"])
        ids={m.group(1) for m in re.finditer(r"\bHTMC[- ]?(\d+)\b",txt,re.I)}
        title=x["title"].lower()
        arm="dex" if "dex" in title else "no-dex" if any(t in title for t in ("vehicle","control","pool-")) else "unknown"
        across[x["gsm"]]=ids;by_condition[arm]|=ids
        row.append({"gsm":x["gsm"],"treatment_guess":arm,"donors_listed":len(ids)})
    return rows,across,by_condition,row
def evaluate(soft_bytes,meta_file):
    samples,orig,by_arm,described=parse_soft(soft_bytes)
    original=set().union(*orig.values()) if orig else set()
    final=defaultdict(set);per_pool=defaultdict(set);counts=Counter()
    with gzip.open(meta_file,"rt",encoding="utf8") as f:
        rd=csv.DictReader(f,delimiter="\t")
        for r in rd:
            d=norm(r.get("donor_id_demuxlet",""))
            if d is None:
                raise ValueError("DEMUX_DONOR_ID_FORMAT_UNMAPPABLE")
            arm=r.get("DEX Treatment","").strip().lower()
            final[arm].add(d)
            pool=r.get("library_id","").split("pool-")[-1].split("_")[0].replace("dex","")
            if pool in "ABCDEFGHI":per_pool[pool].add(d)
            counts[arm]+=1
    both=final["dex"]&final["no-dex"]
    report={"status":"HOLD","dataset":"GSE313736","run_utc":datetime.now(timezone.utc).isoformat(),
      "source":"Official GEO SOFT sample descriptions and final cell_metadata.tsv.gz",
      "paper_original_donors_reported":74,
      "paper_original_singlet_cells_reported":316799,
      "original_libraries":len(samples),
      "original_sample_description_donors":len(original),
      "listed_donors_with_both_treatments":len(by_arm["dex"]&by_arm["no-dex"]),
      "final_metadata_unique_donors":len(set().union(*final.values())),
      "final_paired_donors":len(both),
      "final_metadata_cell_rows":sum(counts.values()),
      "paper_minus_final_singlets":316799-sum(counts.values()),
      "original_described_donors_missing_final":len(original-set().union(*final.values())),
      "final_donors_missing_original_descriptions":len(set().union(*final.values())-original),
      "original_donors_absent_from_SOFT":74-len(original),
      "reported_paper_donors_minus_final":74-len(set().union(*final.values())),
      "by_pool":[],
      "interpretation":["Enrollment and QC attrition cannot be established from these counts alone",
        "Paper's 316799 singlets may represent a different QC stage than final downloaded 169398 rows",
        "No donor identities or phenotypes included in audit evidence"]}
    for i,p in enumerate("ABCDEFGHI"):
        # Pool identity is defined by the exact official GSM accessions.
        # Sample titles are not guaranteed to contain "pool_A/B/etc".
        ctrl=f"GSM{9373276+i}"
        dex=f"GSM{9373285+i}"
        matched=orig.get(ctrl,set())|orig.get(dex,set())
        report["by_pool"].append({"pool":p,
          "SOFT_unique_donors":len(matched),
          "metadata_unique_donors":len(per_pool[p]),
          "metadata_subset_of_SOFT":len(per_pool[p]-matched)==0,
          "SOFT_minus_metadata_donors":len(matched-per_pool[p])})
    if (len(samples)==18 and len(original)>=68 and len(both)==59
        and sum(counts.values())==169398
        and all(r["SOFT_unique_donors"]>0 and r["metadata_subset_of_SOFT"]
                for r in report["by_pool"])):
        report["status"]="PASS_DESCRIPTIVE_ATTRITION_RECONCILIATION"
    return report
def main():
    a=argparse.ArgumentParser();a.add_argument("--out",required=True);args=a.parse_args()
    result={"status":"HOLD"}
    with tempfile.TemporaryDirectory() as td:
        try:
            raw,tx=fetch(root("GSE313736")+"soft/GSE313736_family.soft.gz",limit=9*1024*1024)
            dest=Path(td)/"meta.gz"
            result["metadata_transfer"]=download(
                BASE+"series/GSE313nnn/GSE313736/suppl/GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",
                dest,12*1024*1024)
            result=evaluate(raw,dest)
            result["SOFT_transfer"]=tx
        except Exception as exc:
            result["error"]=type(exc).__name__+": "+str(exc)[:200]
    Path(args.out).write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("P5F_ATTRITION",result["status"],
          "original_sample_donors",result.get("original_sample_description_donors"),
          "final",result.get("final_paired_donors"),
          "error",result.get("error"))
    if result["status"]!="PASS_DESCRIPTIVE_ATTRITION_RECONCILIATION":
        raise SystemExit(2)
if __name__=="__main__":main()
