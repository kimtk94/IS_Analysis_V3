#!/usr/bin/env python3
"""Verify five GIGASTROKE East Asian stroke outcomes with live Catalog v2 metadata.

Records official sample sizes but deliberately does not approve binary N_eff,
cohort independence, variant coverage, or a causal interpretation.
"""
from __future__ import annotations
import argparse,json,re,time,urllib.error,urllib.request
from pathlib import Path
BASE="https://www.ebi.ac.uk/gwas/rest/api/v2/studies/"
FTP=("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
     "GCST90104001-GCST90105000/")
TARGETS={
 "GCST90104544":("AS","MONDO_0005098"),
 "GCST90104545":("AIS","HP_0002140"),
 "GCST90104546":("CES","EFO_1001976"),
 "GCST90104547":("LAS","MONDO_0005490"),
 "GCST90104548":("SVS","EFO_1001504"),
}
def get_json(url,tries=4):
    for attempt in range(tries):
        req=urllib.request.Request(url,headers={"Accept":"application/json",
                            "User-Agent":"TPMI-43321-EAS-outcome-metadata-review/1.0"})
        try:
            with urllib.request.urlopen(req,timeout=40) as res:
                if res.status!=200 or res.geturl()!=url:raise ValueError("API redirected")
                blob=res.read(180000)
            if len(blob)>=180000:raise ValueError("oversized API response")
            return json.loads(blob.decode("utf-8-sig"))
        except urllib.error.HTTPError as exc:
            if exc.code not in (429,502,503,504) or attempt==tries-1:
                raise
            time.sleep(3*(attempt+1))
def parse_cases_controls(s):
    m=re.fullmatch(
      r"([\d,]+) East Asian ancestry cases, ([\d,]+) East Asian ancestry controls",
      s.strip()) if isinstance(s,str) else None
    return (int(m[1].replace(",","")),int(m[2].replace(",",""))) if m else (None,None)
def audit(acc,expect):
    row={"accession":acc,"expected_phenotype_code":expect[0],
         "expected_efo_id":expect[1],"metadata_status":"HOLD"}
    try:
        d=get_json(BASE+acc)
        if d.get("accession_id")!=acc:
            raise ValueError("Study accession mismatch")
        traits=d.get("efo_traits",[])
        identifiers={x.get("efo_id") for x in traits if isinstance(x,dict)}
        if expect[1] not in identifiers:
            raise ValueError("Trait does not match official ontology")
        label=d.get("initial_sample_size")
        cases,controls=parse_cases_controls(label)
        row.update({"reported_initial_sample_size":label,
          "reported_replication_sample_size":d.get("replication_sample_size"),
          "efo_traits":traits,"discovery_cases":cases,"discovery_controls":controls})
        if cases is None or controls is None:
            raise ValueError("No uniquely described EAS discovery cases/controls")
        row["metadata_status"]="PASS_EAS_STUDY_TRAIT_AND_COUNTS"
    except Exception as e:
        row["error"]=f"{type(e).__name__}: {e}"[:220]
    row["file_url"]=FTP+f"{acc}/{acc}_buildGRCh37.tsv.gz"
    row["source_variant_coverage"]="NOT_CHECKED"
    row["sample_independence"]="NOT_CHECKED"
    row["sample_effective_n_approved"]=False
    return row
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    rows=[]
    for acc,expected in TARGETS.items():
        rows.append(audit(acc,expected))
        time.sleep(1.3)
    out={"schema":"GIGASTROKE-EAS-5-outcome-audit-v1",
         "source":"GWAS Catalog REST v2 (2026) and EBI official FTP GRCh37",
         "EAS_phenotype_count":len(rows),
         "EAS_metadata_pass":sum(x["metadata_status"].startswith("PASS") for x in rows),
         "studies":rows,
         "scientific_status":"HOLD",
         "cohort_overlap_not_resolved":True,
         "genomewide_variant_overlap_not_resolved":True,
         "N_eff_not_approved":True,
         "note":"Shared controls make subtype outcomes statistically dependent; cross-ancestry results reuse EAS samples."}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"EAS_metadata_pass":out["EAS_metadata_pass"],
                      "studies":rows,"scientific_status":"HOLD"},ensure_ascii=False))
if __name__=="__main__":main()
