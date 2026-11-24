#!/usr/bin/env python3
"""P6e: check ENA alternate raw-archive paths for P6c shotgun missing FASTQ metadata.

Metadata-only; never download patient samples or reads. A nonempty SRA/submitted
link does not imply analysis-ready paired FASTQ or validated tear-proteome join.
"""
import argparse,csv,io,json,urllib.parse,urllib.request,hashlib
from collections import Counter
from pathlib import Path

TARGET="PRJEB74248"
FIELDS="run_accession,library_strategy,fastq_ftp,submitted_ftp,sra_ftp,fastq_md5,submitted_md5,sra_md5"
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    query=urllib.parse.urlencode({"accession":TARGET,"result":"read_run","fields":FIELDS,"format":"tsv"})
    req=urllib.request.Request("https://www.ebi.ac.uk/ena/portal/api/filereport?"+query,
                              headers={"User-Agent":"ocular-p6e-metadata-path-audit/1.0"})
    with urllib.request.urlopen(req,timeout=90) as h:raw=h.read(2000001)
    if len(raw)>2000000:raise ValueError("ENA_METADATA_TOO_LARGE")
    reader=csv.DictReader(io.StringIO(raw.decode("utf8","replace")),delimiter="\t")
    need={"run_accession","fastq_ftp","submitted_ftp","sra_ftp","library_strategy"}
    if not need.issubset(reader.fieldnames or []):
        raise ValueError("ENA_ALTERNATIVE_FIELD_SCHEMA_INVALID")
    data=list(reader)
    if len(data)!=62 or len({r["run_accession"] for r in data})!=62:
        raise ValueError("PREVIOUS_62_RUN_SOURCE_DRIFT")
    counts=Counter()
    combos=Counter()
    strategies=Counter()
    for r in data:
        flags={k:bool((r.get(k) or "").strip()) for k in ("fastq_ftp","submitted_ftp","sra_ftp")}
        for name,seen in flags.items():
            counts[name+"_present"]+=int(seen)
        counts["at_least_one_raw_route_present"]+=int(any(flags.values()))
        counts["fastq_absent_with_other_route"]+=int(not flags["fastq_ftp"] and (flags["submitted_ftp"] or flags["sra_ftp"]))
        counts["all_routes_absent"]+=int(not any(flags.values()))
        combos[",".join(k for k,v in flags.items() if v) or "NONE"]+=1
        strategies[(r.get("library_strategy") or "UNKNOWN").upper()]+=1
    if counts["fastq_ftp_present"]!=31:
        raise ValueError("PREVIOUS_FASTQ_31_OF_62_DRIFT")
    status="PASS_METADATA_ONLY_ALTERNATIVE_ROUTE_RECONCILED"
    outcome={"status":status,"accession":TARGET,
        "total_runs":len(data),"library_strategy_counts":dict(strategies),
        "run_route_counts":dict(counts),"route_combinations":dict(combos),
        "source_metadata_sha256":hashlib.sha256(raw).hexdigest(),
        "sample_fastq_or_sra_read_downloaded":False,
        "raw_reads_reprocessed":False,
        "clinical_or_tear_proteome_linked":False,
        "missing_fastq_is_not_same_as_missing_archive":True,
        "no_link_proves_downloadability":True,
        "note":"Availability of SRA or submitted FTP metadata does not guarantee paired FASTQ, complete host decontamination, or matched proteomics. Validate URLs, read layout, case/control covariates and low-biomass negative controls before project GO."}
    (args.out/"P6e_shotgun_alternative_access_gate.json").write_text(
        json.dumps(outcome,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("P6E",json.dumps({"status":status,"route_counts":dict(counts),"combinations":dict(combos)},ensure_ascii=False),flush=True)
if __name__=="__main__":
    main()
