#!/usr/bin/env python3
"""Mehrotra 2026 erratum body and newly disclosed PRJNA1243143 ENA manifest.

No FASTQ download. Only original official fullTextXML (PMC13712379)
and ENA project metadata. Zero assumption that N technical runs = N donors.
"""
import argparse,csv,hashlib,io,json,re,time,urllib.parse,urllib.request
import xml.etree.ElementTree as ET
from collections import Counter,defaultdict
from pathlib import Path

PMC="PMC13712379";PROJECT="PRJNA1243143"
ERR="https://www.ebi.ac.uk/europepmc/webservices/rest/"+PMC+"/fullTextXML"
ENA="https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS="study_accession,run_accession,experiment_accession,sample_accession,sample_alias,sample_title,library_layout,library_strategy,fastq_ftp,fastq_md5,fastq_bytes,submitted_ftp,submitted_md5,submitted_bytes"
SMALL_FIELDS="study_accession,run_accession,sample_accession,sample_alias,library_layout,fastq_ftp,fastq_md5,fastq_bytes"

def get(url,cap=9_000_000):
    last=""
    for attempt in range(3):
        try:
            req=urllib.request.Request(url,headers={
              "User-Agent":"ocular-erratum-and-ENA-public-sample-audit/1.0",
              "Accept":"application/xml, text/tab-separated-values, */*"})
            with urllib.request.urlopen(req,timeout=48) as f:
                raw=f.read(cap+1);final=f.url
            if len(raw)>cap or not raw:raise ValueError("RESPONSE_EMPTY_OR_TOO_BIG")
            return raw,final
        except (OSError,ValueError) as e:
            last=type(e).__name__+" "+str(e)[:150]
            if attempt<2:time.sleep(2)
    raise RuntimeError(last)

def xml_node_text(node):
    return re.sub(r"\s+"," "," ".join(node.itertext())).strip()

def classify_erratum(raw):
    root=ET.fromstring(raw)
    bodies=root.findall(".//body")
    if len(bodies)!=1:raise ValueError("ERRATUM_BODY_NOT_UNIQUE")
    paras=[xml_node_text(node) for node in bodies[0].iter() if node.tag.endswith("p")]
    text=" ".join(paras)
    confirms=(PROJECT in text and "Data Availability" in text and
              "Acknowledgments" in text and "raw RNA-seq" in text)
    if not confirms:raise ValueError("ERRATUM_NO_VERIFIABLE_DEPOSITION_STATEMENT")
    # The notice body is compact. It explicitly adds data-availability text;
    # do not claim original supplement file hashes were themselves rechecked.
    is_only=("This Data Availability section has been added" in text and
             PROJECT in text and not any(
       k in text.lower() for k in ["incorrect data","table 2", "supplementary table",
                              "revised values","figure was","p-values","fold changes"]))
    return {"status":"PASS_ERRATUM_DEPOSITION_STATEMENT_VERIFIED",
      "erratum_pmc":PMC,"doi":"10.1167/iovs.67.12.34",
      "xml_sha256":hashlib.sha256(raw).hexdigest(),
      "body_paragraph_count":len(paras),
      "notice_only_data_availability_addition":is_only,
      "notice_evidence_excerpt":text[:950],
      "newly_disclosed_raw_RNAseq_project":PROJECT,
      "explicit_numeric_deq_correction_in_this_notice":False if is_only else "UNKNOWN",
      "prior_supplement_file_bytes_reverified_in_this_step":False,
      "analysis_conclusion":"No DEG/effect-value correction stated in the erratum notice; primary source table unchanged cannot be independently asserted from notice alone."}

def parse_table(raw):
    txt=raw.decode("utf8","replace")
    rdr=csv.DictReader(io.StringIO(txt),delimiter="\t")
    if not rdr.fieldnames or "run_accession" not in rdr.fieldnames:
        raise ValueError("ENA_RUN_METADATA_HEADER_MISSING "+txt[:150])
    data=list(rdr)
    if len(data)>20000:raise ValueError("ENA_EXCESSIVE_RUN_COUNT")
    return data

def file_list(row,name):
    return [x for x in row.get(name,"").split(";") if x]
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path)
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    state={"status":"HOLD","run_accession":PROJECT,"source_kind":"official ENA read_run metadata"}
    try:
        x,final=get(ERR,2_000_000)
        correction=classify_erratum(x)
        (args.out/"MEHROTRA_ERRATUM_ADDENDUM_VERIFIED.json").write_text(
            json.dumps(correction,indent=2,ensure_ascii=False)+"\n")
    except Exception as e:
        correction={"status":"HOLD_OFFICIAL_ERRATUM_XML",
             "error":type(e).__name__+" "+str(e)[:200]}
        (args.out/"MEHROTRA_ERRATUM_ADDENDUM_VERIFIED.json").write_text(
            json.dumps(correction,indent=2,ensure_ascii=False)+"\n")
    attempts=[]
    requested=[FIELDS,SMALL_FIELDS]
    data=None;response=None;used=None
    for fields in requested:
        query=urllib.parse.urlencode({"accession":PROJECT,"result":"read_run",
                      "fields":fields,"format":"tsv"})
        url=ENA+"?"+query
        try:
            raw,final=get(url,12_000_000)
            data=parse_table(raw)
            response={"sha256":hashlib.sha256(raw).hexdigest(),
                      "response_bytes":len(raw),"final_url":final}
            used=fields
            break
        except Exception as e:
            attempts.append({"field_selection":fields,"error":type(e).__name__+" "+str(e)[:200]})
    if data is None:
        state.update({"status":"HOLD_ENA_PROJECT_RUN_REPORT_NOT_ACCESSIBLE",
                      "attempts":attempts})
    else:
        runs=[r for r in data if re.fullmatch(r"[EDS]RR[0-9]+",r.get("run_accession",""))]
        unq=len({r["run_accession"] for r in runs})==len(runs)
        files=0;totalbytes=0;md5covered=0;empty=0;invalidlength=0
        sanitized=[]
        for r in runs:
            fastqs=file_list(r,"fastq_ftp")
            sums=file_list(r,"fastq_md5")
            sizes=file_list(r,"fastq_bytes")
            submitted=file_list(r,"submitted_ftp")
            submitted_sums=file_list(r,"submitted_md5")
            submitted_sizes=file_list(r,"submitted_bytes")
            flen=len(fastqs)
            if not flen:empty+=1
            files+=flen
            if flen and len(sums)==flen and all(re.fullmatch("[0-9a-fA-F]{32}",m) for m in sums):
                md5covered+=flen
            if len(sizes)==flen:
                try:totalbytes+=sum(int(x) for x in sizes)
                except ValueError:invalidlength+=1
            elif flen:invalidlength+=1
            sanitized.append({
              "run_accession":r["run_accession"],
              "sample_accession":r.get("sample_accession",""),
              "sample_alias":r.get("sample_alias","")[:120],
              "sample_title":r.get("sample_title","")[:180],
              "library_layout":r.get("library_layout",""),
              "library_strategy":r.get("library_strategy",""),
              "fastq_file_count":flen,
              "fastq_md5_complete":bool(flen and len(sums)==flen),
              "fastq_bytes":sum(int(z) for z in sizes) if sizes and all(z.isdigit() for z in sizes) else "",
              "submitted_file_count":len(submitted),
              "submitted_md5_complete":len(submitted)==len(submitted_sums) if submitted else "",
              "submitted_bytes":sum(int(z) for z in submitted_sizes) if submitted_sizes and all(z.isdigit() for z in submitted_sizes) else "",
              "file_locations_not_downloaded":True
            })
        state.update({"source_official_URL":response["final_url"],
            "response_sha256":response["sha256"],
            "ENA_query_fields":used,
            "run_count":len(runs),"unique_run_accessions":unq,
            "unique_sample_accessions":len({r.get("sample_accession") for r in runs}),
            "unique_sample_aliases":len({r.get("sample_alias") for r in runs if r.get("sample_alias")}),
            "fastq_files":files,"fastq_ftp_missing_run_count":empty,
            "fastq_file_md5_covered":md5covered,
            "fastq_sum_bytes":totalbytes if invalidlength==0 else None,
            "bytes_information_complete":invalidlength==0,
            "strategy_counts":dict(Counter(r.get("library_strategy","") for r in runs)),
            "layout_counts":dict(Counter(r.get("library_layout","") for r in runs)),
            "metadata_only_no_sequence_download":True,
            "donor_unit_from_runs_NOT_ASSUMED":True,
            "dose_time_group_mapping":"HOLD_UNTIL_EXACT_SOURCE_DESIGN_AND_RUN_SAMPLE_CROSSWALK",
            "attempts":attempts})
        # Count matching sequencing-library pairs only when a shared sample
        # stem is unambiguous. "-1/-2" denotes replicate/run identifiers, NOT
        # proven Vehicle vs DEX treatment assignment or timepoint.
        pairs=defaultdict(list)
        for rec in sanitized:
            alias=rec["sample_alias"]
            stem=re.sub(r"(?i)_(Con-EtOH|Dex)_Trypsin(?:_\d+)?$","",alias)
            stem=re.sub(r"(?i)_(vehicle|dex)$","",stem)
            stem=re.sub(r"-[12]$","",stem)
            if re.search(r"(?i)(con-etoh|vehicle)",alias):
                arm="Vehicle_explicit"
            elif re.search(r"(?i)(Dex_Trypsin|_dex$)",alias):
                arm="DEX_explicit"
            else:
                arm="UNKNOWN_-1_-2"
            pairs[stem].append((rec,arm))
        pair_rows=[]
        for key,records in sorted(pairs.items()):
            arms=[arm for _,arm in records]
            confident=len(records)==2 and sorted(arms)==["DEX_explicit","Vehicle_explicit"]
            category="SCE" if key.upper().startswith("SC") else "TM" if key.upper().startswith(("TM","HTM")) else "UNKNOWN"
            pair_rows.append({
              "source_pair_stem":key,
              "tissue_from_alias":category,
              "run_count":len(records),
              "treatment_pair_status":"EXPLICIT_VEHICLE_DEX" if confident else
                  "HOLD_TREATMENT_UNKNOWN_FOR_-1_-2",
              "run_accessions":";".join(x["run_accession"] for x,_ in records),
              "sample_aliases":";".join(x["sample_alias"] for x,_ in records),
              "fastq_bytes":sum(int(x["fastq_bytes"]) for x,_ in records if isinstance(x["fastq_bytes"],int)),
              "treatment_timepoint_not_verified":True,
              "repeat_samples_may_share_donor_strain":True})
        confident_count=sum(x["treatment_pair_status"]=="EXPLICIT_VEHICLE_DEX" for x in pair_rows)
        uncertain_count=sum(x["treatment_pair_status"]!="EXPLICIT_VEHICLE_DEX" for x in pair_rows)
        state["source_alias_pair_count"]=len(pair_rows)
        state["explicit_vehicle_dex_pair_count"]=confident_count
        state["unassigned_dash1_dash2_pairs"]=uncertain_count
        state["tissue_pair_counts"]=dict(Counter(x["tissue_from_alias"] for x in pair_rows))
        state["sampling_timepoints_verified_from_ENA_alias"]=0
        state["pairing_note"]="Pair alias stem is sequencing-library grouping; cannot assume independent human donors, or deduce 1h/6h/48h without source mapping."
        state["source_exact_23_pairs_expected"]=len(pair_rows)==23 and confident_count==17 and uncertain_count==6 and all(x["run_count"]==2 for x in pair_rows)
        if pair_rows:
            with (args.out/"MEHROTRA_PRJNA1243143_PAIRING_QC.tsv").open("w",newline="") as h:
                writer=csv.DictWriter(h,fieldnames=list(pair_rows[0]),delimiter="\t")
                writer.writeheader();writer.writerows(pair_rows)
        state["status"]="GO_OFFICIAL_ENA_RUN_METADATA_ACCESS_ONLY" if runs and unq else "HOLD_ENA_RUNS_UNUSABLE"
        if sanitized:
            with (args.out/"MEHROTRA_PRJNA1243143_ENA_RUN_MANIFEST.tsv").open("w",newline="") as h:
                writer=csv.DictWriter(h,fieldnames=list(sanitized[0]),delimiter="\t")
                writer.writeheader();writer.writerows(sorted(sanitized,key=lambda r:r["run_accession"]))
    state["erratum_status"]=correction["status"]
    (args.out/"MEHROTRA_PRJNA1243143_METADATA_READINESS.json").write_text(
        json.dumps(state,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("ERRATUM",correction["status"],"deposition",correction.get("newly_disclosed_raw_RNAseq_project"))
    print("ENA",state["status"],"runs",state.get("run_count"),
          "samples",state.get("unique_sample_accessions"),
          "fastqs",state.get("fastq_files"),"bytes",state.get("fastq_sum_bytes"))
if __name__=="__main__":main()
