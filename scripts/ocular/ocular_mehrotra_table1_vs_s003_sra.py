#!/usr/bin/env python3
"""Cross-audit official IOVS Table 1, official s003.xlsx marker columns, and SRA run-pair candidates.

The *paper* lists TM155 as 6h; supplement header has TM151 6h.
Preserve disagreements; never silently assign or correct experiment labels.
Raw gene counts, FASTQ, private donor metadata are not used here.
"""
import argparse,csv,hashlib,json,re,urllib.request,xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
PMC_XML="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12882097/fullTextXML"
PMC_SHA="87a80692f387318037da1d6b7aa45b7a7e546e888abac92a37796e1b274f64c2"
S003_SHA="3c30e93461812d05f42a1978b8e446235709f2283e11f620c212beb26478698d"
TIMEPOINTS=("2Days","6Hrs","1Hr")
def text(e):
    return re.sub(r"\s+"," "," ".join(e.itertext())).strip()
def official_table1():
    req=urllib.request.Request(PMC_XML,headers={"User-Agent":"academic-published-source-table-audit/1.0"})
    with urllib.request.urlopen(req,timeout=45) as response:
        raw=response.read(9_000_001)
    if len(raw)>9_000_000:raise ValueError("PMC_XML_OVER_CAP")
    digest=hashlib.sha256(raw).hexdigest()
    if digest!=PMC_SHA:raise ValueError("PMC_PRIMARY_XML_SHA256_CHANGED")
    root=ET.fromstring(raw)
    matching=[]
    for tw in root.findall(".//table-wrap"):
        heading=" ".join(text(x) for x in tw.findall("./caption"))
        if "Demographic Data of Donors" in heading:matching.append(tw)
    if len(matching)!=1:raise ValueError("ORIGINAL_TABLE1_NOT_UNIQUE")
    table=matching[0].find(".//table")
    if table is None:raise ValueError("TABLE_NO_HTML_TABLE")
    records={}
    for tr in table.findall(".//tr"):
        cells=[text(x) for x in tr if x.tag in ("th","td")]
        if not cells:continue
        m=re.search(r"\b(TM|SC)\s*(\d{2,3})\b",cells[0],re.I)
        if not m:continue
        tag=m.group(1).upper()+m.group(2)
        if tag in records:raise ValueError("DUPLICATE_TABLE1_STRAIN "+tag)
        col=cells[-1]
        which=set()
        for match in re.finditer(r"\b(2\s*d|1\s*h|6\s*h|48\s*h|2\s*days)\b",col,re.I):
            item=match.group(1).lower().replace(" ","")
            which.add("2Days" if item in ("2d","2days","48h") else "1Hr" if item=="1h" else "6Hrs")
        if not which:raise ValueError("PAPER_STRAIN_WITHOUT_TIMES "+tag+repr(cells))
        records[tag]={"strain":tag,"timepoints":sorted(which),"raw_timepoint_text":col}
    if len([x for x in records if x.startswith("TM")])!=9 or len([x for x in records if x.startswith("SC")])!=5:
        raise ValueError("TABLE1_EXPECTED_9_TM_5_SC")
    return records,digest

def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf8"))
def load_tsv(path):
    with Path(path).open(encoding="utf8",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def save_tsv(path,rows):
    if not rows:raise ValueError("ZERO_ROWS")
    with Path(path).open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader();w.writerows(rows)
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-gate",type=Path,required=True)
    p.add_argument("--pairs",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    source=load_json(a.source_gate)
    if source.get("source_s003_sha256")!=S003_SHA or source.get("candidate_48h_TM_library_pairs")!=9:
        raise ValueError("PREVIOUS_SOURCE_GATE_SHA_OR_N_DRIFT")
    design=source["supplement_exact_strains_by_time"]
    table,sha=official_table1()
    comparison=[];setdiff=[]
    for tissue in ("TM","SC"):
        for time in TIMEPOINTS:
            paper={strain for strain,d in table.items() if strain.startswith(tissue) and time in d["timepoints"]}
            supplement=set(design[tissue][time])
            setdiff.append({
              "tissue":tissue,"timepoint":time,
              "primary_Table1_n":len(paper),"publisher_s003_marker_n":len(supplement),
              "primary_only_strains":";".join(sorted(paper-supplement)),
              "supplement_only_strains":";".join(sorted(supplement-paper)),
              "exact_strain_set_match":paper==supplement,
              "source_inconsistency_requires_author_clarification":paper!=supplement})
    pairs=load_tsv(a.pairs)
    if len(pairs)!=23:raise ValueError("EXPECTED_23_NOMINAL_PAIRS")
    for row in pairs:
        tag=row["nominal_strain"]
        time=row["candidate_exposure"]
        tissue=row["tissue_from_source_alias"]
        in_table=time in table[tag]["timepoints"] if tag in table else False
        in_supp=tag in design[tissue].get(time,[])
        comparison.append({
          "source_library_pair_stem":row["source_library_pair_stem"],
          "tissue":tissue,"nominal_strain":tag,
          "candidate_exposure_from_alias":time,
          "paper_Table1_lists_strain_at_time":in_table,
          "supplement_s003_lists_strain_at_time":in_supp,
          "paper_and_s003_agree_for_pair":in_table and in_supp,
          "direct_RUN_timepoint_attribute_proven":False,
          "run_accessions":row["run_accessions"],
          "decision":"REVIEW_THREE_SOURCE_CANDIDATE_NOT_RUN_TIMEPROVEN" if in_table and in_supp else
                     "HOLD_PAPER_TABLE1_VS_SUPPLEMENT_DISCORDANT",
          "never_assign_author_confirmed_exposure_without_RUN_sheet":True})
    unexpected=[r for r in comparison if not r["paper_and_s003_agree_for_pair"]]
    conflicts=[r for r in setdiff if not r["exact_strain_set_match"]]
    good48=[r for r in comparison if r["tissue"]=="TM" and
            r["candidate_exposure_from_alias"]=="2Days" and r["paper_and_s003_agree_for_pair"]]
    if len(good48)!=9 or len(unexpected)!=1 or len(conflicts)!=1:
        raise ValueError("UNEXPECTED_TABLE1_VS_SUPPLEMENT_CONCORDANCE "+str({
           "good48":len(good48),"unexpected":unexpected,"conflicts":conflicts}))
    if unexpected[0]["nominal_strain"]!="TM151" or unexpected[0]["candidate_exposure_from_alias"]!="6Hrs":
        raise ValueError("DISCORDANT_STRAIN_NOT_EXPECTED_TM151_6HRS")
    if conflicts[0]["primary_only_strains"]!="TM155" or conflicts[0]["supplement_only_strains"]!="TM151":
        raise ValueError("CONFLICT_CONTENT_NOT_TM155_VS_TM151")
    save_tsv(a.out/"IOVS_Table1_vs_S003_by_timepoint.tsv",setdiff)
    save_tsv(a.out/"IOVS_Table1_S003_SRA_23pair_triple_source_audit.tsv",comparison)
    summary={
      "status":"PASS_THREE_SOURCE_QC_IDENTIFIED_ONE_TM6HRS_CONFLICT",
      "source_primary_PMC12882097_xml_sha256":sha,
      "source_s003_xlsx_sha256":S003_SHA,
      "source_Table1_named_TM_strains":9,
      "source_Table1_named_SC_strains":5,
      "timepoint_stratum_comparisons":len(setdiff),
      "timepoint_strata_discordant_n":len(conflicts),
      "discordant_stratum":conflicts[0],
      "SRA_candidate_pair_n":len(comparison),
      "nominal_three_source_candidate_matches_n":sum(r["paper_and_s003_agree_for_pair"] for r in comparison),
      "candidate_pair_discordance_n":len(unexpected),
      "candidate_pair_discordance":unexpected,
      "48h_TM_candidate_pairs_supported_by_both_official_published_sources":len(good48),
      "48h_TM_candidate_RUN_count":18,
      "RUN_specific_exposure_time_confirmed_n":0,
      "48h_TM_biological_donor_and_QC_confirmation":"HOLD",
      "conclusion":"Source Table1 lists TM155 at 6h whereas original Supplement 3 Table82 lists TM151 at 6h; this does not affect the nine 2-day TM candidate strain names but direct SRA exposure-time mapping is still missing.",
      "no_gene_level_or_FASTQ_data_used":True}
    (a.out/"IOVS_Table1_versus_S003_SRA_threeway_gate.json").write_text(
       json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("THREE_WAY_SOURCE_AUDIT",json.dumps(summary,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
