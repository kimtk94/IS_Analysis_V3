#!/usr/bin/env python3
"""Constrained candidate matching: official IOVS 2d/6h/1h marker sheet vs ENA.

This algorithm NEVER labels a RUN as definitively 48h without sample-specific
time in its original SRA record. Candidate aliases are inferred by source
naming convention + source paper sample/strain/time set agreement.
"""
import argparse,csv,json,re
from collections import defaultdict,Counter
from pathlib import Path

TM_SHEET="SuppTable82_TM_MarkerExpression"
SC_SHEET="SuppTable83_SC_MarkerExpression"
TIMES=("2Days","6Hrs","1Hr")
def ftsv(path):
    with path.open(encoding="utf8",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))
def write(path,rows):
    if not rows:raise ValueError("NO_OUTPUT_ROWS")
    with path.open("w",encoding="utf8",newline="") as f:
        d=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        d.writeheader();d.writerows(rows)
def source_design(filename):
    audit=json.loads(filename.read_text(encoding="utf8"))
    source=next((x for x in audit if x["filename"]=="iovs-67-2-13_s003.xlsx"),None)
    if not source or source["source_status"]!="PASS_OFFICIAL_SUPPLEMENT_DOWNLOADED_AND_CRC":
        raise ValueError("S003_NOT_GENUINE_OFFICIAL_SOURCE_CRC")
    sheets=source["audit"]["sheet_first_300_rows_summary"]
    sample_tables={}
    for name in (TM_SHEET,SC_SHEET):
        row=next((x for x in sheets if x["workbook_sheet_name_by_ordinal"]==name),None)
        if not row or not row["header_rows"]:raise ValueError("SOURCE_MARKER_TABLE_HEADER_NOT_FOUND "+name)
        sample_tables[name]=row["header_rows"][0]
    originals={}
    for tissue,table in (("TM",TM_SHEET),("SC",SC_SHEET)):
        d=defaultdict(set)
        for column in sample_tables[table]:
            m=re.fullmatch(r"(?i)"+tissue+r"(\d+)_CTRL_(2Days|6Hrs|1Hr)",column)
            if m:d[m.group(2)].add(tissue+m.group(1))
        originals[tissue]={t:sorted(d[t]) for t in TIMES}
    # Exact positive-control counts in original article Table 82/83.
    if {k:len(v) for k,v in originals["TM"].items()}!={"2Days":9,"6Hrs":4,"1Hr":3}:
        raise ValueError("TM_ORIGINAL_EXPOSURE_GROUP_DESIGN_DRIFT")
    if {k:len(v) for k,v in originals["SC"].items()}!={"2Days":4,"6Hrs":4,"1Hr":3}:
        raise ValueError("SC_ORIGINAL_EXPOSURE_GROUP_DESIGN_DRIFT")
    return originals,source["sha256"],source["compressed_bytes"]

def group_alias(name):
    # Only the *pair stems*, not arbitrary row order, may be grouped.
    stem=re.sub(r"(?i)_(?:Con-EtOH|Dex)_Trypsin(?:_\d+)?$","",name)
    stem=re.sub(r"(?i)_(?:vehicle|dex)$","",stem)
    stem=re.sub(r"-[12]$","",stem)
    return stem
def pattern_candidate(stem,aliases):
    m=re.match(r"(?i)^(?:h)?(TM|SC)(\d+)",stem)
    if not m:raise ValueError("UNRECOGNIZED_TISSUE_STRAIN_ALIAS")
    tissue="TM" if m.group(1).upper()=="TM" else "SC"
    strain=tissue+m.group(2)
    # Three different source sample naming generations are aligned against
    # the exact sets of samples *named by timepoint* in Supplementary table.
    # These rules are source-format patterns, NOT an observed time field.
    if all("_Trypsin" in a for a in aliases):
        proposed="2Days"
        rule="TRYPSIN_NAMING_COHORT"
    elif all(re.search(r"(?i)_(?:vehicle|dex)$",a) for a in aliases):
        proposed="6Hrs"
        rule="UNDERSCORE_VEHICLE_DEX_NAMING_COHORT"
    elif all(re.search(r"-[12]$",a) for a in aliases):
        proposed="1Hr"
        rule="DASH_1_2_NAMING_COHORT"
    else:
        proposed="UNRESOLVED";rule="UNMATCHED_SAMPLE_NAME_SCHEMA"
    return tissue,strain,proposed,rule

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--supp-qc",required=True,type=Path)
    p.add_argument("--manifest",required=True,type=Path)
    p.add_argument("--biosample-qc",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    design,sha,xlsx_bytes=source_design(a.supp_qc)
    manifest=ftsv(a.manifest)
    qc=json.loads(a.biosample_qc.read_text())
    if len(manifest)!=46 or qc["pair_timepoint_metadata_verified_count"]!=0:
        raise ValueError("SOURCE_46_RUNS_OR_ZERO_DOCUMENTARY_TIMEPOINTS_DRIFT")
    pair=defaultdict(list)
    for r in manifest:
        pair[group_alias(r["sample_alias"])].append(r)
    if len(pair)!=23 or any(len(z)!=2 for z in pair.values()):
        raise ValueError("EXPECTED_23_MATCHED_PAIRED_SOURCE_LIBRARIES")
    rows=[];matched=Counter()
    for stem,rs in sorted(pair.items()):
        tissue,strain,proposed,rule=pattern_candidate(stem,[r["sample_alias"] for r in rs])
        is_named=proposed in design[tissue] and strain in design[tissue].get(proposed,[])
        status=("REVIEW_NOMINAL_STRAIN_TIMEPOINT_SET_MATCH_NO_DIRECT_RUN_TIME"
                if is_named else "HOLD_SOURCE_STRAIN_TIMEPOINT_MISMATCH")
        if is_named:matched[(tissue,proposed)]+=1
        files=sum(int(r["fastq_file_count"]) for r in rs)
        bytes_sum=sum(int(r["fastq_bytes"]) for r in rs)
        rows.append({"source_library_pair_stem":stem,
           "tissue_from_source_alias":tissue,
           "nominal_strain":strain,"candidate_exposure":proposed,
           "timepoint_proven_by_run_specific_metadata":False,
           "source_supplement_strain_has_said_timepoint":is_named,
           "confidence_tier":status,"candidate_alias_rule":rule,
           "vehicle_dex_pair_originally_verified_by_ENA":True,
           "run_accessions":";".join(sorted(r["run_accession"] for r in rs)),
           "biosample_accessions":";".join(sorted(r["sample_accession"] for r in rs)),
           "sample_aliases":";".join(sorted(r["sample_alias"] for r in rs)),
           "compressed_FASTQ_bytes":bytes_sum,
           "fastq_file_count":files,
           "original_TM_or_SC_study_scale_is_biological_strain_not_run":True,
           "not_author_confirmed_48h_strain_sample_crosswalk":True})
    if sum(matched.values())!=23:
        raise ValueError("NOT_ALL_23_CANDIDATE_PAIRS_AGREE_SOURCE_STRAIN_TIME_SET "+
                         str(dict(matched)))
    # Source supplement contains SC 6-hour experiment strains, but not a
    # single 6h SC library in this 46-run deposited accession.
    exp={("TM","2Days"):9,("TM","6Hrs"):4,("TM","1Hr"):3,
         ("SC","2Days"):4,("SC","6Hrs"):0,("SC","1Hr"):3}
    if any(matched[k]!=v for k,v in exp.items()) or len(matched)!=5:
        raise ValueError("SOURCE_ALIAS_DISTRIBUTION_UNEXPECTED "+str(dict(matched)))
    all_strains={t:sorted(set(g for gs in times.values() for g in gs)) for t,times in design.items()}
    two_day=[r for r in rows if r["tissue_from_source_alias"]=="TM" and
                           r["candidate_exposure"]=="2Days"]
    write(a.out/"MEHROTRA_23_SOURCE_CONSISTENT_TIMEPOINT_CANDIDATES.tsv",rows)
    write(a.out/"MEHROTRA_9_TM_48H_CANDIDATE_RUN_PAIRS_HOLD.tsv",two_day)
    summary={"status":"PASS_STRAIN_TIME_PATTERN_MATCH_HOLD_RUN_SPECIFIC_TIME_PROVENANCE",
         "source_archive":"PMC12882097 official EuropePMC supplement",
         "source_s003_sha256":sha,
         "source_s003_bytes":xlsx_bytes,
         "source_time_specific_marker_sheets":[TM_SHEET,SC_SHEET],
         "supplement_exact_strains_by_time":design,
         "candidate_23_source_pair_counts_by_time":{
            f"{k[0]}_{k[1]}":v for k,v in sorted(matched.items())},
         "original_supplement_SC_6h_4_strains_not_represented_by_this_46run_archive":
            design["SC"]["6Hrs"],
         "official_RUN_timepoint_fields_proven":0,
         "candidate_48h_TM_library_pairs":len(two_day),
         "candidate_48h_TM_total_fastq_bytes":sum(r["compressed_FASTQ_bytes"] for r in two_day),
         "candidate_48h_TM_total_fastq_files":sum(r["fastq_file_count"] for r in two_day),
         "author_original_count_matrix_discovered":False,
         "study_sample_specific_timepoint_crosswalk_proven":False,
         "candidate_pair_classification_is_inference_from_two_original_source_patterns":True,
         "gene_effects_supplement_already_available":True,
         "do_not_promote_to_independent_raw_48h_TM_DE_without_author_confirmation":True}
    (a.out/"MEHROTRA_SOURCE_CONSISTENT_48H_TM_CANDIDATE_GATE.json").write_text(
      json.dumps(summary,indent=2,ensure_ascii=False)+"\n")
    print("MEHROTRA_SOURCE_DESIGN",json.dumps(summary,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
