#!/usr/bin/env python3
"""P8 — orthogonal study-level specificity, BioSample XML, and bounded FASTQ Range probes.

This audit uses public aggregate gene results and external repository metadata.
It uploads no raw sequence, no run/sample identifiers, and no patient data.
"""
import argparse,csv,hashlib,io,json,math,re,ssl,urllib.request,urllib.parse,xml.etree.ElementTree as ET
from pathlib import Path
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed

PROJECT="PRJEB74248"
BASE="https://www.ebi.ac.uk/ena/portal/api/filereport"
XMLBASE="https://www.ebi.ac.uk/ena/browser/api/xml/"
def http(url,cap=1_000_000,timeout=40,headers=None):
    req=urllib.request.Request(url,headers={"User-Agent":"ocular-p8-evidence-gates/1.0",**(headers or {})})
    with urllib.request.urlopen(req,timeout=timeout) as r:
        data=r.read(cap+1)
        status=getattr(r,"status",None)
        hdr={x:r.headers.get(x) for x in ("Content-Type","Content-Range","Content-Length","Accept-Ranges")}
    if len(data)>cap:raise ValueError("REMOTE_RESPONSE_TOO_LARGE")
    return data,status,hdr
def enarun(fields):
    query=urllib.parse.urlencode({"accession":PROJECT,"result":"read_run","fields":",".join(fields),"format":"tsv"})
    raw,_,_=http(BASE+"?"+query,cap=1_800_000,timeout=75)
    rd=csv.DictReader(io.StringIO(raw.decode("utf8","replace")),delimiter="\t")
    if not set(fields).issubset(rd.fieldnames or []):
        raise ValueError("ENA_FIELD_REPORT_SCHEMA_CHANGED")
    rows=list(rd)
    if len(rows)!=62 or len({r["run_accession"] for r in rows})!=62:
        raise ValueError("SHOTGUN_RUN_COUNT_DRIFT")
    return rows,hashlib.sha256(raw).hexdigest()
def write(out,filename,d):
    (out/filename).write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf8")
def tsv(out,filename,rows):
    if not rows:raise ValueError("EMPTY_PUBLIC_SUMMARY")
    with (out/filename).open("w",encoding="utf8",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        wr.writeheader();wr.writerows(rows)
def safe_num(x):
    n=float(x)
    if not math.isfinite(n):raise ValueError("NON_FINITE_PUBLIC_NUMBER")
    return n
def proportion_ci(success,total,z=1.959963984540054):
    p=success/total
    denom=1+z*z/total
    ctr=(p+z*z/(2*total))/denom
    half=z/denom*math.sqrt(p*(1-p)/total+z*z/(4*total*total))
    return [round(ctr-half,6),round(ctr+half,6)]

def emt_compare(out,genes_file,json_file):
    from scipy.stats import fisher_exact
    with genes_file.open(encoding="utf8",newline="") as f:
        genes=list(csv.DictReader(f,delimiter="\t"))
    summary=json.loads(json_file.read_text(encoding="utf8"))
    if len(genes)!=194 or summary.get("EMT_genes")!=194:
        raise ValueError("EXPECTED_P7A_194_EMT_ROWS")
    same=[r for r in genes if r["sign_agrees"]=="True"]
    if len(same)!=141 or summary.get("EMT_same_sign_count")!=141:
        raise ValueError("P7A_141_OF_194_CHANGED")
    nbg=summary["background_nonEMT_n"]
    yesbg=summary["background_nonEMT_same_sign_count"]
    if nbg!=12056 or yesbg!=6810:raise ValueError("P7A_BACKGROUND_DRIFT")
    odds,p=fisher_exact([[141,53],[yesbg,nbg-yesbg]],alternative="greater")
    strong=[x for x in genes if x["both_study_FDR_lt_0_05"]=="True"]
    matched=[x for x in strong if x["sign_agrees"]=="True"]
    if len(strong)!=49 or len(matched)!=46:
        raise ValueError("P7A_49_OF_46_STRONG_EMT_DRIFT")
    threshold_rows=[]
    for threshold in (0,0.25,0.5,0.75,1.0):
        kept=[x for x in genes if min(abs(safe_num(x["primary_500nM_log2fc"])),
          abs(safe_num(x["external_100nM_log2fc"])))>=threshold]
        agree=sum(x["sign_agrees"]=="True" for x in kept)
        both=sum(x["both_study_FDR_lt_0_05"]=="True" for x in kept)
        both_agree=sum(x["both_study_FDR_lt_0_05"]=="True" and x["sign_agrees"]=="True" for x in kept)
        threshold_rows.append({"min_abs_effect_both_cohorts":threshold,
          "EMT_genes_passing":len(kept),"matching_sign":agree,
          "matching_sign_fraction":round(agree/len(kept),6) if kept else "",
          "both_cohorts_significant":both,"both_cohorts_significant_and_same":both_agree,
          "threshold_exploratory_posthoc":True})
    tsv(out,"P8a_EMT_shared_gene_effect_threshold_sensitivity.tsv",threshold_rows)
    targets=[]
    for x in strong:
        targets.append({"gene_symbol":x["gene_symbol"],"gene_id":x["gene_id"],
          "primary_log2fc":safe_num(x["primary_500nM_log2fc"]),
          "external_log2fc":safe_num(x["external_100nM_log2fc"]),
          "primary_FDR":safe_num(x["primary_FDR"]),
          "external_FDR":safe_num(x["external_FDR"]),
          "sign_agrees":x["sign_agrees"]=="True",
          "is_one_of_12_original_discordant":x["in_original_12_targets"]=="True",
          "not_independently_verified_at_protein_level":True})
    tsv(out,"P8a_49_double_FDR_significant_EMT_genes.tsv",
        sorted(targets,key=lambda r:(not r["sign_agrees"],-min(abs(r["primary_log2fc"]),abs(r["external_log2fc"])),r["gene_symbol"])))
    outcome={"status":"PASS_EMT_VS_MATCHED_BACKGROUND_DESCRIPTIVE_ASSOCIATION",
      "p7a_source_run":41179413138,"n_EMT":194,"n_EMT_agree":141,
      "n_other":nbg,"n_other_agree":yesbg,
      "EMT_sign_concordance":141/194,"nonEMT_sign_concordance":yesbg/nbg,
      "absolute_percentage_point_advantage":round(100*(141/194-yesbg/nbg),3),
      "odds_ratio_vs_nonEMT":round(odds,6),
      "gene_level_fisher_one_sided_p":p,
      "EMT_wilson_95_interval":proportion_ci(141,194),
      "nonEMT_wilson_95_interval":proportion_ci(yesbg,nbg),
      "EMT_49_both_study_FDR_significant":len(strong),
      "EMT_49_both_study_FDR_significant_same_sign":len(matched),
      "null_gene_independence_assumption_is_unreliable":True,
      "fisher_p_not_valid_pathway_replication_p":True,
      "no_matched_dose_or_patient_IOP":True,
      "note":"Fisher treats genes as independent units despite co-expression and gene-set curation; descriptive ranking/CI must NOT serve as formal independent biological replication or pathway significance."}
    write(out,"P8a_emt_background_specificity_gate.json",outcome)
    print("P8A",json.dumps(outcome,ensure_ascii=False),flush=True)

def xmlprobe(out):
    rows,manifest_sha=enarun(["run_accession","sample_accession"])
    samples=sorted({r["sample_accession"] for r in rows})
    if len(samples)!=62:raise ValueError("P7B_EXPECTED_62_UNIQUE_SAMPLES")
    # No raw sample names, values or accession strings in uploaded artifact.
    def one(sample):
        try:
            blob,_,_=http(XMLBASE+urllib.parse.quote(sample,safe=""),cap=500_000,timeout=30)
            root=ET.fromstring(blob)
            attributes=[]
            for node in root.findall(".//SAMPLE_ATTRIBUTE"):
                tag=node.findtext("TAG","").strip()
                values=[v.text.strip() for v in node.findall("VALUE") if v.text]
                if tag:attributes.append((tag,values))
            categorytags=[]
            categories=[]
            blank_match=False
            for tag,vals in attributes:
                key=tag.lower()
                if any(k in key for k in ("disease","phenotype","health status","health state","diagnos","condition","treatment","medication","control")):
                    categorytags.append(tag)
                    textvals=" ".join(vals).lower()
                    # Explicit attribute value, no inference from shared titles.
                    diseased=bool(re.search(r"\b(?:glaucoma|poag|pexg)\b",textvals))
                    healthy=bool(re.search(r"\b(?:healthy|normal|control)\b",textvals))
                    if diseased and not healthy:categories.append("CASE")
                    elif healthy and not diseased:categories.append("HEALTHY")
                    elif diseased and healthy:categories.append("AMBIGUOUS")
                if "blank" in key or "negative control" in key:
                    blank_match=True
            category=set(categories)
            if category=={"CASE"}:assigned="CASE"
            elif category=={"HEALTHY"}:assigned="HEALTHY"
            elif not category:assigned="UNKNOWN"
            else:assigned="AMBIGUOUS"
            return {"ok":True,"tags":categorytags,"category":assigned,
              "blank_attribute":blank_match,"xml_bytes":len(blob),
              "attribute_count":len(attributes)}
        except Exception as e:
            return {"ok":False,"error":type(e).__name__}
    details=[]
    with ThreadPoolExecutor(max_workers=8) as ex:
        fs={ex.submit(one,s):s for s in samples}
        for fut in as_completed(fs):
            details.append(fut.result())
    passed=[r for r in details if r["ok"]]
    fails=Counter(r["error"] for r in details if not r["ok"])
    cats=Counter(r["category"] for r in passed)
    tags=Counter(t for r in passed for t in set(r["tags"]))
    good=(len(passed)==62)
    result={"status":"PASS_62_ENA_BIOSAMPLE_XML_ATTR_AUDIT" if good else "HOLD_PARTIAL_ENA_XML_QUALITY",
      "source":PROJECT,"expected_unique_samples":62,
      "sample_xml_download_success":len(passed),
      "sample_xml_errors_by_type":dict(fails),
      "explicit_clinical_attribute_categories":dict(cats),
      "attribute_keys_containing_clinical_terms":dict(tags.most_common(25)),
      "explicit_control_or_blank_attribute_count":sum(r["blank_attribute"] for r in passed),
      "zero_cases_in_xml_does_not_prove_no_case_labels_elsewhere":True,
      "patient_to_run_crosswalk_verified":False,
      "glaucoma_group_can_be_assigned_from_sufficient_explicit_xml_fields":
       good and cats["CASE"]>0 and cats["HEALTHY"]>0 and cats["AMBIGUOUS"]==0 and cats["UNKNOWN"]==0,
      "negative_blank_sequencing_read_validated":False,
      "raw_sample_values_or_ids_published":False,
      "run_sample_manifest_sha256":manifest_sha,
      "interpretation":"Explicit XML sample attributes audited, no sample data published. If labels are missing, hold differential microbiome analysis; if labels present, validate donor/eye units and clinical covariates before inference."}
    write(out,"P8b_ENA_62_BioSample_XML_clinical_field_gate.json",result)
    print("P8B",json.dumps({k:result[k] for k in ("status","sample_xml_download_success","explicit_clinical_attribute_categories","glaucoma_group_can_be_assigned_from_sufficient_explicit_xml_fields")},ensure_ascii=False),flush=True)

def rangeprobe(out):
    fields=["run_accession","fastq_ftp","submitted_ftp","submitted_bytes","submitted_md5"]
    rows,sha=enarun(fields)
    a=sorted((r for r in rows if r["fastq_ftp"] and r["submitted_ftp"]),key=lambda x:x["run_accession"])
    b=sorted((r for r in rows if not r["fastq_ftp"] and r["submitted_ftp"]),key=lambda x:x["run_accession"])
    if len(a)!=31 or len(b)!=31:raise ValueError("P7C_31_AND_31_CATEGORIES_DRIFT")
    picks=[("fastq_ftp_also_available",a[0]),("submitted_ftp_only",b[0])]
    result=[]
    for category,run in picks:
        filename=run["submitted_ftp"].split(";")[0]
        if not filename.startswith(("ftp.sra.ebi.ac.uk/","ftp.ebi.ac.uk/")):
            result.append({"category":category,"status":"HOLD_UNEXPECTED_ENA_ARCHIVE_HOST",
              "downloaded_sequence_bytes":0,"raw_sequence_published":False})
            continue
        url="https://"+filename
        # Remote hosts may ignore Range. We read at most 64 KiB and close; never persist reads.
        try:
            req=urllib.request.Request(url,headers={"Range":"bytes=0-65535","User-Agent":"ocular-p8-gzip-range-QC/1.0"})
            with urllib.request.urlopen(req,timeout=35) as response:
                status=response.status
                headers={"content_range":response.headers.get("Content-Range") or "",
                         "content_type":response.headers.get("Content-Type") or ""}
                data=response.read(65536)
            magic=data.startswith(b"\x1f\x8b")
            result.append({"category":category,"status":"PASS_BOUNDED_GZIP_MAGIC" if magic else "HOLD_NOT_GZIP_MAGIC",
                "http_status":status,"range_header_returned":bool(headers["content_range"]),
                "gzip_magic":magic,"bytes_read_max_65536":len(data),
                "raw_sequence_published":False,"checksum_verified":False,
                "complete_fastq_file_downloaded":False})
        except Exception as e:
            result.append({"category":category,"status":"HOLD_HTTP_READ_FAILED",
                "error_type":type(e).__name__,"bytes_read_max_65536":0,
                "raw_sequence_published":False,"checksum_verified":False})
    result_summary={"status":"PASS_BOUNDED_PUBLIC_FILE_PROBES_AUDITED_NOT_WHOLE_FILE_VERIFIED",
       "accession":PROJECT,"study_total_submitted_runs":62,"study_total_submitted_fastq_files":155,
       "max_content_read_per_probe":65536,"probes":result,
       "complete_file_checksum_validation_performed":False,
       "paired_fastq_record_integrity_tested":False,
       "host_read_filtering_tested":False,"clinical_case_control_verified":False,
       "metadata_sha256":sha,"no_fastq_sequences_saved_or_uploaded":True}
    write(out,"P8c_two_bounded_shotgun_fastq_readability_probes.json",result_summary)
    print("P8C",json.dumps(result_summary,ensure_ascii=False),flush=True)
def main():
    a=argparse.ArgumentParser()
    a.add_argument("--mode",required=True,choices=("emt","biosample","range"))
    a.add_argument("--p7-emt",type=Path)
    a.add_argument("--p7-summary",type=Path)
    a.add_argument("--out",required=True,type=Path)
    v=a.parse_args();v.out.mkdir(parents=True,exist_ok=True)
    if v.mode=="emt":
        if not v.p7_emt or not v.p7_summary:a.error("emt requires --p7-emt and --p7-summary")
        emt_compare(v.out,v.p7_emt,v.p7_summary)
    elif v.mode=="biosample":xmlprobe(v.out)
    else:rangeprobe(v.out)
if __name__=="__main__":main()
