#!/usr/bin/env python3
"""Only 12 public BioSample XML records behind six ambiguous -1/-2 pairs.

Do not infer sample-level steroid timepoint or control/dex from position.
No FASTQ, counts, human genotype or identifiable patient records requested.
"""
import argparse,csv,hashlib,io,json,re,urllib.request,xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict,Counter

FIELDS=("treatment","drug","culture","time","exposure","condition","control",
        "source","sample","tissue","cell","strain","library")
def direct(url):
    req=urllib.request.Request(url,headers={"User-Agent":"ocular-public-biosample-crosswalk/1.0"})
    with urllib.request.urlopen(req,timeout=23) as f:
        data=f.read(650_000)
    if not data or len(data)>600_000:raise ValueError("EMPTY_OR_EXCESSIVE_XML")
    return data
def candidate_treatment(snippets):
    string=" ".join(snippets).lower()
    dex=bool(re.search(r"\b(dexamethasone|dex)\b",string))
    veh=bool(re.search(r"\b(vehicle|ethanol control|untreated control)\b",string))
    return "DEX_EXPLICIT" if dex and not veh else "VEHICLE_EXPLICIT" if veh and not dex else "UNKNOWN_OR_MIXED"
def target_body(root):
    attrs=[]
    for a in root.findall(".//SAMPLE_ATTRIBUTE"):
        tag=(a.findtext("TAG") or "").strip()
        value=(a.findtext("VALUE") or "").strip()
        if any(x in tag.lower() for x in FIELDS):
            attrs.append({"tag":tag[:80],"value":value[:140]})
    titles=[e.text.strip()[:160] for e in root.findall(".//TITLE") if e.text and e.text.strip()]
    descriptions=[e.text.strip()[:160] for e in root.findall(".//DESCRIPTION") if e.text and e.text.strip()]
    return attrs,titles[:2],descriptions[:2]
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    with args.manifest.open(encoding="utf8",newline="") as f:
        raw=list(csv.DictReader(f,delimiter="\t"))
    ambiguous=[r for r in raw if re.search(r"-[12]$",r["sample_alias"])]
    if len(raw)!=46 or len(ambiguous)!=12:
        raise ValueError("SOURCE_SRA_46RUN_OR_12_UNASSIGNED_DRIFT")
    rows=[];fetched=0;failed=0
    for r in ambiguous:
        accession=r["sample_accession"]
        if not re.fullmatch(r"SAM[END][A-Z]?[0-9]+",accession):
            raise ValueError("INVALID_BIOSAMPLE_ACCESSION")
        url="https://www.ebi.ac.uk/ena/browser/api/xml/"+accession
        try:
            content=direct(url);root=ET.fromstring(content)
            attrs,titles,descs=target_body(root)
            # Assign treatment only from explicit sample-local label/attribute;
            # do not derive -1=vehicle, -2=DEX or infer time from publication.
            data=[x["value"] for x in attrs if any(y in x["tag"].lower() for y in
                    ("treatment","drug","exposure","condition","control"))]
            tr=candidate_treatment(data)
            times=[]
            for x in attrs:
                if any(y in x["tag"].lower() for y in ("time","exposure","duration")):
                    times.append(x["tag"]+"="+x["value"])
            out={"run_accession":r["run_accession"],
                 "biosample_accession":accession,
                 "sample_alias":r["sample_alias"],
                 "status":"FETCHED",
                 "biosample_xml_sha256":hashlib.sha256(content).hexdigest(),
                 "sample_title":";".join(titles),
                 "sample_description":";".join(descs),
                 "relevant_attributes":attrs,
                 "explicit_treatment_from_BioSample_attributes":tr,
                 "explicit_exposure_duration_attributes":times[:5]}
            fetched+=1
        except Exception as e:
            out={"run_accession":r["run_accession"],
                 "biosample_accession":accession,
                 "sample_alias":r["sample_alias"],
                 "status":"HOLD_XML_UNAVAILABLE",
                 "reason":type(e).__name__+" "+str(e)[:160]}
            failed+=1
        rows.append(out)
    assignments=Counter(x.get("explicit_treatment_from_BioSample_attributes","UNAVAILABLE") for x in rows)
    summary={"status":"GO_PUBLIC_BIOSAMPLE_PROVENANCE_ONLY" if fetched==12 else "HOLD_INCOMPLETE_BIOSAMPLE_XML",
             "project":"PRJNA1243143",
             "ambiguous_library_files_examined":12,
             "fetched":fetched,"unavailable":failed,
             "metadata_attribute_based_treatment_assignments":dict(assignments),
             "sample_timepoint_assignments_verified_from_attributes":
                 sum(bool(x.get("explicit_exposure_duration_attributes")) for x in rows),
             "file_order_not_a_valid_Dex_vehicle_label":True,
             "author_original_sample_sheet_still_required_if_ambiguous":True,
             "not_a_raw_FASTQ_or_count_analysis":True,
             "not_12_independent_human_donors":True}
    (args.out/"MEHROTRA_SRA_12_AMBIGUOUS_BIOSAMPLES.json").write_text(
        json.dumps(rows,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    (args.out/"MEHROTRA_SRA_12_BIOSAMPLE_MAPPING_QC.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("AMBIGUOUS_BIOSAMPLE_QC",json.dumps(summary,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
