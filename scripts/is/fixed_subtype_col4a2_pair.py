#!/usr/bin/env python3
"""Recover only the previously excluded, externally prespecified rs9515201 A/C contrast.
Keep existing five-marker family and original reference exclusions; no alternate-allele fishing.
"""
import argparse,json,re,time
from pathlib import Path
import fixed_subtype_indexed as base

RS="rs9515201"
SOURCE="https://pmc.ncbi.nlm.nih.gov/articles/PMC5664302/"
def reference():
    raw,_=base.get("https://rest.ensembl.org/variation/human/"+RS+"?content-type=application/json",2000000)
    x=json.loads(raw)
    base.require(x["name"]==RS,"Wrong rsID")
    m=[v for v in x["mappings"] if v.get("assembly_name")=="GRCh38" and str(v.get("seq_region_name"))=="13" and v.get("strand")==1 and v.get("start")==v.get("end")==110388451]
    base.require(len(m)==1,"Unique primary coordinate not confirmed")
    alleles=m[0]["allele_string"].split("/")
    base.require({"A","C"}<=set(alleles) and alleles[0] in {"A","C"},"Prior published pair not in GRCh38 mapping")
    return {RS:{"status":"MATCHED_REFERENCE","chr":"13","pos":110388451,"ref":"C","alt":"A","genome_reference":alleles[0],"observed_ensembl_alleles":alleles,"ensembl_json_sha256":base.sha(raw),"prespecified_pair":"A/C","effect_orientation":"A relative to C, regardless of genomic REF","prior_source":SOURCE,"prior_table":"Table 3; 2017 lacunar stroke lead, A vs C"}}
def extract(study,ref):
    base.require(set(ref)=={RS} and ref[RS]["prespecified_pair"]=="A/C","Locked prior pair")
    base.require(ref[RS]["chr"]=="13" and ref[RS]["pos"]==110388451 and ref[RS]["ref"]=="C" and ref[RS]["alt"]=="A","Locked coordinate/allele pair")
    x=base.extract(study,ref)
    row=x.pop("rows")[RS]
    if row["status"]=="MATCHED_HARMONISED_ROW":
        # beta_alt is A-oriented because of the explicitly locked input contrast.
        row["beta_A"]=row.pop("beta_alt")
        row["A_eaf"]=row.pop("alt_eaf")
        row["effect_allele"]="A";row["other_allele"]="C"
    x["row"]=row;x["reference"]=ref[RS]
    # Reject a metadata phenotype or ancestry disagreement rather than relabeling.
    expected={"AIS":{"ischemic stroke","ischaemic stroke","any ischemic stroke","any ischaemic stroke"},
              "CES":{"cardioembolic stroke","cardioembolic ischemic stroke","cardioembolic ischaemic stroke"},
              "LAS":{"large artery stroke","large artery atherosclerotic stroke","large artery atherosclerosis stroke","large artery ischemic stroke","large artery ischaemic stroke"},
              "SVS":{"small vessel stroke","small vessel ischemic stroke","small vessel ischaemic stroke","small vessel occlusion","lacunar stroke"}}[study.split("_")[1]]
    if row["status"]=="MATCHED_HARMONISED_ROW":
        meta=x["metadata"];trait=meta.get("trait_description") or []
        if isinstance(trait,str):trait=[trait]
        norm=lambda s:re.sub(r"[-\s]+"," ",str(s).strip().lower())
        ancestries={norm(v) for s in meta.get("samples",[]) or [] for v in s.get("sample_ancestry_category",[])}
        want="european" if study.startswith("EUR") else "east asian"
        if len(trait)!=1 or norm(trait[0]) not in expected or ancestries!={want}:
            x["row"]={"status":"HOLD_METADATA_PHENOTYPE_OR_ANCESTRY","observed_trait":trait,"observed_ancestries":sorted(ancestries)}
            x["matched"]=0
    x["limits"]=["Single externally specified A/C pair in a multi-allelic reference record, not general exclusion relaxation.","No other rsID or allele substitution; preserve original fixed family of 40 contrasts.","Indexed harmonised readback only, not whole-file checksum validation or new causal inference."]
    return x
def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["reference","extract"]);p.add_argument("--study",choices=base.STUDIES);p.add_argument("--ref",type=Path);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    if a.mode=="reference":
        try:x=reference()
        except (OSError,ValueError,KeyError,TypeError) as e:x={RS:{"status":"HOLD_REFERENCE","reason":str(e)[:180]}}
    else:
        ref=json.loads(a.ref.read_text())
        if ref[RS]["status"]!="MATCHED_REFERENCE":x={"study":a.study,"accession":base.STUDIES[a.study],"row":ref[RS]}
        else:x=extract(a.study,ref)
    base.write(a.out,x)
if __name__=="__main__":main()
