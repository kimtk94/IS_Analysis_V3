#!/usr/bin/env python3
"""Two-way GRCh38 -> GRCh37 exact SNP locus map with source REF allele check.

Reads public one-variant original TPMI GWAS extraction; uses official
Ensembl REST assembly map, reverse-map and each reference assembly base.
Records HOLD if mapping unsupported/ambiguous or reference bases differ.
No source GWAS re-download and no GIGASTROKE outcome effect claim.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import urllib.request
from urllib.parse import urlsplit

GRCH38="rest.ensembl.org"
GRCH37="grch37.rest.ensembl.org"
SHA="ce14daef23f9fadf921030d9ed2429383d927f6da391c620d9abf870045b2e62"
COMPLEMENT=str.maketrans("ACGT","TGCA")
def complement(base):
    return base.translate(COMPLEMENT)

def source_gate(source):
    if (source.get("schema")!="ckm-TPMI43321-pinned-source-genomewide-significant-lead-v1"
        or source.get("technical_status")!="PASS_PINNED_1739048_SOURCE_ROWS_ONE_SIGNIFICANT_VARIANT"
        or source.get("pinned_source_sha256_verified")!=SHA
        or source.get("all_original_rows_verified")!=1739048
        or source.get("source_genomewide_significant_count")!=1
        or source.get("original_source_genome_build")!="GRCh38"
        or source.get("scientific_status")!="HOLD"
        or source.get("ready_for_ldsc") is not False):
        raise ValueError("unverified published TPMI one-SNP original source")
    rows=source.get("original_significant_associations")
    if not isinstance(rows,list) or len(rows)!=1:
        raise ValueError("one original genomewide hit not provided")
    item=rows[0]
    if not(str(item.get("chrom_GRCh38")).isdecimal() and
           1<=int(item["chrom_GRCh38"])<=22 and
           type(item.get("position_GRCh38")) is int and item["position_GRCh38"]>0 and
           item.get("A1_effect_allele") in (item.get("REF_GRCh38"),item.get("ALT_GRCh38"))):
        raise ValueError("invalid source chrom/coordinate/source-effect allele")
    return item

def official_ensembl(url,host,opener=urllib.request.urlopen):
    p=urlsplit(url)
    if p.scheme!="https" or p.hostname!=host:
        raise ValueError("unexpected official Ensembl source host")
    req=urllib.request.Request(url,headers={"Content-Type":"application/json",
            "User-Agent":"CKM-Stroke-TPMI-Lead-Build-Map/1.0"})
    with opener(req,timeout=45) as resp:
        redirect=urlsplit(resp.geturl())
        if resp.status!=200 or redirect.scheme!="https" or redirect.hostname!=host:
            raise ValueError("official Ensembl REST redirect/status changed")
        body=resp.read(40_001)
        if len(body)>40_000:raise ValueError("oversized Ensembl single-variant mapping")
    return json.loads(body.decode("utf-8","strict"))

def parse_map(data,source_chrom,source_pos,from_build,to_build):
    arr=data.get("mappings")
    if not isinstance(arr,list) or len(arr)!=1:
        raise ValueError("no unique liftover between genome assemblies")
    mapping=arr[0]
    ori=mapping.get("original",{});target=mapping.get("mapped",{})
    if (ori.get("assembly")!=from_build or target.get("assembly")!=to_build
        or str(ori.get("seq_region_name"))!=str(source_chrom)
        or ori.get("start")!=source_pos or ori.get("end")!=source_pos
        or target.get("start")!=target.get("end")
        or type(target.get("start")) is not int or target["start"]<=0
        or target.get("strand") not in (-1,1)
        or ori.get("strand") not in (-1,1)
        or not str(target.get("seq_region_name","")).isdecimal()
        or not 1<=int(target.get("seq_region_name"))<=22):
        raise ValueError("assembly map source/target precision or assembly drift")
    return {"source_build":from_build,
            "source_chrom":str(source_chrom),"source_pos":source_pos,
            "target_build":to_build,"target_chrom":str(target["seq_region_name"]),
            "target_pos":target["start"],
            "mapping_orientation":target["strand"]*ori["strand"]}

def validate_roundtrip(fwd,back):
    if (fwd["source_build"]!="GRCh38" or fwd["target_build"]!="GRCh37"
        or back["source_build"]!="GRCh37" or back["target_build"]!="GRCh38"
        or fwd["target_chrom"]!=back["source_chrom"]
        or fwd["target_pos"]!=back["source_pos"]
        or fwd["source_chrom"]!=back["target_chrom"]
        or fwd["source_pos"]!=back["target_pos"]
        or fwd["mapping_orientation"]*back["mapping_orientation"]!=1):
        raise ValueError("GRCh38 -> GRCh37 -> GRCh38 roundtrip not unique/consistent")

def check_alleles(hit,fwd,ref38,ref37):
    if len(hit["REF_GRCh38"])!=1 or len(hit["ALT_GRCh38"])!=1:
        raise ValueError("not a SNV; normalized indel liftover requires separate allele context")
    if hit["REF_GRCh38"]==hit["ALT_GRCh38"]:
        raise ValueError("original source REF equals ALT")
    if ref38.get("seq")!=hit["REF_GRCh38"]:
        raise ValueError("GRCh38 original REF does not match original Ensembl reference")
    trans=complement if fwd["mapping_orientation"]==-1 else (lambda x:x)
    ori_ref=trans(hit["REF_GRCh38"])
    ori_alt=trans(hit["ALT_GRCh38"])
    a1=trans(hit["A1_effect_allele"])
    if ref37.get("seq")!=ori_ref:
        raise ValueError("GRCh37 mapped reference base disagrees after strand conversion")
    return {
      "chrom_GRCh37":fwd["target_chrom"],"pos_GRCh37":fwd["target_pos"],
      "REF_GRCh37":ori_ref,"ALT_GRCh37":ori_alt,
      "A1_effect_allele_GRCh37":a1,
      "A1_effect_direction_logOR_preserved_under_reference_orientation":True,
      "allele_source_GRCh38_reference_verified":True,
      "mapped_GRCh37_reference_verified":True,
      "result_not_yet_harmonized_with_GIGASTROKE_beta_OR_SE":True
    }

def source_to_reference(source,reader=official_ensembl):
    hit=source_gate(source)
    chrom=hit["chrom_GRCh38"];pos=hit["position_GRCh38"]
    url38=f"https://{GRCH38}/map/human/GRCh38/{chrom}:{pos}..{pos}:1/GRCh37?content-type=application/json"
    fwd=parse_map(reader(url38,GRCH38),chrom,pos,"GRCh38","GRCh37")
    chrom37=fwd["target_chrom"];p37=fwd["target_pos"]
    url37=f"https://{GRCH37}/map/human/GRCh37/{chrom37}:{p37}..{p37}:1/GRCh38?content-type=application/json"
    back=parse_map(reader(url37,GRCH37),chrom37,p37,"GRCh37","GRCh38")
    validate_roundtrip(fwd,back)
    source38=reader(f"https://{GRCH38}/sequence/region/human/{chrom}:{pos}..{pos}:1?content-type=application/json",GRCH38)
    source37=reader(f"https://{GRCH37}/sequence/region/human/{chrom37}:{p37}..{p37}:1?content-type=application/json",GRCH37)
    attrs=check_alleles(hit,fwd,source38,source37)
    return {
       "schema":"ckm-TPMI43321-lead-variant-Ensembl-GRCh38-GRCh37-v1",
       "technical_status":"PASS_BIDIRECTIONAL_ENSEMBL_BUILD_MAP_AND_REF_ALLELE" ,
       "original_source_SHA256_verified":SHA,
       "lead_GRCh38":hit,
       "mapped_GRCh37":attrs,
       "liftover_first_forward":fwd,
       "liftover_reverse_verified":True,
       "reference_sources":[url38,url37],
       "TPMI_original_2025_PLINK2_is_not_2026_Neurology_paper":True,
       "GIGASTROKE_EUR_source_effect_at_target_observed":False,
       "scientific_replication_GO":False,
       "case_control_Neff_available":False,
       "scientific_status":"HOLD","ready_for_ldsc":False,
       "LDSC_h2_rg_MR_coloc_executed":False,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--lead",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    try:
        result=source_to_reference(json.loads(a.lead.read_text("utf-8")))
    except Exception as e:
        result={"technical_status":"HOLD_ENSEMBL_GRCH37_38_VARIANT_MAPPING_UNVERIFIED",
            "error_type":type(e).__name__,"error_detail":str(e)[:280],
            "scientific_status":"HOLD","ready_for_ldsc":False,
            "LDSC_h2_rg_MR_coloc_executed":False}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n","utf-8")
    print(json.dumps({"status":result["technical_status"],
        "mapped":result.get("mapped_GRCh37"),
        "error":result.get("error_detail"),
        "scientific":result["scientific_status"]}))
    if result["technical_status"]!="PASS_BIDIRECTIONAL_ENSEMBL_BUILD_MAP_AND_REF_ALLELE":
        raise SystemExit(1)
if __name__=="__main__":
    main()
