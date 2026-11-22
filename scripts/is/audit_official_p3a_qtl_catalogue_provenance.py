#!/usr/bin/env python3
"""Verify 3 official eQTL Catalogue GE study identities for all 36 P3A assays.

Science HOLD: eQTL Catalogue nominal cis records are GRCh38 (ALT-effect).
The GIGASTROKE P3A GWAS is GRCh37. Dataset sample size is not automatically
effective per-gene N. No direct SNP matching, liftover, MR, or Coloc here.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path

PINNED_COMMIT="07918a267afdc1c84b5444bea9bd902e7201728f"
FILES={
 "paths":("tabix/tabix_ftp_paths.tsv","2a5aea0b8fe7f75ace257dc3b4ff59ec18b97a68"),
 "meta":("data_tables/dataset_metadata_r7.tsv","df3163482aca99c481840caf47691f6da4acb76d"),
 "columns":("tabix/Columns.md","b3fafeaaf4b8faff9f2b4ad905b570fa74049b25"),
}
EXPECTED={
 "QTD000131":{"study_id":"QTS000015","study_label":"GTEx","sample_group":"artery_aorta","sample_size":"387","quant_method":"ge","tissue_id":"UBERON_0001496","pmid":"32913098","group_n":10},
 "QTD000171":{"study_id":"QTS000015","study_label":"GTEx","sample_group":"brain_cortex","sample_size":"205","quant_method":"ge","tissue_id":"UBERON_0001870","pmid":"32913098","group_n":11},
 "QTD000434":{"study_id":"QTS000025","study_label":"ROSMAP","sample_group":"brain_naive","sample_size":"560","quant_method":"ge","tissue_id":"UBERON_0009834","pmid":"28869584","group_n":15},
}
BASE="https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/"


def require(expr:bool,reason:str)->None:
    if not expr: raise ValueError(reason)


def git_blob(raw:bytes)->str:
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()


def get_official_sources(fetch=urllib.request.urlopen)->dict:
    found={}
    for label,(path,sha) in FILES.items():
        url=BASE+PINNED_COMMIT+"/"+path
        req=urllib.request.Request(url,headers={"User-Agent":"IS-P3A-QTL-Source-Audit/1.0"})
        with fetch(req,timeout=35) as response:
            raw=response.read(950_000)
        require(git_blob(raw)==sha,"pinned official metadata file changed: "+path)
        found[label]=raw.decode("utf-8-sig")
    return found


def table(raw:str)->dict:
    rows=list(csv.DictReader(io.StringIO(raw),delimiter="\t"))
    require(rows and all(None not in x for x in rows),"invalid official TSV schema")
    return {x["dataset_id"]:x for x in rows}


def validate(manifest:dict,sources:dict)->dict:
    require(manifest.get("schema")=="is-p3a-36-real-group-keys-from-legacy-393-v1"
            and manifest.get("origin",{}).get("commit_sha")==
                 "0bf646fd6f29322d17da878253f37def24c01c74"
            and manifest.get("scientific_status")=="HOLD","changed archived QTL registry")
    groups=manifest.get("groups",[])
    require(len(groups)==36 and len({(g["region_id"],g["QTL_dataset"],g["gene_id"])
                                      for g in groups})==36
            and len({g["gene_id"] for g in groups})==20
            and all(g["region_id"]=="IS_XDATA_G0015" and
                    g["approved_for_new_coloc"] is False for g in groups),
            "incorrect original 36-group manifest")

    expected_counts={k:v["group_n"] for k,v in EXPECTED.items()}
    require(Counter(g["QTL_dataset"] for g in groups)==expected_counts,
            "unexpected QTL dataset assignment to original candidate groups")
    require(all(x in sources for x in FILES),"official catalogue metadata not loaded")
    columns=sources["columns"]
    for fragment in ("Based on GRCh38 coordinates","GRCh38 reference allele",
                     "alt - GRCh38 alternative allele (also the effect allele)",
                     "beta - Regression coefficient", "se - Standard error",
                     "maf - Minor allele frequency"):
        require(fragment in columns,"official GRCh38 QTL columns schema absent: "+fragment)
    meta=table(sources["meta"])
    paths=table(sources["paths"])
    assessments={}
    by_group=[]
    for qtl,exp in EXPECTED.items():
        require(qtl in meta and qtl in paths,"dataset missing in official Release 7: "+qtl)
        a,b=meta[qtl],paths[qtl]
        for key in ("study_id","dataset_id","study_label","sample_group","sample_size",
                    "quant_method","tissue_id","tissue_label","condition_label"):
            require(a.get(key)==b.get(key),"official tabix/Release 7 metadata conflict: "+qtl+" "+key)
        for key in ("study_id","study_label","sample_group","sample_size",
                    "quant_method","tissue_id","pmid"):
            require(a.get(key)==exp[key],"unexpected study/tissue/sample metadata: "+qtl+" "+key)
        require(a.get("study_type")=="bulk","QTL type differs from known bulk study")
        ftp=b["ftp_path"]
        canonical=f"ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/{exp['study_id']}/{qtl}/{qtl}.all.tsv.gz"
        require(ftp==canonical,"not the official full nominal gene expression QTL: "+qtl)
        assessments[qtl]={
           "study_id":a["study_id"],"study_label":a["study_label"],
           "tissue":a["tissue_label"],"sample_group":a["sample_group"],
           "nominal_catalogue_sample_size":int(a["sample_size"]),
           "sample_size_role":"QTL_DATASET_METADATA_NOT_VALIDATED_PER_GENE_EFFECTIVE_N",
           "source_original_ftp":ftp,
           "source_https_probe":ftp.replace("ftp://ftp.ebi.ac.uk/","https://ftp.ebi.ac.uk/"),
           "study_pmid":a["pmid"],
           "QTL_effect_allele":"ALT_GRCh38",
           "QTL_build":"GRCh38",
           "GWAS_build":"GRCh37",
           "build_liftover_required":True,
           "source_file_variant_rows_verified":False,
           "cis_full_variants_gene_matched":False,
           "official_individual_N_approved":False,
           "scientific_status":"HOLD_BUILD_AND_PER_GENE_SOURCE_N",
        }
    for g in sorted(groups,key=lambda x:(x["QTL_dataset"],x["gene_id"])):
        a=assessments[g["QTL_dataset"]]
        by_group.append({
            "region_id":g["region_id"],
            "QTL_dataset":g["QTL_dataset"],
            "gene_id":g["gene_id"],"gene_symbol":g["gene_symbol"],
            "study_label":a["study_label"],"tissue":a["tissue"],
            "catalogue_study_n":a["nominal_catalogue_sample_size"],
            "gene_effective_n_approved":False,
            "eQTL_nominal_fullcis_QTL_GRCh38_source_proven":False,
            "GRCh38_to_GRCh37_liftover_independently_verified":False,
            "GRCh37_GWAS_signed_alleles_matched":False,
            "coloc_abf_eligible":False,"coloc_susie_eligible":False,
            "status":"HOLD_QTL_FULL_CIS_SOURCE_BUILD_GWAS_N",
        })
    require(len(by_group)==36 and len({x["gene_id"] for x in by_group})==20,
            "QTL group coverage changed")
    return {
       "schema":"is-p3a-official-qtl-study-metadata-v1",
       "technical_status":"PASS_3_SOURCE_IDENTITIES_IN_IMMUTABLE_CATALOGUE",
       "scientific_status":"HOLD",
       "official_source_commit":PINNED_COMMIT,
       "official_source_git_blob_sha1":{k:v[1] for k,v in FILES.items()},
       "QTL_datasets":assessments,
       "n_original_assay_groups":36,
       "n_distinct_genes":20,
       "n_studies":2,
       "n_tissue_datasets":3,
       "metadata_sample_size_automatically_approved_per_gene":False,
       "qtl_reference_build":"GRCh38",
       "gwas_reference_build":"GRCh37",
       "liftover_complete":False,
       "coloc_abf_eligible":0,"coloc_abf_executed":0,
       "coloc_susie_eligible":0,"coloc_susie_executed":0,
       "group_rows":by_group,
    }


def probe_head(result:dict,fetch=urllib.request.urlopen)->None:
    # Three sequential HEAD requests max. Never download the large sumstats files.
    for qtl,row in result["QTL_datasets"].items():
        url=row["source_https_probe"]
        req=urllib.request.Request(url,method="HEAD",
                    headers={"User-Agent":"IS-Masters-QTL-Catalogue-ReadOnly/1.0"})
        try:
            with fetch(req,timeout=24) as r:
                status=int(getattr(r,"status",0))
                length=r.headers.get("Content-Length")
                row["source_HEAD_probe"]={
                    "http_status":status,"content_length_header":length,
                    "status":"REMOTE_FILE_REACHABLE_NOT_DOWNLOADED"
                              if status==200 else "HOLD_HTTP_STATUS"}
        except (urllib.error.URLError,TimeoutError,OSError,ValueError) as exc:
            row["source_HEAD_probe"]={
                "status":"HOLD_NETWORK_OR_REMOTE_SOURCE",
                "error_type":type(exc).__name__}
        row["source_file_variant_rows_verified"]=False
    result["any_source_fullcis_associations_downloaded"]=False
    result["scientific_status"]="HOLD"


def main()->int:
    parser=argparse.ArgumentParser()
    parser.add_argument("--manifest",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--probe-head",action="store_true")
    args=parser.parse_args()
    try:
        manifest=json.loads(args.manifest.read_text(encoding="utf-8"))
        result=validate(manifest,get_official_sources())
        if args.probe_head:probe_head(result)
    except (OSError,ValueError,KeyError,TypeError,urllib.error.URLError) as exc:
        result={"schema":"is-p3a-official-qtl-study-metadata-v1",
          "technical_status":"HOLD_ORIGINAL_OFFICIAL_METADATA_UNVERIFIED",
          "scientific_status":"HOLD","coloc_abf_eligible":0,
          "coloc_abf_executed":0,"error":str(exc)[:280]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:result.get(k) for k in
      ("technical_status","scientific_status","n_original_assay_groups",
       "n_tissue_datasets","qtl_reference_build","gwas_reference_build","error")},
      ensure_ascii=False))
    return 0 if result["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
