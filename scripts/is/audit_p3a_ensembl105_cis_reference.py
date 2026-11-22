#!/usr/bin/env python3
"""Compare all 36 authentic P3A QTL assay keys with a pinned official
eQTL Catalogue phenotype-position reference (Ensembl release 105).
This is a cis-window COVERAGE AUDIT, not cohort-/assay-proven Coloc input.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import urllib.request
from collections import Counter
from pathlib import Path

RECORD="https://zenodo.org/records/4715946"
FILENAME="gene_counts_Ensembl_105_phenotype_metadata.tsv.gz"
URL=f"{RECORD}/files/{FILENAME}?download=1"
EXPECTED_MD5="b999b21dcffeecbfe0bf886a5a91997b"
MAX_SOURCE_BYTES=4_000_000
WINDOW_BP=1_000_000
EXPECTED_GROUPS={"QTD000131":10,"QTD000171":11,"QTD000434":15}
REGION={"chrom":"10","start":102403919,"end":104460593}


def require(condition:bool,reason:str)->None:
    if not condition:raise ValueError(reason)


def fetch_official(open_url=urllib.request.urlopen)->tuple[bytes,dict]:
    req=urllib.request.Request(URL,headers={
         "User-Agent":"IS-Masters-P3A-gene-phenotype-source-audit/1.0"})
    with open_url(req,timeout=90) as r:
        raw=r.read(MAX_SOURCE_BYTES+1)
    require(100_000<len(raw)<=MAX_SOURCE_BYTES,
            "Official phenotype metadata size outside expected bounds")
    actual=hashlib.md5(raw,usedforsecurity=False).hexdigest()
    require(actual==EXPECTED_MD5,"Published phenotype MD5 mismatch: "+actual)
    return raw,{
       "source_record":RECORD,"source_url":URL,
       "release":"eQTL Catalogue Ensembl 105 v22.03.1",
       "publisher_MD5":EXPECTED_MD5,"observed_MD5":actual,
       "observed_SHA256":hashlib.sha256(raw).hexdigest(),
       "compressed_bytes":len(raw),
       "source_reference_gene_positions_only":True,
    }


def load_genes(raw:bytes)->dict:
    reader=csv.DictReader(
        io.StringIO(gzip.decompress(raw).decode("utf-8-sig")),
        delimiter="\t")
    columns=set(reader.fieldnames or [])
    require({"gene_id","chromosome","phenotype_pos"}.issubset(columns),
            "no official phenotype ID/GRCh38 position header")
    rows={}
    for row in reader:
        gene=(row.get("gene_id") or "").split(".")[0]
        if not gene.startswith("ENSG"):
            continue
        chromosome=(row.get("chromosome") or "").removeprefix("chr")
        position_s=row.get("phenotype_pos","")
        try: pos=int(position_s)
        except (TypeError,ValueError):continue
        if chromosome!="10":continue
        require(0<pos<140_000_000,"invalid reference chr10 gene phenotype")
        record={"gene_id":gene,"chromosome":"10","phenotype_pos_GRCh38":pos,
           "gene_name":row.get("gene_name"),
           "gene_strand":row.get("strand"),
           "original_phenotype_id":row.get("phenotype_id")}
        previous=rows.get(gene)
        require(previous is None or previous["phenotype_pos_GRCh38"]==pos,
                "ambiguous phenotype_pos from official reference: "+gene)
        rows[gene]=record
    require(len(rows)>300,"unexpectedly empty official gene annotation")
    return rows


def audit(manifest:dict,qtl:dict,positions:dict,source_info:dict)->dict:
    groups=manifest.get("groups",[])
    require(manifest.get("schema")==
                "is-p3a-36-real-group-keys-from-legacy-393-v1"
            and manifest.get("scientific_status")=="HOLD"
            and len(groups)==36
            and len({(g["region_id"],g["QTL_dataset"],g["gene_id"])
                      for g in groups})==36
            and len({g["gene_id"] for g in groups})==20
            and Counter(g["QTL_dataset"] for g in groups)==EXPECTED_GROUPS
            and all(g["region_id"]=="IS_XDATA_G0015"
                    and not g["approved_for_new_coloc"] for g in groups),
            "source 36-gene manifest changed")
    require(qtl.get("schema")==
                "is-p3a-36-original-eqtl-catalogue-indexed-chr10-audit-v1"
            and qtl.get("technical_status")==
                "PASS_3_OFFICIAL_INDEXED_QTL_CONTEXT_REPORTS"
            and qtl.get("scientific_status")=="HOLD"
            and len(qtl.get("groups",[]))==36
            and qtl.get("converted_GRCh38_queried_context")==
                {**REGION,"is_lower_bound_of_full_cis_with_guaranteed_coverage":False}
            and qtl.get("gene_full_cis_window_coverage_verified") is False
            and qtl.get("coloc_abf_eligible")==0,
            "QTL context source proof incompatible with this build")
    qr={(g["region_id"],g["QTL_dataset"],g["gene_id"]):g
        for g in qtl["groups"]}
    require(set(qr)=={(g["region_id"],g["QTL_dataset"],g["gene_id"])
                  for g in groups},"assay keys changed between source runs")
    unique=sorted({g["gene_id"] for g in groups})
    group_rows=[]
    rows={}
    for gene in unique:
        x=positions.get(gene)
        if x is None:
            rows[gene]={"gene_id":gene,"official_gene_reference_found":False,
                         "status":"HOLD_MISSING_IN_ENSEMBL_105_REFERENCE"}
            continue
        pos=x["phenotype_pos_GRCh38"]
        start=max(1,pos-WINDOW_BP)
        end=pos+WINDOW_BP
        seen=max(0,min(REGION["end"],end)-max(REGION["start"],start)+1)
        total=end-start+1
        rows[gene]={**x,"official_gene_reference_found":True,
            "predicted_cis_1Mb_left":start,"predicted_cis_1Mb_right":end,
            "predicted_cis_window_basepairs":total,
            "intersected_context_basepairs":seen,
            "reference_cis_overlap_fraction":round(seen/total,7),
            "missing_left_bp":max(0,REGION["start"]-start),
            "missing_right_bp":max(0,end-REGION["end"]),
            "reference_expected_window_fully_inside_common_context":
                start>=REGION["start"] and end<=REGION["end"],
            "reference_window_is_study_specific_QTL_nominal_cis_proof":False,
            "status":"HOLD_STUDY_SPECIFIC_CIS_VARIANT_COVERAGE",
        }
    for a in groups:
        key=(a["region_id"],a["QTL_dataset"],a["gene_id"])
        t=rows[a["gene_id"]]
        old=qr[key]
        require((old.get("GRCh38_context_gene_associations_seen") or 0)>0,
                "real source QTL associations missing in group: "+str(key))
        group_rows.append({
            "region_id":key[0],"QTL_dataset":key[1],"gene_id":key[2],
            "gene_symbol":a["gene_symbol"],
            "observed_context_candidate_rows":
                old["GRCh38_context_gene_associations_seen"],
            "observed_context_QC_nonpal":
                old["GRCh38_context_valid_nonpal_unique"],
            "reference_annotation_gene_105_pos_GRCh38":
                t.get("phenotype_pos_GRCh38"),
            "expected_plus_minus_1Mb_cis":
                [t["predicted_cis_1Mb_left"],t["predicted_cis_1Mb_right"]]
                    if t.get("official_gene_reference_found") else None,
            "observed_context_fraction_of_reference_cis_span":
                t.get("reference_cis_overlap_fraction"),
            "reference_cis_fully_within_previous_context":
                t.get("reference_expected_window_fully_inside_common_context",False),
            "source_study_exact_phenotype_metadata_matched":False,
            "per_SNP_nominal_full_cis_verified":False,
            "GWAS_per_SNP_harmonization_verified":False,
            "coloc_abf_eligible":False,
            "status":"HOLD_STUDY_CIS_AND_GWAS_N_HARMONIZATION",
        })
    require(len(group_rows)==36,"36 source groups missing")
    n_found=sum(t["official_gene_reference_found"] for t in rows.values())
    n_full=sum(t.get("reference_expected_window_fully_inside_common_context",False)
               for t in rows.values())
    return {
        "schema":"is-p3a-official-105-phenotype-cis-window-context-audit-v1",
        "technical_status":("PASS_REF_POS_20_GENES_AND_36_KEYS" if n_found==20
                             else "HOLD_MISSING_GENE_REFERENCE"),
        "scientific_status":"HOLD",
        "official_Ensembl_105_phenotype_reference":source_info,
        "cis_window_inferred_from_qtlmap_publication_bp_each_side":WINDOW_BP,
        "full_gene_specific_nominal_tested_cis_proven":False,
        "reference_105_positions_cannot_approve_dataset_specific_TSS":True,
        "context_GRCh38":REGION,
        "n_unique_genes":20,"n_unique_ref_position_found":n_found,
        "n_reference_cis_entirely_in_prior_context":n_full,
        "n_reference_cis_not_fully_covered":n_found-n_full,
        "n_real_groups":36,
        "n_group_coloc_abf_eligible":0,"coloc_abf_executed":0,
        "n_group_coloc_susie_eligible":0,"coloc_susie_executed":0,
        "genes":rows,"groups":group_rows,
    }


def main()->int:
    parser=argparse.ArgumentParser()
    parser.add_argument("--manifest",type=Path,required=True)
    parser.add_argument("--qtl-context",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    try:
        m=json.loads(args.manifest.read_text(encoding="utf8"))
        q=json.loads(args.qtl_context.read_text(encoding="utf8"))
        content,info=fetch_official()
        result=audit(m,q,load_genes(content),info)
    except (ValueError,KeyError,TypeError,OSError,UnicodeError) as e:
        result={"schema":"is-p3a-official-105-phenotype-cis-window-context-audit-v1",
            "technical_status":"HOLD_SOURCE_METADATA",
            "scientific_status":"HOLD","n_group_coloc_abf_eligible":0,
            "coloc_abf_executed":0,"error":str(e)[:300]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",
                        encoding="utf8")
    print(json.dumps({k:result.get(k) for k in
           ("technical_status","n_unique_ref_position_found",
            "n_reference_cis_entirely_in_prior_context",
            "n_reference_cis_not_fully_covered",
            "n_real_groups","n_group_coloc_abf_eligible","error")},
           ensure_ascii=False))
    return 0 if result["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
