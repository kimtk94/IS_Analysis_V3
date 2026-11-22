#!/usr/bin/env python3
"""Primary eQTL Catalogue full nominal vs independently published permutation n_variants.

For 36 authentic IS-P3A QTL gene×tissue assays from 3 studies, compare *unique*
variants in full original Tabix-indexed cis nominal associations against the
QTLtools permutation-pass n_variants (per-gene tested cis count). Validate
molecular_trait_id vs the pinned Ensembl105 phenotype metadata and source
genotype allele count summaries. Never infer regression model effective N or
declare formal Coloc eligibility from matching counts alone.
"""
from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import importlib
import io
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

BASE=importlib.import_module("scripts.is.audit_p3a_indexed_qtl_real_chr10")
UNION=importlib.import_module("scripts.is.audit_p3a_official_cis_union_36")
PERMUTED_HEADER={"molecular_trait_id","n_variants"}
MAX_PERM_ZIP=30_000_000
MAX_PERM_PLAIN=120_000_000
MAX_NOMINAL_ROWS=1_200_000
QTL_N={"QTD000131":387,"QTD000171":205,"QTD000434":560}


def require(value,msg):
    if not value:raise ValueError(msg)


def pin_inputs(manifest,union,reference,official,allele_n):
    orig=manifest.get("groups",[])
    require(manifest.get("schema")=="is-p3a-36-real-group-keys-from-legacy-393-v1"
            and manifest.get("scientific_status")=="HOLD" and len(orig)==36
            and len({(g["region_id"],g["QTL_dataset"],g["gene_id"])
                     for g in orig})==36 and len({g["gene_id"] for g in orig})==20,
            "original 36 assay identities drift")
    require(union.get("schema")=="is-p3a-36-expanded-original-cis-reference-union-v1"
            and union.get("technical_status")==
                "PASS_THREE_REAL_INDEXED_REFERENCE_UNION_SURVEYS"
            and union.get("scientific_status")=="HOLD"
            and union.get("n_genuine_assay_groups")==36
            and union.get("n_approved_coloc_abf")==0
            and union.get("cis_union_reference_Ensembl105_GRCh38_1based")==
                UNION.EXPECTED_WINDOW,
            "original source union missing or scientifically misapproved")
    require(reference.get("schema")==
                "is-p3a-official-105-phenotype-cis-window-context-audit-v1"
            and reference.get("technical_status")=="PASS_REF_POS_20_GENES_AND_36_KEYS"
            and len(reference.get("genes",{}))==20
            and reference.get("scientific_status")=="HOLD",
            "Ensembl105 official 20 gene metadata missing")
    require(official.get("technical_status")==
            "PASS_3_SOURCE_IDENTITIES_IN_IMMUTABLE_CATALOGUE"
            and set(official["QTL_datasets"])==set(QTL_N),
            "original source release metadata missing")
    require(allele_n.get("schema")=="is-p3a-primary-gwas-yaml-QTL-allele-N-36-v1"
            and allele_n.get("technical_status")==
                "PASS_PRIMARY_PUBLISHER_METADATA_AND_36_QTL_AN_SOURCE_READ"
            and len(allele_n.get("groups",[]))==36
            and allele_n.get("formal_coloc_abf_eligible")==0,
            "real per-SNP QTL an/2 count baseline not validated")
    groups={(g["QTL_dataset"],g["gene_id"]) for g in orig}
    for tag,collection,columns in (
        ("union",union["groups"],("QTL_dataset","gene_id")),
        ("per-SNP an",allele_n["groups"],("QTL_dataset","gene_id")),
        ("Ensembl105",reference["groups"],("QTL_dataset","gene_id"))):
        keys={(x[columns[0]],x[columns[1]]) for x in collection}
        require(len(keys)==36 and keys==groups,tag+" candidate group mismatch")
    require(all(official["QTL_datasets"][q]["nominal_catalogue_sample_size"]==n
                for q,n in QTL_N.items()),
            "publisher tissue/sample metadata changed")
    return {(g["QTL_dataset"],g["gene_id"]):g for g in orig}


def permutation_url(q,official):
    stem=official["QTL_datasets"][q]["source_https_probe"]
    require(stem.endswith("/"+q+".all.tsv.gz") and
            stem.startswith("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"),
            "permutation path not derived from pinned official original")
    return stem.removesuffix(".all.tsv.gz")+".permuted.tsv.gz"


def fetch_permuted(q,official,fetch=urllib.request.urlopen):
    url=permutation_url(q,official)
    req=urllib.request.Request(url,headers={
        "User-Agent":"IS-Masters-P3A-primary-permutation-cis-tested-count-audit/1.0"})
    with fetch(req,timeout=95) as r:
        require(getattr(r,"status",200)==200,"publisher permutation HTTP status not 200")
        parsed=urllib.parse.urlsplit(r.geturl())
        require(parsed.scheme=="https" and parsed.hostname=="ftp.ebi.ac.uk",
                "permutation source unexpected redirect")
        raw=r.read(MAX_PERM_ZIP+1)
    require(200<len(raw)<=MAX_PERM_ZIP,
            "original permutation gzip empty or outside bounded size")
    try:
        plain=gzip.decompress(raw)
    except (EOFError,OSError) as exc:
        raise ValueError("original permutation gzip invalid") from exc
    require(100<len(plain)<=MAX_PERM_PLAIN,
            "original permutation expanded bytes outside cap")
    rows=csv.DictReader(io.StringIO(plain.decode("utf-8-sig")),delimiter="\t")
    require(PERMUTED_HEADER.issubset(set(rows.fieldnames or [])),
            "original permutation n_variants or phenotype column absent")
    found=collections.defaultdict(list)
    total=0
    for row in rows:
        total+=1
        require(total<200_000 and None not in row,
                "publisher permutation source row bound/physical layout")
        trait=(row.get("molecular_trait_id") or "").split(".")[0]
        count=row.get("n_variants","")
        if not trait.startswith("ENSG"):continue
        try:
            n=int(count)
            require(n>0,"permutation variant count not positive")
        except (TypeError,ValueError):
            raise ValueError("publisher permutation n_variants invalid")
        found[trait].append({"n_variants":n,
            "molecular_trait_id_original":row.get("molecular_trait_id"),
            "molecular_trait_object_id":row.get("molecular_trait_object_id"),
            "gene_id":row.get("gene_id"),
            "p_beta":row.get("p_beta")})
    return found,{"source_url":url,"gzip_bytes":len(raw),
                  "gzip_SHA256":hashlib.sha256(raw).hexdigest(),
                  "number_original_permutation_rows":total,
                  "original_columns":list(rows.fieldnames or []),
                  "original_phenotypes_with_permutations":len(found)}


def count_original_nominal(q,official,genes,reference,open_tabix=None):
    if open_tabix is None:
        import pysam
        open_tabix=pysam.TabixFile
    reg=UNION.EXPECTED_WINDOW
    src=official["QTL_datasets"][q]["source_https_probe"]
    groups={g:{"gene":g,
               "nominal_gene_rows_physical":0,
               "all_unique_variant_keys":set(),
               "observed_molecular_trait_ids":set(),
               "observed_outside_reference_1Mb":0,
               "source_variant_class_counts":collections.Counter(),
               "duplicate_rsid_source_rows":0}
            for g in genes}
    digest=hashlib.sha256(("\t".join(BASE.OFFICIAL_COLUMN_ORDER)+"\n").encode())
    n=0
    with open_tabix(src) as tabix:
        chrom="10" if "10" in tabix.contigs else "chr10" if "chr10" in tabix.contigs else None
        require(chrom is not None,"original nominal Tabix missing chr10")
        for line in tabix.fetch(chrom,reg["start"]-1,reg["end"]):
            if isinstance(line,bytes):line=line.decode("utf-8")
            line=line.rstrip("\r\n")+"\n"
            digest.update(line.encode("utf-8"))
            n+=1
            require(n<=MAX_NOMINAL_ROWS,
                    "original nominal bounded region exceeded review cap")
            cells=line.rstrip("\n").split("\t")
            require(len(cells)==len(BASE.OFFICIAL_COLUMN_ORDER),
                    "original 19 physical nominal columns changed")
            x=dict(zip(BASE.OFFICIAL_COLUMN_ORDER,cells))
            gene=x["gene_id"]
            if gene not in genes:continue
            p=int(x["position"])
            require(x["chromosome"] in ("10","chr10")
                    and reg["start"]<=p<=reg["end"],
                    "source nominal variant outside pinned union")
            g=groups[gene]
            g["nominal_gene_rows_physical"]+=1
            trait=x["molecular_trait_id"].split(".")[0]
            g["observed_molecular_trait_ids"].add(trait)
            if not (reference[gene]["predicted_cis_1Mb_left"]<=p<=
                    reference[gene]["predicted_cis_1Mb_right"]):
                g["observed_outside_reference_1Mb"]+=1
            require(x["variant"]==f"chr10_{p}_{x['ref']}_{x['alt']}",
                    "unexpected nominal canonical variant ID")
            key=(p,x["ref"],x["alt"],x["molecular_trait_id"])
            if key in g["all_unique_variant_keys"]:
                g["duplicate_rsid_source_rows"]+=1
                continue
            g["all_unique_variant_keys"].add(key)
            kind=x["type"]
            g["source_variant_class_counts"][kind]+=1
    for gene,v in groups.items():
        v["nominal_all_variant_count_unique"]=len(v.pop("all_unique_variant_keys"))
        v["observed_molecular_trait_ids"]=sorted(v["observed_molecular_trait_ids"])
        v["source_variant_class_counts"]=dict(v["source_variant_class_counts"])
    return groups,{"query_source":src,"indexed_region_rows":n,
                   "tabix_region_source_stream_SHA256":digest.hexdigest()}


def evaluate(manifest,union,reference,official,allele_n,perm_rows,nominal_rows,provenance):
    by_original=pin_inputs(manifest,union,reference,official,allele_n)
    previous={(x["QTL_dataset"],x["gene_id"]):x for x in union["groups"]}
    observed_n={(x["QTL_dataset"],x["gene_id"]):
                x["original_variant_allele_N_summary"] for x in allele_n["groups"]}
    rows=[]
    for key,g in sorted(by_original.items()):
        q,gene=key
        nom=nominal_rows[q][gene]
        perms=perm_rows[q].get(gene,[])
        require(nom["nominal_gene_rows_physical"]==
                previous[key]["extended_context_gene_association_rows"],
                "source full nominal physical rows disagree with prior main run")
        require(nom["observed_outside_reference_1Mb"]==0,
                "an original source QTL row lies outside candidate gene reference cis")
        expected=observed_n[key]["count_valid_nonpal_original_QTL"]
        # Nonpal SNP count is a strict subset of all-source nominal variants.
        require(nom["nominal_all_variant_count_unique"]>=expected,
                "permutation comparison source nominal count below prior valid SNP QC")
        ids=nom["observed_molecular_trait_ids"]
        all_genes_equivalent=(ids==[gene])
        perm_unique=len(perms)==1
        source_nvariants=perms[0]["n_variants"] if perm_unique else None
        matching=(perm_unique and all_genes_equivalent and
                  source_nvariants==nom["nominal_all_variant_count_unique"])
        rows.append({
            "region_id":g["region_id"],"QTL_dataset":q,
            "gene_id":gene,"gene_symbol":g["gene_symbol"],
            "reference_gene_Ensembl105_GRCh38_phenotype_pos":
                reference["genes"][gene]["phenotype_pos_GRCh38"],
            "nominal_gene_rows_physical":nom["nominal_gene_rows_physical"],
            "nominal_distinct_all_variant_records":
                nom["nominal_all_variant_count_unique"],
            "nominal_strict_nonpal_SNP_source_records":expected,
            "original_variant_class_counts":nom["source_variant_class_counts"],
            "original_duplicate_rsid_rows":nom["duplicate_rsid_source_rows"],
            "nominal_molecular_trait_ids":ids,
            "molecular_trait_id_matches_Ensembl_gene_id":all_genes_equivalent,
            "original_permutation_rows_for_trait":len(perms),
            "source_original_cis_n_variants_from_permutation":source_nvariants,
            "published_permutation_vs_recovered_nominal_unique_agree":bool(matching),
            "full_tested_cis_variant_COUNT_independently_reconciled":bool(matching),
            "original_per_phenotype_tested_variant_membership_independently_verified":
                False,
            "original_dataset_specific_effective_regression_N_approved":False,
            "primary_GWAS_per_SNP_effective_N_approved":False,
            "formal_coloc_abf_eligible":False,
            "formal_coloc_susie_eligible":False,
            "status":("PASS_SOURCE_CIS_TESTED_NVAR_COUNT_STILL_MODEL_N_HOLD"
                      if matching else "HOLD_CIS_NVAR_OR_PHENOTYPE_MISMATCH"),
        })
    require(len(rows)==36 and len({x["gene_id"] for x in rows})==20,
            "original 36 gene tissue keys missing in permutation comparison")
    total=sum(x["published_permutation_vs_recovered_nominal_unique_agree"] for x in rows)
    return {
        "schema":"is-p3a-36-original-eQTL-permuted-n_variants-nominal-completeness-v1",
        "technical_status":"PASS_3_ORIGINAL_PERMUTATION_AND_NOMINAL_REGION_SOURCES_READ",
        "scientific_status":"HOLD",
        "n_original_gene_tissue_groups":36,
        "n_unique_genes":20,"n_stroke_loci":1,
        "n_gene_tissue_groups_matching_full_cis_tested_variant_COUNT":total,
        "n_groups_without_matching_source_n_variants":36-total,
        "all_gene_tissue_model_effective_N_approved":False,
        "n_formal_coloc_abf_eligible":0,"n_formal_coloc_abf_executed":0,
        "n_formal_coloc_susie_eligible":0,"n_formal_coloc_susie_executed":0,
        "source_provenance":provenance,"groups":rows,
    }


def main():
    p=argparse.ArgumentParser()
    for key in ("manifest","union","reference","official","allele-n","out"):
        p.add_argument("--"+key,type=Path,required=True)
    args=p.parse_args()
    out={"schema":"is-p3a-36-original-eQTL-permuted-n_variants-nominal-completeness-v1",
         "technical_status":"HOLD_UNVERIFIED_SOURCE",
         "scientific_status":"HOLD","n_formal_coloc_abf_eligible":0,
         "n_formal_coloc_abf_executed":0,
         "n_formal_coloc_susie_eligible":0,"n_formal_coloc_susie_executed":0}
    try:
        def load(key):
            return json.loads(getattr(args,key.replace("-","_")).read_text(encoding="utf-8"))
        manifest,union,reference,official,allele_n=[
            load(x) for x in ("manifest","union","reference","official","allele-n")]
        groups=pin_inputs(manifest,union,reference,official,allele_n)
        perms={};nominals={};provenance={}
        for i,q in enumerate(QTL_N):
            if i:time.sleep(9)
            raw,pr=fetch_permuted(q,official)
            target={gene for dataset,gene in groups if dataset==q}
            nominal,source=count_original_nominal(q,official,target,
                                                   reference["genes"])
            previous=union["QTL_datasets"][q]
            require(source["indexed_region_rows"]==
                    previous["interval_rows"]["region_rows_scanned"]
                    and source["tabix_region_source_stream_SHA256"]==
                    previous["region_stream_SHA256"],
                    "complete original indexed nominal source stream bytes drift")
            perms[q]=raw;nominals[q]=nominal
            provenance[q]={"permuted":pr,"nominal":source}
        out=evaluate(manifest,union,reference,official,allele_n,
                     perms,nominals,provenance)
    except (ValueError,OSError,KeyError,TypeError,UnicodeError,TimeoutError,
            ImportError,OverflowError) as e:
        out["technical_status"]="HOLD_PRIMARY_PERMUTATION_SOURCE_OR_CIS_REFERENCE"
        out["error_type"]=type(e).__name__
        out["error"]=str(e)[:400]
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:out.get(k) for k in (
       "technical_status","scientific_status","n_original_gene_tissue_groups",
       "n_gene_tissue_groups_matching_full_cis_tested_variant_COUNT",
       "n_groups_without_matching_source_n_variants",
       "n_formal_coloc_abf_eligible","error")},ensure_ascii=False))
    return 0 if out["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
