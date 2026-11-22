#!/usr/bin/env python3
"""IS P3A: one genuine indexed QTL cis-window union survey per public dataset.

Official Ensembl105 phenotype positions yield a *reference union* of 20 gene
cis windows. This does not independently prove the original study's tested
phenotype centers, all-variant genome liftover, GWAS harmonization, sample N,
or any coloc.abf eligibility.

Only public summary QC counts are uploaded. No raw SNP/effect rows or patients.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import tempfile
import time
from collections import Counter
from pathlib import Path

BASE=importlib.import_module("scripts.is.audit_p3a_indexed_qtl_real_chr10")
CIS=importlib.import_module("scripts.is.audit_p3a_ensembl105_cis_reference")
UCSC=importlib.import_module("scripts.is.audit_p3a_chr10_ucsc_liftover_fasta")
EXPECTED_GROUPS={"QTD000131":10,"QTD000171":11,"QTD000434":15}
EXPECTED_WINDOW={"chrom":"10","start":101837472,"end":105353833}
OLD_WINDOW={"chrom":"10","start":102403919,"end":104460593}
SOURCE_FIELDS=BASE.OFFICIAL_COLUMN_ORDER
MAX_READS_PER_DATASET=BASE.MAX_LINES


def require(ok:bool,reason:str)->None:
    if not ok:
        raise ValueError(reason)


def check_inputs(manifest:dict,baseline:dict,reference:dict,official:dict,liftover:dict):
    old=BASE.verify_inputs(manifest,official,liftover)
    require({k:old[k] for k in ("chrom","start","end")}==OLD_WINDOW,
            "original verified 2Mb context changed")
    groups=manifest.get("groups",[])
    require(len(groups)==36
            and len({(g["region_id"],g["QTL_dataset"],g["gene_id"])
                     for g in groups})==36
            and Counter(g["QTL_dataset"] for g in groups)==EXPECTED_GROUPS
            and len({g["gene_id"] for g in groups})==20,
            "authentic 36 QTL group identities missing")
    require(reference.get("schema")==
                "is-p3a-official-105-phenotype-cis-window-context-audit-v1"
            and reference.get("technical_status")=="PASS_REF_POS_20_GENES_AND_36_KEYS"
            and reference.get("scientific_status")=="HOLD"
            and reference.get("n_unique_ref_position_found")==20
            and reference.get("n_real_groups")==36
            and reference.get("n_group_coloc_abf_eligible")==0
            and reference.get("n_reference_cis_entirely_in_prior_context")==1
            and reference.get("full_gene_specific_nominal_tested_cis_proven") is False
            and reference.get("reference_105_positions_cannot_approve_dataset_specific_TSS") is True
            and reference.get("official_Ensembl_105_phenotype_reference",{}).get("publisher_MD5")==
                CIS.EXPECTED_MD5,
            "Ensembl 105 official reference MD5 or source-science gate incomplete")
    require(baseline.get("schema")==
                "is-p3a-36-original-eqtl-catalogue-indexed-chr10-audit-v1"
            and baseline.get("technical_status")==
                "PASS_3_OFFICIAL_INDEXED_QTL_CONTEXT_REPORTS"
            and baseline.get("scientific_status")=="HOLD"
            and baseline.get("coloc_abf_eligible")==0
            and baseline.get("gene_full_cis_window_coverage_verified") is False
            and baseline.get("source_datasets_with_actual_indexed_rows_read")==3
            and len(baseline.get("groups",[]))==36,
            "previous real 2Mb original-source QTL baseline invalid")
    originals={(g["region_id"],g["QTL_dataset"],g["gene_id"])
               for g in groups}
    for name,source in (("reference",reference),("baseline",baseline)):
        rows=source.get("groups",[])
        keys={(x["region_id"],x["QTL_dataset"],x["gene_id"]) for x in rows}
        require(len(keys)==36 and keys==originals,
                name+" candidate assay identity mismatch")
    return bounds(reference["genes"]),groups


def bounds(genes:dict)->dict:
    require(len(genes)==20,"not 20 unique Ensembl105 reference genes")
    require(all(g["official_gene_reference_found"] is True
                and g["chromosome"]=="10"
                and isinstance(g["phenotype_pos_GRCh38"],int)
                and g.get("reference_window_is_study_specific_QTL_nominal_cis_proof") is False
                and g["predicted_cis_1Mb_left"]==
                    g["phenotype_pos_GRCh38"]-1_000_000
                and g["predicted_cis_1Mb_right"]==
                    g["phenotype_pos_GRCh38"]+1_000_000
                for g in genes.values()),
            "original reference positions/cis declarations inconsistent")
    reg={"chrom":"10",
         "start":min(g["predicted_cis_1Mb_left"] for g in genes.values()),
         "end":max(g["predicted_cis_1Mb_right"] for g in genes.values())}
    require(reg==EXPECTED_WINDOW,"expected source reference union changed; HOLD")
    require(reg["start"]<OLD_WINDOW["start"] and
            reg["end"]>OLD_WINDOW["end"],
            "union did not contain entire old original context")
    return reg


def update_position_stats(pos:int,gene:str,genes:dict,summary:dict)->None:
    bounds_gene=genes[gene]
    outside=(pos < bounds_gene["predicted_cis_1Mb_left"]
             or pos > bounds_gene["predicted_cis_1Mb_right"])
    item=summary[gene]
    item["source_rows"]+=1
    item["observed_min_GRCh38_pos"]=min(item["observed_min_GRCh38_pos"] or pos,pos)
    item["observed_max_GRCh38_pos"]=max(item["observed_max_GRCh38_pos"] or pos,pos)
    if outside:
        item["observed_rows_outside_Ensembl105_reference_cis"]+=1


def scan_one(url:str,reg:dict,sequence:bytes,genes:dict,tabix_factory=None)->dict:
    if tabix_factory is None:
        import pysam
        tabix_factory=pysam.TabixFile
    per_gene={g:{"source_rows":0,"observed_rows_outside_Ensembl105_reference_cis":0,
                 "observed_min_GRCh38_pos":None,"observed_max_GRCh38_pos":None}
              for g in genes}
    with tabix_factory(url) as indexed:
        contigs=set(indexed.contigs)
        chrom="10" if "10" in contigs else "chr10" if "chr10" in contigs else None
        require(chrom is not None,"original .tbi chromosome 10 missing")
        def source_rows():
            yield "\t".join(SOURCE_FIELDS)+"\n"
            for line in indexed.fetch(chrom,reg["start"]-1,reg["end"]):
                row=(line.decode("utf8") if isinstance(line,bytes) else line)
                cells=row.rstrip("\r\n").split("\t")
                require(len(cells)==len(SOURCE_FIELDS),
                        "original indexed QTL malformed physical 19-column row")
                gene=cells[SOURCE_FIELDS.index("gene_id")]
                if gene in genes:
                    pos=int(cells[SOURCE_FIELDS.index("position")])
                    update_position_stats(pos,gene,genes,per_gene)
                yield row+"\n"
        record=BASE.read_rows(source_rows(),sequence,set(genes),reg)
    require(record["interval_rows"]["region_rows_scanned"]>0,
            "extended actual remote source empty")
    require(sum(p["source_rows"] for p in per_gene.values())==
             record["interval_rows"]["candidate_gene_rows_seen"],
            "original QTL source rows disagree between independent parser tracks")
    record["gene_study_reference_cis_checks"]=per_gene
    record["original_source_tabix_chrom"]=chrom
    record["source_nominal_header_from_pinned_release"]=list(SOURCE_FIELDS)
    return record


def join_results(manifest,baseline,reference,reg,datasets):
    groups=manifest["groups"]
    old={(g["QTL_dataset"],g["gene_id"]):g for g in baseline["groups"]}
    cand=[]
    for g in groups:
        d=g["QTL_dataset"]
        gid=g["gene_id"]
        row=datasets.get(d,{})
        gene=row.get("gene_rows",{}).get(gid,{})
        observed=row.get("gene_study_reference_cis_checks",{}).get(gid,{})
        prior=old[(d,gid)]
        first=prior.get("GRCh38_context_gene_associations_seen")
        current=gene.get("candidate_rows")
        before=prior.get("GRCh38_context_valid_nonpal_unique")
        now=gene.get("valid_nonpal_SNPs")
        if row.get("technical_status")=="PASS_REAL_INDEXED_SOURCE":
            require(all(isinstance(x,int) and x>=0
                       for x in (first,current,before,now)),
                    "QTL source count cannot compare against frozen original source")
            require(current>=first and now>=before and
                    observed.get("source_rows")==current,
                    "extended QTL query is not reproducible superset of old chr10 query")
        cand.append({
            "region_id":g["region_id"],"QTL_dataset":d,
            "gene_id":gid,"gene_symbol":g["gene_symbol"],
            "old_context_gene_association_rows":first,
            "extended_context_gene_association_rows":current,
            "newly_recovered_gene_association_rows":
                current-first if isinstance(current,int) else None,
            "old_context_valid_nonpal_SNP_assay_rows":before,
            "extended_context_valid_nonpal_SNP_assay_rows":now,
            "extra_valid_nonpal_SNP_assay_rows":
                now-before if isinstance(now,int) else None,
            "source_rows_outside_Ensembl105_gene_plus_minus_1Mb":
                observed.get("observed_rows_outside_Ensembl105_reference_cis"),
            "source_observed_min_max_GRCh38_pos":
                [observed["observed_min_GRCh38_pos"],
                 observed["observed_max_GRCh38_pos"]]
                    if observed and observed["source_rows"]>0 else None,
            "reference_study_specific_phenotype_ID_confirmed":False,
            "full_nominal_tested_variants_independently_approved":False,
            "QTL_GWAS_signed_effect_alleles_matched":False,
            "gene_specific_QTL_effective_N_approved":False,
            "GWAS_phenotype_case_control_N_approved":False,
            "formal_coloc_abf_eligible":False,
            "status":"HOLD_PER_ASSAY_STUDY_CIS_AND_GWAS_HARMONIZATION",
        })
    n=sum(x.get("technical_status")=="PASS_REAL_INDEXED_SOURCE"
          for x in datasets.values())
    total=lambda key:sum(x.get(key) or 0 for x in cand)
    out={
        "schema":"is-p3a-36-expanded-original-cis-reference-union-v1",
        "technical_status":("PASS_THREE_REAL_INDEXED_REFERENCE_UNION_SURVEYS"
                            if n==3 else "HOLD_INCOMPLETE_INDEXED_QTL_SOURCES"),
        "scientific_status":"HOLD",
        "cis_union_reference_Ensembl105_GRCh38_1based":reg,
        "reference_cis_range_covers_20_reference_gene_windows":True,
        "study_original_gene_specific_cis_tested_variants_proven":False,
        "complete_QTL_GWAS_per_variant_harmonization_proven":False,
        "single_locus_not_36_independent_stroke_loci":True,
        "n_unique_genes":20,"n_genuine_assay_groups":36,
        "n_approved_coloc_abf":0,"n_executed_coloc_abf":0,
        "n_approved_coloc_susie":0,"n_executed_coloc_susie":0,
        "source_datasets_with_actual_union_query":n,
        "original_common_context_GRCh38_1based":OLD_WINDOW,
        "total_old_candidate_rows":total("old_context_gene_association_rows"),
        "total_extended_candidate_rows":total("extended_context_gene_association_rows"),
        "total_new_candidate_rows":total("newly_recovered_gene_association_rows"),
        "total_old_valid_nonpal_assay_rows":total("old_context_valid_nonpal_SNP_assay_rows"),
        "total_extended_valid_nonpal_assay_rows":total("extended_context_valid_nonpal_SNP_assay_rows"),
        "total_new_valid_nonpal_assay_rows":total("extra_valid_nonpal_SNP_assay_rows"),
        "QTL_datasets":datasets,"groups":cand,
    }
    require(len(cand)==36 and len({x["gene_id"] for x in cand})==20,
            "36 authentic groups not returned")
    return out


def main():
    p=argparse.ArgumentParser()
    for name in ("manifest","baseline","cis-reference","official","liftover","out"):
        p.add_argument("--"+name,type=Path,required=True)
    args=p.parse_args()
    out={"schema":"is-p3a-36-expanded-original-cis-reference-union-v1",
         "technical_status":"HOLD_SOURCE_OR_REFERENCE_UNAVAILABLE",
         "scientific_status":"HOLD","n_approved_coloc_abf":0,
         "n_executed_coloc_abf":0,"n_approved_coloc_susie":0,
         "n_executed_coloc_susie":0,
         "groups":[],"QTL_datasets":{}}
    try:
        read=lambda k:json.loads(getattr(args,k.replace("-","_")).read_text(encoding="utf8"))
        manifest=read("manifest")
        baseline=read("baseline")
        reference=read("cis-reference")
        official=read("official")
        lift=read("liftover")
        reg,groups=check_inputs(manifest,baseline,reference,official,lift)
        datasets={}
        with tempfile.TemporaryDirectory(prefix="is_p3a_36_cis_union") as temp:
            root=Path(temp)
            official_ref=UCSC.fetch_checked("hg38",root/"hg38_chr10.fa.gz")
            require(official_ref["observed_SHA256"]==
                    lift["source_provenance"]["hg38"]["observed_SHA256"],
                    "UCSC complete GRCh38 FASTA SHA256 drift since prior approved CI")
            seq=UCSC.chromosome_fasta_sequence(root/"hg38_chr10.fa.gz")
            require(hashlib.sha256(seq).hexdigest()==
                    lift["chr10_uncompressed_fasta_base_SHA256"]["GRCh38_chr10"],
                    "GRCh38 reference raw chromosome sequence drift")
            for i,(dataset,expected) in enumerate(EXPECTED_GROUPS.items()):
                if i:time.sleep(9)  # Limit repeat HTTP range requests to EBI
                source=official["QTL_datasets"][dataset]
                genes={g["gene_id"]:reference["genes"][g["gene_id"]]
                       for g in groups if g["QTL_dataset"]==dataset}
                url=source["source_https_probe"]
                try:
                    result=scan_one(url,reg,seq,genes)
                    result["technical_status"]="PASS_REAL_INDEXED_SOURCE"
                    result["original_url"]=url
                    result["catalogue_study_n_not_approved_per_gene_effective_N"]=(source["nominal_catalogue_sample_size"])
                    datasets[dataset]=result
                except (ValueError,OSError,OverflowError,TimeoutError,UnicodeError) as e:
                    datasets[dataset]={"technical_status":"HOLD_REMOTE_INDEXED_SOURCE",
                                       "error_type":type(e).__name__,
                                       "error":str(e)[:250],"gene_rows":{},
                                       "gene_study_reference_cis_checks":{}}
                    for later in list(EXPECTED_GROUPS)[i+1:]:
                        datasets[later]={"technical_status":"HOLD_NOT_QUERIED_AFTER_SOURCE_ERROR",
                                         "gene_rows":{},"gene_study_reference_cis_checks":{}}
                    break
        out=join_results(manifest,baseline,reference,reg,datasets)
    except (ValueError,TypeError,KeyError,OSError,ImportError) as e:
        out["technical_status"]="HOLD_REFERENCE_OR_GROUP_VALIDATION_FAILED"
        out["error"]=str(e)[:300]
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:out.get(k) for k in (
         "technical_status","scientific_status","n_genuine_assay_groups",
         "source_datasets_with_actual_union_query","total_old_candidate_rows",
         "total_extended_candidate_rows","total_new_valid_nonpal_assay_rows",
         "n_approved_coloc_abf","error")},ensure_ascii=False))
    return 0 if out["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
