#!/usr/bin/env python3
"""Re-check all 36 recovered P3A follow-up group identifiers from PINNED old Git.

Never replace the former server's 25/11 *source-stat* evidence with 25/11
independent eGene BH overlap. Read actual 2026-10-09 branch TSV files in CI.
Report all 36 authentic groups, 20 unique genes, and precise source blockers.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import re
from collections import Counter
from pathlib import Path

NAMES={
 "readiness":"IS_EAS5_GENE_TISSUE_ALLELE_HARMONIZED_READINESS.tsv",
 "egene":"IS_EAS5_98_BH05_FDR_EGENE_FULLCIS_GWAS_ALLELE_READINESS.tsv",
 "samesnp":"IS_EAS5_SAME_SNP_ASSOCIATION_GENE_TISSUE_QC.tsv",
}
BLOB={
 "readiness":"a1b8d13c9460f31788cde11b1c6e3f65340af749",
 "egene":"a2a77a93df04496cf88a2a0d95fd23d8c79d88d1",
 "samesnp":"ed9a5ddb5df956f861038eef36e1189205fbfeff",
}
KEY=("region_id","QTL_dataset","gene_id")
EXPECTED={"QTD000131":10,"QTD000171":11,"QTD000434":15}

def require(p,reason):
    if not p: raise ValueError(reason)

def key(x):
    return tuple(x.get(k) for k in KEY)

def get_source(root,what):
    path=root/NAMES[what]
    raw=path.read_bytes()
    require(len(raw)<400_000 and len(raw)>200, "source file size invalid: "+what)
    sha=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    require(sha==BLOB[what],"pinned original Git blob SHA differs: "+what)
    r=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")),delimiter="\t")
    rows=list(r)
    require(all(None not in x for x in rows),"broken TSV schema: "+what)
    require(len({key(x) for x in rows})==len(rows),"duplicate source assay key: "+what)
    return rows

def validate(manifest, readiness, egene, samesnp):
    require(manifest.get("schema")=="is-p3a-36-real-group-keys-from-legacy-393-v1"
            and manifest.get("scientific_status")=="HOLD"
            and manifest.get("formal_new_P3A_coloc_executed")==0,
            "source manifest scientific state promoted")
    origin=manifest.get("origin",{})
    require(origin.get("commit_sha")=="0bf646fd6f29322d17da878253f37def24c01c74"
            and [x.get("git_blob_sha") for x in origin.get("files",[])]==list(BLOB.values())
            and origin.get("branch")=="research/is-broad-discovery-20261009",
            "source Git commit and blob provenance mismatch")
    require(len(readiness)==393 and len(egene)==98 and len(samesnp)==393,
            "readiness / BH05 / same SNP input sizes drift")
    by_key={key(x):x for x in readiness}
    by_egene={key(x):x for x in egene}
    by_samesnp={key(x):x for x in samesnp}
    require(len(by_key)==393 and len(by_egene)==98 and len(by_samesnp)==393,
            "duplicate ID in original source")
    follow=[x for x in readiness if x["possible_single_signal_ABF_followup"]=="True"]
    require(len(follow)==36
            and sum(x["possible_single_signal_ABF_followup"]=="False" for x in readiness)==357
            and len({x["gene_id"] for x in follow})==20
            and all(x["region_id"]=="IS_XDATA_G0015" for x in follow),
            "original P3A flag universe drift")
    require(Counter(x["QTL_dataset"] for x in follow)==EXPECTED,
            "QTL source distribution drift")
    require(len([x for x in follow if key(x) in by_egene])==25,
            "public FDR BH05 overlap drift")
    require(all(key(x) in by_samesnp for x in follow),
            "missing original same-SNP association checks")

    saved=manifest.get("groups",[])
    require(isinstance(saved,list) and len(saved)==36
            and len({key(x) for x in saved})==36,"snapshot 36-group keys missing/duplicated")
    require({key(x) for x in saved}=={key(x) for x in follow},
            "snapshot 36 group identifiers not equal to original archived source")
    require(manifest.get("denominators",{}).get("assays_with_public_studywide_BH05_eGene")==25
            and manifest.get("denominators",{}).get("per_group_original_server_source_stats_positive_verified")==0
            and manifest.get("denominators",{}).get("unique_gene_ids_in_followup")==20,
            "25 public FDR group overlaps incorrectly asserted server source-stat proof")
    require(manifest.get("independent_allele_and_QTL_N_approved") is False
            and manifest.get("matched_LD_approved") is False,
            "false biological scientific source approval")
    results=[]
    for row in saved:
        k=key(row);old=by_key[k];r=by_egene.get(k);s=by_samesnp[k]
        require(all(row.get(k)==old[k] for k in KEY)
                and row["gene_symbol"]==old["gene_symbol"]
                and row["archived_ABF_assays"]==int(old["original_archived_ABF_assays"])
                and row["archived_nonpalindromic_BBJ_EAS_SNP_overlap"]==
                    int(old["exact_NONPALINDROMIC_SNP_in_BOTH_BBJ_AIS"]),
                "source group identity or SNP overlap has changed")
        require(row.get("official_studywide_eGene_BH05_assay") is (r is not None),
                "BH05 eGene assay membership mismatch")
        if r:
            require(float(r["official_studywide_BH_q"])<=0.05
                    and abs(row["official_QTL_studywide_BH_q"]-
                        float(r["official_studywide_BH_q"]))<1e-13
                    and row["official_full_cis_tested_variant_count"]==
                        int(r["official_full_cis_tested_n_variants"])
                    and row["fullcis_BBJ_exact_SNPs"]==int(r["NEW_fullcis_BBJ_exact_n"])
                    and row["fullcis_EAS_exact_SNPs"]==int(r["NEW_fullcis_EAS_AIS_exact_n"])
                    and r["any_eGene_to_AIS_coloc_computed"]=="False"
                    and r["source_QTL_genotype_matched_LD_verified"]=="False"
                    and r["BBJ_and_GIGA_independent_stroke_replication"]=="False",
                    "BH05 QTL or original sample LD source drift")
        else:
            require(all(row.get(k) is None for k in (
                "official_QTL_studywide_BH_q","official_full_cis_tested_variant_count",
                "fullcis_BBJ_exact_SNPs","fullcis_EAS_exact_SNPs")),
                "missing eGene was incorrectly assigned summary values")
        require(int(s["same_SNP_QTL_le_1e4_and_AIS_lt_5e8"])==
                row["same_SNP_QTL_1e4_AIS_5e8"]
                and s["QTL_study_ancestry_independently_verified"]=="False"
                and s["genotype_LD_matched_across_QTL_and_GWAS"]=="False"
                and s["colocalization_or_causal_gene_proven"]=="False"
                and row["approved_for_new_coloc"] is False
                and row["archived_GWAS_QTL_coloc_executed"] is False
                and row["original_P3A_source_stat_positive_proven"] is False,
                "scientific stop/ancestry/allele source status falsely approved")
        results.append({
            "region_id":k[0],"QTL_dataset":k[1],"gene_id":k[2],
            "gene_symbol":row["gene_symbol"],
            "source_BH05_eGene_assay":"YES" if r else "NO_SOURCE_BH05_MATCH",
            "official_BH_q":row["official_QTL_studywide_BH_q"],
            "nonpal_SNP_original_BBJ_EAS_region":row["archived_nonpalindromic_BBJ_EAS_SNP_overlap"],
            "complete_cis_tested_SNP_from_BH05_source":row["official_full_cis_tested_variant_count"],
            "QTL_matched_LD_verified":False,
            "GWAS_effect_SE_casefraction_N_verified":False,
            "new_coloc_executed":False,
            "source_stat_group_membership":"UNVERIFIED_25_OF_36_AGGREGATE_ONLY",
            "science_status":"HOLD_SOURCE_EFFECT_N_AND_MATCHED_LD"
        })
    return {
        "schema":"is-p3a-restored-36-pinned-legacy-source-audit-v1",
        "technical_status":"PASS_36_REAL_GROUP_IDS_CROSSCHECKED_TO_THREE_PINNED_ARCHIVED_TABLES",
        "scientific_status":"HOLD",
        "pinned_branch_commit":origin["commit_sha"],
        "restored_group_key_count":36,
        "distinct_gene_id_count":20,
        "by_QTL_dataset":EXPECTED,
        "public_BH05_eGene_assay_count":25,
        "not_in_public_BH05_assay_count":11,
        "original_server_per_group_source_stats_positive_membership_restored":False,
        "source_stats_positive_group_aggregate_only":25,
        "source_stats_absent_group_aggregate_only":11,
        "fully_approved_new_P3A_coloc_groups":0,
        "new_P3A_coloc_computations":0,
        "row_blocker_count":36,
        "source_provenance_is_not_independent_replication":True,
        "groups":results
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-root",type=Path,required=True)
    ap.add_argument("--manifest",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    try:
        original={name:get_source(a.source_root,name) for name in NAMES}
        result=validate(json.loads(a.manifest.read_text(encoding="utf-8")),
                        original["readiness"],original["egene"],original["samesnp"])
    except (OSError,ValueError,TypeError,KeyError) as e:
        result={"schema":"is-p3a-restored-36-pinned-legacy-source-audit-v1",
                "technical_status":"HOLD_ORIGINAL_GROUP_SOURCE_MISSING_OR_CHANGED",
                "scientific_status":"HOLD","new_P3A_coloc_computations":0,
                "error":str(e)[:260]}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k in
        ("technical_status","restored_group_key_count","distinct_gene_id_count",
         "public_BH05_eGene_assay_count","scientific_status","error")},ensure_ascii=False))
    return 0 if result["technical_status"].startswith("PASS_") else 1

if __name__=="__main__":
    raise SystemExit(main())
