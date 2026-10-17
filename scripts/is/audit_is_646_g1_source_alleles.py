#!/usr/bin/env python3
"""IS original 646 GTEx/GWAS variant source G1 structural harmonization audit.

Exact GRCh38 variant identity/GTEx fields and GWAS-MAF internal arithmetic.
Does not prove GWAS effect allele orientation or GRCh37->38 liftover.
No changes to input data. Reports missing local source, never PASS.
"""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path

def data(path):
    with path.open(encoding="utf8",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def filehash(p):
    sha=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            sha.update(b)
    return sha.hexdigest()

def parse_record(z):
    fields=("variant","chromosome","position","ref","alt","match_key",
            "variant_id_hg38","gwas_eaf","gwas_maf","eqtl_maf","maf",
            "eqtl_beta","beta","eqtl_se","se","harmonization","variant_id")
    if any(k not in z or z[k]=="" for k in fields):
        raise ValueError("SOURCE_ALLELE_REQUIRED_FIELD_MISSING")
    ch=z["chromosome"].removeprefix("chr")
    pos=int(z["position"])
    ref=z["ref"].upper()
    alt=z["alt"].upper()
    if ref==alt or any(x not in "ACGT" for x in ref+alt):
        raise ValueError("QTL_ALLELE_NONSTANDARD_OR_IDENTICAL")
    target=f"{ch}:{pos}:{ref}:{alt}"
    v=f"chr{ch}_{pos}_{ref}_{alt}"
    if z["variant"]!=v or z["match_key"]!=target or z["variant_id_hg38"]!=target:
        raise ValueError("GTEX_GRCH38_ALLELE_IDENTITY_MISMATCH")
    if z["harmonization"]!="EXACT_REF_ALT_GRCH38":
        raise ValueError("SOURCE_HARMONIZATION_TAG_UNEXPECTED")
    if abs(float(z["eqtl_beta"])-float(z["beta"]))>1e-9 or abs(float(z["eqtl_se"])-float(z["se"]))>1e-9:
        raise ValueError("EQTL_BETA_SE_NOT_SOURCE")
    if abs(float(z["eqtl_maf"])-float(z["maf"]))>1e-9:
        raise ValueError("EQTL_MAF_NOT_SOURCE")
    eaf=float(z["gwas_eaf"])
    gmaf=float(z["gwas_maf"])
    qmaf=float(z["eqtl_maf"])
    if any(not math.isfinite(x) for x in (eaf,gmaf,qmaf)) or not(0<eaf<1 and 0<gmaf<=.5 and 0<qmaf<=.5):
        raise ValueError("INVALID_MAF_OR_EAF")
    if abs(gmaf-min(eaf,1-eaf))>1e-9:
        raise ValueError("GWAS_EAF_MAF_INCONSISTENT")
    parts=z["variant_id"].split(":")
    if len(parts)!=4 or not parts[1].isdigit() or parts[2] not in ("A","C","G","T") or parts[3] not in ("A","C","G","T"):
        raise ValueError("GRCH37_GWAS_VARIANT_ID_INVALID")
    return {
        "is_palindromic":int({ref,alt} in ({"A","T"},{"C","G"})),
        "is_palindromic_ambig_eaf":int({ref,alt} in ({"A","T"},{"C","G"}) and .4<=eaf<=.6),
        "is_gwas_gtEx_maf_diff_gt_0_1":int(abs(gmaf-qmaf)>.1),
        "is_gwas_gtEx_maf_diff_gt_0_2":int(abs(gmaf-qmaf)>.2),
    }

def run(index,inputs,out):
    if out.exists():raise FileExistsError("REFUSE_OLD_OUTPUT")
    rows=data(index)
    keys=set()
    if len(rows)!=646:raise ValueError("EXPECTED_ORIGINAL_646_SOURCE_INDEX")
    results=[]
    for rec in rows:
        locus,dataset,gene=tuple(rec[k] for k in ("locus","dataset_key","gene_base"))
        key=(locus,dataset,gene)
        if key in keys:raise ValueError("DUPLICATED_ORIGINAL_COMPOSITE_KEY")
        keys.add(key)
        name="__".join(key)+".tsv"
        if Path(rec["file"]).name!=name:raise ValueError("INDEX_NAME_MISMATCH")
        p=inputs/name
        item=dict(locus=locus,dataset_key=dataset,gene_base=gene,
                  source_filename=name,status="MISSING_INPUT",reason="LOCAL_SNP_TSV_NOT_CACHED",
                  n_rows=0,n_unique_match_keys=0,n_palindromic=0,
                  n_palindromic_eaf_in_0_4_to_0_6=0,
                  n_maf_diff_gt_0_1=0,n_maf_diff_gt_0_2=0,
                  source_sha256="",
                  exact_GRCh38_GTeX_ID_tag="NOT_ASSESSABLE",
                  GWAS_effect_allele_direction="EXTERNAL_ORIGINAL_GWAS_SCHEMA_REQUIRED",
                  GRCh37_GRCh38_liftover="NOT_INDEPENDENTLY_VERIFIED",
                  GTEx_ancestry_matched_LD="NOT_PROVIDED")
        if p.exists():
            try:
                snps=data(p)
                if len(snps)!=int(rec["nsnps"]):
                    raise ValueError("ROW_COUNT_NOT_INDEX")
                if any((z["locus"],z["dataset_key"],z["gene_base"])!=key for z in snps):
                    raise ValueError("INTERNAL_SOURCE_KEY_MISMATCH")
                if len({z["match_key"] for z in snps})!=len(snps):
                    raise ValueError("DUPLICATED_GRCH38_SNP_KEY")
                flags=[parse_record(z) for z in snps]
                if len(flags)!=len(snps):raise ValueError("INCOMPLETE_SOURCE_AUDIT")
                item.update(status="PASS_STRUCTURAL_ONLY",reason="",
                            n_rows=len(snps),n_unique_match_keys=len(snps),
                            n_palindromic=sum(z["is_palindromic"] for z in flags),
                            n_palindromic_eaf_in_0_4_to_0_6=sum(z["is_palindromic_ambig_eaf"] for z in flags),
                            n_maf_diff_gt_0_1=sum(z["is_gwas_gtEx_maf_diff_gt_0_1"] for z in flags),
                            n_maf_diff_gt_0_2=sum(z["is_gwas_gtEx_maf_diff_gt_0_2"] for z in flags),
                            source_sha256=filehash(p),
                            exact_GRCh38_GTeX_ID_tag="EXACT_REF_ALT_GRCH38_REPRODUCED")
            except Exception as e:
                item.update(status="FAIL",reason=str(e)[:300])
        results.append(item)
    if len(results)!=646:raise ValueError("DENOMINATOR_LOST")
    counts=dict(collections.Counter(z["status"] for z in results))
    assessed=[z for z in results if z["status"]=="PASS_STRUCTURAL_ONLY"]
    overall=dict(
        schema="IS_PHASE10_11_G1_SOURCE_ALLELE_STRUCTURAL_V1",
        total_indexed=646,
        status_counts=counts,
        total_gene_snp_records=sum(z["n_rows"] for z in assessed),
        total_palindromic=sum(z["n_palindromic"] for z in assessed),
        total_palindromic_in_eaf_0_4_to_0_6=sum(z["n_palindromic_eaf_in_0_4_to_0_6"] for z in assessed),
        total_gwas_gtEx_maf_diff_gt_0_1=sum(z["n_maf_diff_gt_0_1"] for z in assessed),
        total_gwas_gtEx_maf_diff_gt_0_2=sum(z["n_maf_diff_gt_0_2"] for z in assessed),
        qtl_GRCh38_variant_internal_structure="PASS" if not counts.get("FAIL") else "FAIL",
        gwas_effect_allele_direction="UNVERIFIED_NEEDS_ORIGINAL_BBJ_GWAS_EFFECT_ALLELE_CODEBOOK",
        gwas_GRCh37_to_qtl_GRCh38_liftover="UNVERIFIED_NEEDS_INDEPENDENT_REFERENCE_MAPPING",
        gtEx_cohort_matched_LD="UNVERIFIED",
        causal_IS_mechanism="NOT_ESTABLISHED",
        interpretation="Repeated gene-SNP rows, not unique SNP count or distinct association signals")
    out.mkdir(parents=True)
    with (out/"IS_G1_SOURCE_STRUCTURAL_STATUS_646.tsv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter="\t")
        w.writeheader();w.writerows(results)
    (out/"IS_G1_SOURCE_STRUCTURAL_SUMMARY.json").write_text(json.dumps(overall,indent=2,ensure_ascii=False)+"\n")
    (out/"IS_G1_METHOD_BOUNDARY.md").write_text(
        "# IS GWAS/GTEx source G1 structural harmonization QC\n\n"
        f"- Assessed {counts.get('PASS_STRUCTURAL_ONLY',0)} of 646 source assays from local cache; unassessed {counts.get('MISSING_INPUT',0)} are not interpreted as PASS.\n"
        f"- Input gene–SNP rows assessed: {overall['total_gene_snp_records']} (SNPs repeated across tissues/genes).\n"
        "- GRCh38 match_key equals original GTEx variant encoding and GRCh38 REF/ALT fields; GTEx reported effect/SE/MAF preserved; GWAS MAF is internally consistent with GWAS EAF.\n"
        "- Nonzero number of palindromic or ancestry-differing MAF variants is NOT itself evidence of allele flip. Neither ref/alt identity nor GWAS MAF demonstrates original GWAS effect-allele coding.\n"
        "- Original BBJ GWAS variant-specific effect allele and independently verified GRCh37->38 liftover remain external missing evidence. Joint GTEx ancestry matched molecular LD also remains unavailable.\n"
        "- SOURCE G1 STRUCTURAL PASS does not imply causal gene validity, reciprocal allelic effects, or cohort genotype harmonization.\n",
        encoding="utf8")
    print(json.dumps(overall,ensure_ascii=False))
    return overall

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--index",type=Path,required=True)
    p.add_argument("--inputs",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.index,a.inputs,a.out)
