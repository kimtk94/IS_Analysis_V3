#!/usr/bin/env python3
"""Ischemic-stroke original 646 coloc panel numerical/allele LD QC (read-only).

Read source SNPs available in the LOCAL cache only. Matrix is PLINK 2.0
--r-unphased square bin ref-based, <f8 matrix in .vars order. This checks
reference-matrix numerics and GRCh37 allele IDs, NOT QTL-matched LD.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("REQUIRES_NUMPY: use validated MasterOmics Python venv") from exc

def sha256(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for z in iter(lambda:f.read(4*1024*1024),b""):
            h.update(z)
    return h.hexdigest()

def table(p):
    with open(p,encoding="utf8",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def parse_pvar(path):
    data={}
    with open(path,encoding="utf8") as f:
        for line in f:
            if line.startswith("##"):continue
            if line.startswith("#CHROM"):continue
            c=line.rstrip("\n").split("\t")
            if len(c)<5:raise ValueError("PVAR_NOT_STANDARD")
            chrom,pos,id_,ref,alt=c[:5]
            if id_ in data:raise ValueError("PVAR_DUPLICATE_ID")
            data[id_]=(chrom,pos,ref,alt)
    return data

def run(input_dir, ldroot, out):
    if out.exists():raise FileExistsError("REFUSE_OUTPUT_OVERWRITE")
    sources=sorted(input_dir.glob("BBJ_IS_L00*__GTEx_V8__*.tsv"))
    if not sources:raise ValueError("NO_SOURCE_TSV_INPUTS")
    matrices={}
    srcrows=[]
    genetic_rows=[]
    for path in sources:
        r=table(path)
        if not r:raise ValueError("EMPTY_ASSAY_"+path.name)
        locus=r[0]["locus"]
        if any(z["locus"]!=locus for z in r):
            raise ValueError("SOURCE_MULTILOCUS_INVALID")
        if locus not in matrices:
            varsfile=ldroot/f"{locus}.unphased.vcor1.bin.vars"
            binfile=ldroot/f"{locus}.unphased.vcor1.bin"
            buildlog=ldroot/f"{locus}.log"
            pvar=ldroot.parent/"ld"/f"{locus}.matched.pvar"
            psam=ldroot.parent/"ld"/f"{locus}.matched.psam"
            files=(varsfile,binfile,buildlog,pvar,psam)
            if any(not p.exists() for p in files):
                raise FileNotFoundError(f"LD_SOURCE_FILES_MISSING:{locus}")
            if "504 samples" not in buildlog.read_text(errors="replace"):
                raise ValueError("LD_504_REFERENCE_UNVERIFIED")
            npsam=sum(bool(x.strip()) and not x.startswith("#") for x in psam.read_text().splitlines())
            if npsam!=504:raise ValueError("PSAM_504_SAMPLE_COUNT_MISMATCH")
            v=varsfile.read_text().splitlines()
            n=len(v)
            if n<100 or len(v)!=len(set(v)) or binfile.stat().st_size!=8*n*n:
                raise ValueError("LD_PANEL_LAYOUT_INVALID")
            mat=np.memmap(binfile,dtype="<f8",mode="r",shape=(n,n))
            if not np.isfinite(mat).all():raise ValueError("NONFINITE_LD_MATRIX")
            diag=np.diag(mat)
            if not np.allclose(diag,1.0,rtol=0,atol=1e-9):
                raise ValueError("BAD_LD_DIAGONAL")
            symerr=float(np.max(np.abs(mat-mat.T)))
            if symerr>1e-8:raise ValueError("LD_NOT_SYMMETRIC")
            lmin,lmax=float(mat.min()),float(mat.max())
            if lmin < -1.00000001 or lmax > 1.00000001:
                raise ValueError("LD_RANGE_OUTSIDE_MINUS1_PLUS1")
            # Deterministic spread across full matrix. Not a whole-matrix PSD certification.
            sample_ids=np.linspace(0,n-1,min(128,n),dtype=int)
            eigs=np.linalg.eigvalsh(mat[np.ix_(sample_ids,sample_ids)])
            eigmin=float(eigs.min())
            if eigmin< -1e-6:
                raise ValueError("LD_SAMPLED_PSD_FAIL")
            pvars=parse_pvar(pvar)
            if len(pvars)!=n:raise ValueError("PVAR_VARS_COUNT_MISMATCH")
            for key in v:
                if key not in pvars:raise ValueError("PVAR_VARS_ID_MISMATCH")
                chrom,pos,ref,alt=pvars[key]
                if key!=f"{chrom}:{pos}:{ref}:{alt}":
                    raise ValueError("PVAR_SOURCE_REF_ALT_ID_MISMATCH")
            matrices[locus]=dict(
                locus=locus,reference="1000G_EAS_504_NOT_COHORT_MATCHED",
                sample_count=504,n_ld_vars=n,diag_min=float(diag.min()),
                diag_max=float(diag.max()),max_abs_asymmetry=symerr,
                ld_min=lmin,ld_max=lmax,
                n_sampled_eigenvalues=len(eigs),min_sampled_eigenvalue=eigmin,
                pvar_refalt_alignment="PASS_FOR_GRCH37_REF_ALT_ID",
                matrix_binary_sha256=sha256(binfile),
                matrix_vars_sha256=sha256(varsfile),
                genotype_pvar_sha256=sha256(pvar),
                genotype_psam_sha256=sha256(psam))
            # Keep lookups memory-local, do not store in auditable report.
            matrices[locus]["_lookup"]=set(v)
            matrices[locus]["_pvar"]=pvars
        snps=matrices[locus]["_lookup"]
        pvars=matrices[locus]["_pvar"]
        found=sum(x["variant_id"] in snps for x in r)
        direct_qtl=sum(x["match_key"] in snps for x in r)
        if found!=len(r):
            raise ValueError("NOT_ALL_SNP_IDS_IN_REFERENCE_LD")
        bad_alleles=0
        palindromic=0
        gwas_maf_invalid=0
        input_harmonization=set()
        for v in r:
            id_=v["variant_id"]
            info=pvars[id_]
            expected=f"{info[0]}:{info[1]}:{info[2]}:{info[3]}"
            if id_!=expected:bad_alleles+=1
            if {info[2],info[3]} in ({"A","T"},{"C","G"}):
                palindromic+=1
            gwas_maf=float(v["gwas_maf"])
            gwas_eaf=float(v["gwas_eaf"])
            if abs(min(gwas_eaf,1-gwas_eaf)-gwas_maf)>1e-7:
                gwas_maf_invalid+=1
            input_harmonization.add(v["harmonization"])
        if bad_alleles or gwas_maf_invalid:
            raise ValueError("GRCH37_ALLELE_OR_EAF_INTERNAL_QC_FAIL")
        genetic_rows.append(dict(
            source_file=path.name,locus=locus,
            dataset=r[0]["dataset_key"],gene=r[0]["gene_base"],
            n_source_snps=len(r),n_matching_GRCh37_LD=found,
            n_matching_GRCh38_keys=direct_qtl,
            n_GRCh37_pvar_REF_ALT_ID_mismatch=bad_alleles,
            n_GWAS_MAF_vs_EAF_mismatch=gwas_maf_invalid,
            palindromic_source_snps=palindromic,
            n_harmonization_label_types=len(input_harmonization),
            allele_effect_orientation_status="UNVERIFIED_UPSTREAM",
            matched_QTL_LD_status="NOT_AVAILABLE"))
        srcrows.append(dict(file=path.name,sha256=sha256(path),
                            n_snps=len(r)))
    panels=[]
    for k,v in sorted(matrices.items()):
        del v["_lookup"],v["_pvar"]
        panels.append(v)
    result={
        "schema":"IS_LEGACY_ORIGINAL_SNP_CACHE_NUMERIC_LD_QC_V1",
        "panel_reference":"1000G_EAS_504_NOT_BBJ_OR_JPT_ONLY",
        "source_file_count":len(srcrows),
        "n_source_SNP_rows":sum(x["n_source_snps"] for x in genetic_rows),
        "n_locus_panels":len(panels),
        "LD_numerical_status":"PASS_FINITE_SYMMETRIC_UNIT_DIAGONAL_SAMPLED_PSD",
        "GRCh37_variant_REF_ALT_status":"PASS",
        "GWAS_effect_allele_alignment_status":"NOT_INDEPENDENTLY_VERIFIED",
        "GTEx_cohort_LD_status":"NOT_OBTAINED",
        "scientific_mechanism_status":"NOT_ESTABLISHED",
        "panels":panels,
        "assays":genetic_rows,
        "source_sha256":srcrows
    }
    out.mkdir(parents=True)
    (out/"IS_646_CACHE_LD_NUMERIC_AUDIT.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    for name,records in (("IS_646_LD_PANELS.tsv",panels),
                         ("IS_646_CACHED_ASSAYS_LD_MATCH.tsv",genetic_rows)):
        with (out/name).open("w",newline="",encoding="utf8") as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t")
            w.writeheader();w.writerows(records)
    (out/"IS_646_LD_NUMERIC_INTERPRETATION.md").write_text(
        "# Direct SNP source cache / 1000G EAS LD numerical QC\n\n"
        f"- {len(srcrows)} cached GTEx original assay TSVs; "
        f"{sum(x['n_source_snps'] for x in genetic_rows)} gene–SNP rows; "
        f"{len(panels)} EAS LD region matrices.\n"
        "- All sources match GRCh37 variant ID and local REF/ALT labelled pvar IDs; "
        "this is an ID agreement check, not verified effect-allele orientation.\n"
        "- Full matrix finite, symmetric, correlation diagonal 1, range [-1,1]. "
        "128 deterministic spread-out eigenvalues show no substantial negative eigenvalue. "
        "This is a sampled, not global PSD guarantee.\n"
        "- Underlying panel is 1000 Genomes EAS504. Not original BBJ GWAS participants, "
        "not Japan-only JPT and not GTEx molecular-QTL donor LD.\n"
        "- Do not use these diagnostics as matched-ancestry, multi-signal SuSiE readiness. "
        "Stroke-specific or causal-regulatory conclusions remain blocked.\n",
        encoding="utf8")
    print("IS_EAS_LD_NUMERICAL_QC_PASS",
          "CACHED_ASSAYS",len(srcrows),
          "SNP_ROWS",result["n_source_SNP_rows"],
          "PANELS",len(panels))
    return result

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--input-dir",type=Path,required=True)
    ap.add_argument("--ld-root",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    run(a.input_dir,a.ld_root,a.out)
