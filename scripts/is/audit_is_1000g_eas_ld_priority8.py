#!/usr/bin/env python3
"""Audit 1000G EAS504 GWAS LD variant matching for original 8 BBJ locus assays.

LD file names 'BBJ_IS_*' are locus labels, NOT evidence the genotype
reference is BBJ or JPT. The PLINK log describes 504 reference samples.
Do NOT claim LD matrix validation for GTEx donor molecular QTL.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

def read(path):
    with path.open(newline="",encoding="utf8") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(4*1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def run(input_dir,priority8,ld_root,out):
    if out.exists():
        raise FileExistsError("NO_OVERWRITE")
    summary=read(priority8)
    if len(summary)!=8 or len(set(x["gene"] for x in summary))!=8:
        raise ValueError("PRIORITY_8_AUDIT_EXPECTED")
    output=[]
    source_hash={}
    ld_hash={}
    all_rows=[]
    for a in summary:
        gene,locus,tissue=a["gene"],a["locus"],a["dataset_key"]
        ensg={"FGF5":"ENSG00000138675","CALHM2":"ENSG00000138172",
              "NEURL1":"ENSG00000107954","C4orf22":"ENSG00000197826",
              "INA":"ENSG00000148798","SH3PXD2A":"ENSG00000107957",
              "COL4A2":"ENSG00000134871","COL4A1":"ENSG00000187498"}[gene]
        file=input_dir/f"{locus}__{tissue}__{ensg}.tsv"
        snps=read(file)
        if len(snps)!=int(a["nsnps"]):raise ValueError("SOURCE_SNP_COUNT_CHANGED_"+gene)
        varsfile=ld_root/f"{locus}.unphased.vcor1.bin.vars"
        matfile=ld_root/f"{locus}.unphased.vcor1.bin"
        plink_log=ld_root/f"{locus}.log"
        psam=ld_root.parent/f"ld/{locus}.matched.psam"
        for p in (varsfile,matfile,plink_log,psam):
            if not p.exists():raise FileNotFoundError(p)
        vars_=varsfile.read_text(encoding="utf8").splitlines()
        if len(vars_)!=len(set(vars_)) or len(vars_)<100:
            raise ValueError("LD_DUPLICATE_OR_SHORT")
        n=len(vars_)
        if matfile.stat().st_size!=n*n*8:
            raise ValueError("LD_MATRIX_BINARY_DIMENSION_FAIL")
        if "504 samples" not in plink_log.read_text(encoding="utf8",errors="replace"):
            raise ValueError("REFERENCE_SAMPLE_LOG_NOT_504")
        npsam=sum(1 for line in psam.read_text(encoding="utf8").splitlines() if line and not line.startswith("#"))
        if npsam!=504:
            raise ValueError("REFERENCE_PSAM_NOT_504")
        index={v:i+1 for i,v in enumerate(vars_)}
        s37=[x["variant_id"] for x in snps]
        s38=[x["match_key"] for x in snps]
        if len(s37)!=len(set(s37)) or len(s38)!=len(set(s38)):
            raise ValueError("SOURCE_SNP_IDS_DUPLICATED")
        in37=[x in index for x in s37]
        in38=[x in index for x in s38]
        if not all(in37):
            raise ValueError("BBJ_REGION_EAS_LD_SNP_MISSING_"+gene)
        if any(in38):
            raise ValueError("GRCH38_ISSUE_UNEXPECTED_MIXED_LD_KEYS")
        if any(x["variant_id_hg38"]!=x["match_key"] for x in snps):
            raise ValueError("GRCH38_MATCH_KEY_INCONSISTENT")
        g_lead=min(snps,key=lambda x:float(x["gwas_p"]))
        q_lead=min(snps,key=lambda x:float(x["eqtl_p"]))
        if g_lead["variant_id"] not in index or q_lead["variant_id"] not in index:
            raise ValueError("GWAS_OR_EQTL_LEAD_NOT_IN_REFERENCE_LD")
        labels=sorted({x["harmonization"] for x in snps})
        if "" in labels:
            raise ValueError("MISSING_HARMONIZATION_STATUS")
        prev=ld_hash.get(locus)
        if prev is None:
            ld_hash[locus]=dict(variants_sha256=sha(varsfile),matrix_sha256=sha(matfile),
                                sample_psam_sha256=sha(psam),n_variants=n,n_reference_samples=504)
        source_hash[gene]=sha(file)
        output.append(dict(gene=gene,locus=locus,assay=tissue,
            n_input=len(snps),n_reference_ld_variants=n,n_EAS_LD_match_GRCh37=sum(in37),
            n_EAS_LD_match_GRCh38=sum(in38),
            gwas_lead_in_EAS_LD=True,eqtl_lead_in_EAS_LD=True,
            ld_ancestry="1000G_EAS_504_REFERENCE",
            GWAS_LD_input_status="EAS_LD_SNP_COVERAGE_PASS",
            GTEx_QTL_LD_status="GTEX_MATCHED_LD_NOT_PROVIDED_IN_THIS_AUDIT",
            multi_signal_coloc_status="NOT_READY_WITHOUT_QTL_LD_AND_ADDITIONAL_QC",
            harmonization_labels=";".join(labels)))
        for i,x in enumerate(snps):
            all_rows.append(dict(gene=gene,locus=locus,
                source_row=i+1,
                GWAS_LD_id_GRCh37=x["variant_id"],
                QTL_match_key_GRCh38=x["match_key"],
                EAS_LD_matrix_one_based_index=index[x["variant_id"]],
                source_harmonization=x["harmonization"]))
    out.mkdir(parents=True)
    with (out/"IS_PRIORITY8_EAS_LD_READINESS.tsv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(output[0]),delimiter="\t")
        w.writeheader();w.writerows(output)
    with (out/"IS_PRIORITY8_EAS_LD_VARIANT_ORDER_MAP.tsv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(all_rows[0]),delimiter="\t")
        w.writeheader();w.writerows(all_rows)
    report=dict(schema="IS_PRIORITY8_EAS_LD_SOURCE_MATCH_V1",status="GRCH37_EAS504_SNP_COVERAGE_PASS",
        n_gene_tissue_tests=8,n_locus_panels=len(ld_hash),
        n_source_rows=len(all_rows),
        reference_ancestry="1000_GENOMES_EAS_504_NOT_BBJ_COHORT_SPECIFIC_NOT_JPT_ONLY",
        GWAS_ld_panel_status="SNP_COVERAGE_VERIFIED_ALLELE_ORIENTATION_AND_LD_COMPATIBILITY_NOT_FULLY_AUDITED",
        QTL_ld_panel_status="GTEX_MOLECULAR_COHORT_LD_NOT_PROVIDED",
        matched_ancestry_joint_susie_status="BLOCKED",
        source_input_sha256=source_hash,
        ld_panel_sha256=ld_hash,
        results=output)
    (out/"IS_PRIORITY8_EAS_LD_READINESS.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    (out/"IS_PRIORITY8_LD_INTERPRETATION.md").write_text(
      "# EAS LD matching for eight original GTEx coloc tests\n\n"
      "- Verified 1000 Genomes **EAS 504**, from PLINK logs/PSAM. Files are named by BBJ loci but are not BBJ cohort genotype data nor Japan-only JPT.\n"
      "- All eight shared SNP inputs map completely to EAS matrix by original GWAS **GRCh37 variant_id**, not the GRCh38 match_key.\n"
      "- Input REF/ALT or source harmonization field alignment remains as recorded by upstream workflow; this audit does not independently validate allele-flip semantics.\n"
      "- Each per-locus binary matrix size agrees with its square double-precision variant count. Binary numerical contents are not inspected here.\n"
      "- GTEx molecular-cohort matched LD matrix has not been provided/verified; therefore matched GWAS-QTL LD multi-signal SuSiE **BLOCKED**.\n"
      "- EAS reference SNP coverage alone does not verify matched-ancestry population LD, causal variant concordance, or tissue-specific gene mediation.\n",
      encoding="utf8")
    print("IS_PRIORITY8_EAS_LD_MATCH_PASS",
          "CASES=",len(output),"VARIANTS=",len(all_rows),
          "PANELS=",len(ld_hash),"GTEX_QTL_LD_STATUS=NOT_VERIFIED")
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for t in ("input_dir","priority8","ld_root","out"):
        p.add_argument("--"+t,type=Path,required=True)
    a=p.parse_args()
    run(a.input_dir,a.priority8,a.ld_root,a.out)
