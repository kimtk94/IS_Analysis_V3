#!/usr/bin/env python3
"""Audit actual nominal eQTL Catalogue BLUEPRINT public data rows (GRCh38).

Read a small bounded prefix only, validate ALT-as-effect reference schema,
tested SNP alleles/position, beta/se/MAF and readable nonsignificant rows.
NEVER call this full cis region completeness, GWAS harmonization, or coloc.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import io
import json
import math
import os
from pathlib import Path
import re
import statistics
import zlib

from audit_eur_eqtl_public_bytes import EBI_BLUEPRINT_URL,remote_prefix

REQUIRED=(
    "variant","pvalue","molecular_trait_id","maf","gene_id","beta","se",
    "chromosome","position","ref","alt","an","ac",
)
VALID_NUCLEOTIDES={"A","C","G","T"}


def parse_nominal_rows(gzip_bytes, max_rows=300):
    if not gzip_bytes.startswith(b"\x1f\x8b"):
        raise ValueError("Not original gzip eQTL bytes")
    inflater=zlib.decompressobj(16+zlib.MAX_WBITS)
    text=inflater.decompress(gzip_bytes,3*1024*1024).decode("utf-8","strict")
    lines=text.split("\n")
    if len(lines)<2:
        raise ValueError("No full nominal row decoded in source prefix")
    columns=[v.lstrip("#") for v in lines[0].rstrip("\r").split("\t")]
    if len(set(columns))!=len(columns):
        raise ValueError("Duplicate summary fields")
    missing=sorted(set(REQUIRED)-set(columns))
    if missing:
        raise ValueError("Missing original eQTL nominal fields: "+str(missing))
    positions={name:columns.index(name) for name in REQUIRED}
    rows=[]
    for line in lines[1:-1]:
        if len(rows)>=max_rows:
            break
        if not line.strip():
            continue
        fields=line.rstrip("\r").split("\t")
        if len(fields)!=len(columns):
            raise ValueError("Malformed source summary row")
        def fld(key):return fields[positions[key]]
        alleles=(fld("ref"),fld("alt"))
        if any(not x or x.upper() not in VALID_NUCLEOTIDES for x in alleles):
            # Count other variant types separately rather than assuming all SNP
            continue
        try:
            position=int(fld("position"))
            beta=float(fld("beta"))
            se=float(fld("se"))
            p=float(fld("pvalue"))
            maf=float(fld("maf"))
            an=float(fld("an"))
            ac=float(fld("ac"))
        except ValueError as exc:
            raise ValueError("Actual nominal test has invalid numeric value") from exc
        if (position<1 or not all(map(math.isfinite,(beta,se,p,maf,an,ac)))
                or se<=0 or not 0<=p<=1 or not 0<=maf<=0.5 or an<=0 or ac<0 or ac>an):
            raise ValueError("Actual nominal SNP row failed statistical field QC")
        chrom=fld("chromosome")
        expected=f"chr{chrom}_{position}_{fld('ref')}_{fld('alt')}"
        if fld("variant")!=expected:
            raise ValueError("Actual nominal chromosome:position:REF:ALT mismatch")
        if not fld("gene_id").startswith("ENSG"):
            raise ValueError("Non-Ensembl gene identifier")
        rows.append({
            "pvalue":p,"se":se,"beta":beta,
            "variant":fld("variant"),"gene":fld("gene_id"),
        })
    if len(rows)<1:
        raise ValueError("No complete SNP association rows validated in bounded source")
    n_nominal_non_sig=sum(x["pvalue"]>0.05 for x in rows)
    return {
        "assessed_complete_SNP_rows":len(rows),
        "first_prefix_rows_with_p_gt_0_05":n_nominal_non_sig,
        "has_nonsignificant_association_in_prefix":n_nominal_non_sig>0,
        "header":columns,
        "all_checked_rows_have_beta_se_maf_and_ref_alt":True,
        "sample_allele_genome_assembly":"GRCh38 by EBI documentation, independently coordinate unproven",
        "alternate_allele_effect_direction":"ALT per EBI documentation",
        "full_cis_distribution_tested":"NOT_VERIFIED",
        "source_crc_full_file":"NOT_VERIFIED",
        "first_variant_example_only":rows[0]["variant"],
        "first_gene_example_only":rows[0]["gene"],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",default="artifacts/is_ebi_nominal_real_rows")
    args=p.parse_args()
    out=Path(args.out_dir)
    out.mkdir(parents=True,exist_ok=True)
    report={
        "source":EBI_BLUEPRINT_URL,
        "checked_utc":datetime.now(timezone.utc).isoformat(),
        "status":"SOURCE_OR_ROW_QC_HOLD",
        "ancestry_of_BLUEPRINT_donors_verified":"NOT_INDEPENDENTLY_AUDITED",
        "GIGASTROKE_GRCh37_liftover":"NOT_RUN",
        "full_locus_cis_pair_completeness":"NOT_VERIFIED",
        "GWAS_effect_alleles_harmonized":False,
        "coloc":"NOT_RUN",
    }
    try:
        source,meta=remote_prefix(EBI_BLUEPRINT_URL)
        result=parse_nominal_rows(source)
        report.update({"status":"PUBLIC_NOMINAL_QTL_FIRST_ROWS_STRUCTURE_PASS",
                       "download":meta,"actual_source_rows":result})
    except (OSError,ValueError,UnicodeDecodeError,zlib.error) as e:
        report["source_issue"]=f"{type(e).__name__}: {str(e)[:200]}"
    (out/"real_nominal_row_qc.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    lines=["# eQTL Catalogue BLUEPRINT real nominal association record validation","",
           "**"+report["status"]+"**",""]
    if item:=report.get("actual_source_rows"):
        lines.extend([
            f"- SNP rows (bounded beginning of public file): {item['assessed_complete_SNP_rows']}",
            f"- Nonsignificant p>0.05 examples: {item['first_prefix_rows_with_p_gt_0_05']}",
            f"- beta/SE/MAF/REF/ALT present: {item['all_checked_rows_have_beta_se_maf_and_ref_alt']}",
            f"- Original source header: {', '.join(item['header'])}",
        ])
    if e:=report.get("source_issue"):
        lines.append("- issue: "+e)
    lines.extend(["",
        "This is source row-level integrity only. GRCh38 BLUEPRINT ≠ GRCh37 GIGASTROKE; no direct variant/effect allele join.",
        "The complete locus, LD matrix, prior tissue/donor ancestry, and full nominal coverage remain HOLD."
    ])
    md="\n".join(lines)+"\n"
    (out/"summary.md").write_text(md,encoding="utf-8")
    print(md)
    if name:=os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(name,"a",encoding="utf-8") as f:f.write(md)


if __name__=="__main__":
    main()
