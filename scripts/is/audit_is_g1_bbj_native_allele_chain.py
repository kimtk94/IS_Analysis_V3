#!/usr/bin/env python3
"""Read-only independent BBJ native GWAS Allele2/BETA/EAF + 37->38 chain audit.

Scopes ONLY locally available original-646 coloc input files (currently 30).
No claim that GRCh38 GTEx effect-allele coding is experimentally verified.
The original BBJ .zip contains a gz stream; stream and filter by SNP IDs.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
import sys
import zipfile
from pathlib import Path

FIELDS=["locus","dataset_key","gene_base","variant_id","match_key",
        "source_gwas_beta","native_gwas_beta","source_gwas_eaf",
        "native_af_allele2","native_allele1","native_allele2",
        "native_effect_allele","gwas_id_ref","gwas_id_alt",
        "gwas_effect_vs_id","chain_target","chain_strand",
        "chain_mapping_count","target_qtl_refalt","status","reason"]

def float_eq(a,b,tol=1e-9):
    return math.isfinite(a) and math.isfinite(b) and abs(a-b)<=tol*max(1,abs(a),abs(b))

def unique_sources(input_dir):
    ids={}
    repeated=0
    for p in sorted(input_dir.glob("BBJ_IS_L00*__GTEx_V8__*.tsv")):
        with p.open(newline="",encoding="utf8") as h:
            for r in csv.DictReader(h,delimiter="\t"):
                ident=r["variant_id"]
                if ident in ids:
                    q=ids[ident]
                    if any(r[k]!=q[k] for k in ("gwas_beta","gwas_eaf","match_key")):
                        raise ValueError("SAME_GWAS_SNP_INCONSISTENT_CACHED_SOURCE:"+ident)
                    repeated+=1
                else:
                    ids[ident]=r
    if not ids:raise ValueError("NO_CACHED_SOURCE_INPUTS")
    return ids,repeated

def load_native(zip_file,targets):
    found={}
    stream_rows=0
    with zipfile.ZipFile(zip_file) as arc:
        matches=[n for n in arc.namelist() if n.endswith(".auto.txt.gz")]
        if len(matches)!=1:raise ValueError("EXPECTED_SINGLE_AUTOSOME_GWAS_SOURCE")
        with arc.open(matches[0]) as compressed:
            with gzip.open(compressed,"rt") as handle:
                reader=csv.DictReader(handle,delimiter="\t")
                required=("v","CHR","POS","Allele1","Allele2","AF_Allele2","BETA")
                if any(k not in reader.fieldnames for k in required):
                    raise ValueError("NATIVE_BBJ_GWAS_SCHEMA_INVALID")
                for r in reader:
                    stream_rows+=1
                    ident=r["v"]
                    if ident not in targets:continue
                    if ident in found:raise ValueError("DUPLICATED_NATIVE_BBJ_SNP:"+ident)
                    found[ident]={k:r[k] for k in required}
    return found,stream_rows

def score(source,native,lo):
    key=source["variant_id"]
    chrom,pos,ref,alt=key.split(":")
    ref,alt=ref.upper(),alt.upper()
    out={k:"" for k in FIELDS}
    for k in ("locus","dataset_key","gene_base","variant_id","match_key"):
        out[k]=source[k]
    out.update(source_gwas_beta=source["gwas_beta"],
               source_gwas_eaf=source["gwas_eaf"],
               gwas_id_ref=ref,gwas_id_alt=alt)
    if native is None:
        out.update(status="NOT_FOUND_NATIVE",reason="No variant_id in original autosome GWAS")
        return out
    out.update(native_gwas_beta=native["BETA"],
               native_af_allele2=native["AF_Allele2"],
               native_allele1=native["Allele1"],
               native_allele2=native["Allele2"],
               native_effect_allele=native["Allele2"])
    # Allele2 is the documented native effect allele: BETA and AF_Allele2.
    if native["CHR"]!=chrom or native["POS"]!=pos:
        out.update(status="FAIL",reason="NATIVE_BBJ_CHROM_POS_NOT_SOURCE_ID")
        return out
    a1,a2=native["Allele1"].upper(),native["Allele2"].upper()
    if {a1,a2}!={ref,alt}:
        out.update(status="FAIL",reason="NATIVE_BBJ_ALLELE_SET_NOT_GRCH37_ID")
        return out
    if a2==alt:out["gwas_effect_vs_id"]="ALT"
    elif a2==ref:out["gwas_effect_vs_id"]="REF"
    else:out["gwas_effect_vs_id"]="NEITHER"
    if not float_eq(float(source["gwas_beta"]),float(native["BETA"]),tol=1e-9):
        out.update(status="FAIL",reason="GWAS_BETA_DIVERGES_FROM_NATIVE_ALLELE2")
        return out
    if not float_eq(float(source["gwas_eaf"]),float(native["AF_Allele2"]),tol=1e-9):
        out.update(status="FAIL",reason="GWAS_EAF_DIVERGES_FROM_NATIVE_ALLELE2")
        return out
    mapped=lo.convert_coordinate("chr"+chrom,int(pos)-1)
    out["chain_mapping_count"]=str(len(mapped) if mapped else 0)
    if not mapped or len(mapped)!=1:
        out.update(status="FAIL",reason="NOT_UNIQUE_GRCH37_TO_GRCH38_CHAIN")
        return out
    c,p,strand,_=mapped[0]
    out["chain_target"]=f"{c.removeprefix('chr')}:{p+1}"
    out["chain_strand"]=strand
    chr38,pos38,ref38,alt38=source["match_key"].split(":")
    if c.removeprefix("chr")!=chr38 or p+1!=int(pos38):
        out.update(status="FAIL",reason="CHAIN_GRCH38_POSITION_NOT_QTL_MATCH_KEY")
        return out
    if strand=="+":
        er,ea=ref,alt
    elif strand=="-":
        comp=str.maketrans("ACGT","TGCA")
        er,ea=ref.translate(comp)[::-1],alt.translate(comp)[::-1]
    else:
        out.update(status="FAIL",reason="UNSUPPORTED_CHAIN_STRAND")
        return out
    out["target_qtl_refalt"]=f"{ref38}:{alt38}"
    if er!=ref38 or ea!=alt38:
        out.update(status="FAIL",reason="QTL_GRCH38_REF_ALT_NOT_CHAIN_STRAND_CORRESPONDING")
        return out
    out.update(status="PASS_NATIVE_EFFECT_ALLELE_PLUS_CHAIN",
               reason="Native BBJ effect Allele2 BETA/EAF and 37_to_38 chain concordant")
    return out

def run(inputs,zip_file,chain,library,out):
    if out.exists():raise FileExistsError("REFUSE_EXISTING_OUTPUT")
    if not all(p.is_file() for p in (zip_file,chain)):
        raise FileNotFoundError("NATIVE_GWAS_OR_CHAIN_MISSING")
    ids,repeated=unique_sources(inputs)
    sys.path.insert(0,str(library))
    from pyliftover import LiftOver
    native,scanned=load_native(zip_file,set(ids))
    lo=LiftOver(str(chain))
    rows=[score(v,native.get(k),lo) for k,v in sorted(ids.items())]
    statuses=collections.Counter(x["status"] for x in rows)
    effect=collections.Counter(x["gwas_effect_vs_id"] for x in rows)
    strand=collections.Counter(x["chain_strand"] for x in rows)
    d={
       "schema":"IS_G1_BBJ_NATIVE_ALLELE_AND_CHAIN_INDEPENDENT_V1",
       "scope":"ORIGINAL_646_LOCAL_CACHE_DISTINCT_GRCH37_GWAS_SNP",
       "unique_GRCh37_SNP_count":len(ids),
       "duplicate_gene_tissue_SNP_records":repeated,
       "native_GWAS_full_stream_row_count":scanned,
       "native_SNP_found_count":len(native),
       "status_counts":dict(statuses),
       "native_effect_allele_vs_GRCh37_id":dict(effect),
       "chain_strand_counts":dict(strand),
       "source_zip_bytes":zip_file.stat().st_size,
       "source_zip_name":zip_file.name,
       "chain_file_sha256":hashlib.sha256(chain.read_bytes()).hexdigest(),
       "gwas_direction_status":"NATIVE_BBJ_ALLELE2_BETA_EAF_CHECKED" if all(x["status"]=="PASS_NATIVE_EFFECT_ALLELE_PLUS_CHAIN" for x in rows) else "SOME_NOT_VERIFIED",
       "reference_base_status":"GRCh38_ACTUAL_REFERENCE_BASES_NOT_INDEPENDENTLY_CHECKED",
       "GTEx_effect_beta_allele_status":"REQUIRES_OFFICIAL_GTEX_V8_EFFECT_ALLELE_CONVENTION",
       "sample_ancestry_LD":"NOT_GTEX_COHORT_MATCHED",
       "causal_verdict":"NOT_ESTABLISHED",
    }
    out.mkdir(parents=True)
    with (out/"IS_G1_NATIVE_BBJ_ALLELE_CHAIN_SNP.tsv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter="\t")
        w.writeheader();w.writerows(rows)
    (out/"IS_G1_NATIVE_BBJ_ALLELE_CHAIN_SUMMARY.json").write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n")
    (out/"IS_G1_NATIVE_BBJ_CHAIN_BOUNDARY.md").write_text(
        "# Native BBJ GWAS effect allele + chain source audit\n\n"
        "- The original hum0197 GWAS IS autosomal native archive is parsed directly: Allele2 is the source BETA/AF_Allele2 allele.\n"
        "- Each unique cached GRCh37 variant is checked for native beta and AF values, allele-set correspondence, single-chain GRCh37->38 coordinate mapping and oriented REF/ALT correspondence to GTEx match_key.\n"
        "- Native GWAS effect allele can differ from GRCh37 ALT, and this distinction must not be interpreted as a beta sign flip without checking the original harmonization conventions.\n"
        "- Genomic coordinate chain is not a GRCh38 reference-sequence allele validation; official GTEx v8 beta effect-allele convention is a separate validation gate.\n"
        "- This audit of cached original-646 source variants does not verify all 646 separate gene-tissue source TSVs.\n",
        encoding="utf8")
    print(json.dumps(d,ensure_ascii=False))
    return d

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for arg in ("inputs","zip_file","chain","library","out"):
        p.add_argument("--"+arg,type=Path,required=True)
    a=p.parse_args()
    run(a.inputs,a.zip_file,a.chain,a.library,a.out)
