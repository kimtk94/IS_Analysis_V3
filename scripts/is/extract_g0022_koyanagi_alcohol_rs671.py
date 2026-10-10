#!/usr/bin/env python3
"""Extract rs671 / 4 G0022 CS variants from VERIFIED Japanese alcohol GWAS.

Strict official Zenodo v1 MD5, CHR/POS hg19, EA/NEA, EAF, BETA, SE, P, N.
No allele inference from SNP rsID; explicit ALT alignment only. This returns
genetic exposure associations, NEVER mediated stroke effect or MR proof.
"""
import argparse,csv,gzip,hashlib,json,math
from pathlib import Path
ROOT=Path("/srv/is-analysis/data/is/gwas_exposure/japanese_alcohol_2024")
OUT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
VARIANTS={
 ("12","112168009"):("G","A"),
 ("12","112241766"):("G","A"),
 ("12","112468206"):("C","T"),
 ("12","112736118"):("A","G"),
}
EXPECTED=("SNP","CHR","POS","EA","NEA","EAF","BETA","SE","P","HetP","N")
FILENAMES={
 "alcohol_intake":"1_Alcohol_intake_Unstratified.tsv.gz",
 "drinking_status":"5_Drinking_Unstratified.tsv.gz",
}
def source_md5(path):
    h=hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1048576),b""):h.update(block)
    return h.hexdigest()
def read_source(path,expected,minimum_rows=100000):
    if not path.is_file() or path.stat().st_size!=expected["size"]:
        raise ValueError("FULL_SOURCE_SIZE_NOT_VERIFIED")
    if source_md5(path)!=expected["checksum"][4:]:
        raise ValueError("FULL_SOURCE_MD5_NOT_VERIFIED")
    matching={};n=0
    with gzip.open(path,"rt") as f:
        reader=csv.DictReader(f,delimiter="\t")
        if tuple(reader.fieldnames or ())!=EXPECTED:
            raise ValueError("Unexpected primary alcohol GWAS header; manual schema audit")
        for row in reader:
            n+=1
            ch=row["CHR"].removeprefix("chr")
            pos=row["POS"]
            key=(ch,pos)
            if key not in VARIANTS:continue
            pair=VARIANTS[key];ea=row["EA"].upper();nea=row["NEA"].upper()
            if (ea,nea)==(pair[1],pair[0]):sign=1
            elif (ea,nea)==(pair[0],pair[1]):sign=-1
            else:raise ValueError("rs671/G0022 exposure source allele discordance")
            if key in matching:raise ValueError("Duplicate chromosome-position exposure row")
            beta=float(row["BETA"]);se=float(row["SE"])
            ef=float(row["EAF"]);p=float(row["P"]);sample=int(float(row["N"]))
            hp=row["HetP"]
            if not (math.isfinite(beta) and math.isfinite(se) and
                math.isfinite(ef) and math.isfinite(p) and
                0<ef<1 and se>0 and 0<=p<=1 and sample>0):
                raise ValueError("Invalid source exposure GWAS statistics")
            if hp not in ("",".","NA","nan") and not (0<=float(hp)<=1):
                raise ValueError("Invalid HetP field")
            matching[key]={"variant_grch37":f"{ch}:{pos}:{pair[0]}:{pair[1]}",
                "source_SNP":row["SNP"],"source_effect_allele":ea,
                "source_non_effect_allele":nea,"harmonized_effect_allele_ALT":pair[1],
                "beta_ALT":sign*beta,"se":se,"p":p,"ALT_EAF":ef if sign==1 else 1-ef,
                "metaGWAS_variant_n":sample,"heterogeneity_p":hp,
                "study_genome_build":"GRCh37_hg19",
                "allele_direction":"EXACT_EA_NEA"}
    if n<minimum_rows:raise ValueError("Unexpectedly truncated source row count")
    return matching,n
def extract(root,out,source_type,minimum_rows=100000):
    meta=json.loads((root/"ZENODO_10038152_METADATA.json").read_text())
    filename=FILENAMES[source_type]
    index={x["key"]:x for x in meta["files"]}
    expected=index[filename]
    selected,n=read_source(root/filename,expected,minimum_rows)
    records=[]
    for key,pair in VARIANTS.items():
        if key not in selected:continue
        record=dict(selected[key]);record["exposure"]=source_type
        record["source"]="Koyanagi_ScienceAdvances2024_open_Zenodo_v1"
        record["source_md5"]=expected["checksum"]
        record["BBJ_overlap_risk"]="YES_BBJ_PARTICIPANTS_INCLUDED_IN_SOURCE_META"
        record["validation_status"]="EXPOSURE_ALLELE_ASSOCIATION_ONLY"
        records.append(record)
    out.mkdir(parents=True,exist_ok=True)
    name=out/f"G0022_4CS_KOYANAGI2024_{source_type.upper()}_EXPOSURE.tsv"
    if records:
        with name.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t")
            w.writeheader();w.writerows(records)
    result={"exposure":source_type,"source":filename,
      "official_MD5_verified":True,"GWAS_meta_summary_rows":n,
      "4CS_GWAS_variants_harmonized":len(records),
      "rs671_present":any(x["variant_grch37"]=="12:112241766:G:A" for x in records),
      "trait_units":"log2(g/day+1)" if source_type=="alcohol_intake" else "binary_status_log_odds_CASE_ENCODING_NEEDS_VALIDATION",
      "effect_allele":"ALT",
      "causal_stroke_alcohol_mediation_computed":False}
    (out/f"G0022_KOYANAGI2024_{source_type.upper()}_AUDIT.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    return result
if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--root",type=Path,default=ROOT)
    parser.add_argument("--out",type=Path,default=OUT)
    parser.add_argument("--exposure",required=True,choices=list(FILENAMES))
    parser.add_argument("--min-source-rows",type=int,default=100000)
    args=parser.parse_args()
    extract(args.root,args.out,args.exposure,args.min_source_rows)
