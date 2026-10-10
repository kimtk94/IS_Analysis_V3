#!/usr/bin/env python3
"""Recover Korean KCPS2 rs671 in selectively acquired SBP / DBP GWAS.

Read actual SAIGE source, validate member CRC32 + stored SHA256, effect of
Allele2 (SAIGE convention), harmonize to rs671-A; source imputed build
inferred only at unambiguous reciprocal GRCh37/GRCh38 rs671 coordinate.
BP effects are inverse-rank-normal-transformed traits, not raw mmHg.
"""
import argparse,csv,gzip,json,math
from pathlib import Path
from extract_kcps2_korean_rs671_alcohol import DATA,RESULT,EXPECTED,LOOKUP,sha256
FILES={"sbp":"SAIGE_SBP_INFO.txt.gz","dbp":"SAIGE_DBP_INFO.txt.gz"}
def extract(root,out,trait,minrows=100000):
    name=FILES[trait];source=root/name;manifest=root/(name+".member_crc_verified.json")
    x=json.loads(manifest.read_text())
    if x["archive_member_name"]!=name or x.get("status")!="SELECTIVE_ZIP_MEMBER_CRC_VERIFIED_ARCHIVE_MD5_NOT_COMPUTED":
        raise ValueError("KCPS2 original ZIP member CRC gate not met")
    if not source.is_file() or source.stat().st_size!=x["ZIP_member_size_bytes"] or sha256(source)!=x["inner_gzip_sha256"]:
        raise ValueError("KCPS2 member integrity mismatch")
    observations=[];n=0
    with gzip.open(source,"rt") as f:
        reader=csv.DictReader(f,delimiter="\t")
        if tuple(reader.fieldnames or ())!=EXPECTED:raise ValueError("KCPS2 SAIGE schema mismatch")
        for r in reader:
            n+=1
            key=(r["CHR"].removeprefix("chr"),r["POS"])
            if key not in LOOKUP:continue
            a=r["Allele1"].upper();b=r["Allele2"].upper()
            if set((a,b))!={"G","A"}:raise ValueError("rs671 alleles not matched")
            beta=float(r["BETA"]);se=float(r["SE"]);p=float(r["p.value"]);af=float(r["AF_Allele2"])
            if not all(map(math.isfinite,(beta,se,p,af))) or se<=0 or not(0<af<1 and 0<=p<=1):
                raise ValueError("Invalid Korean BP GWAS effect")
            sign=1 if b=="A" else -1
            observations.append({"trait":trait.upper(),
                "variant":"12:112241766:G:A",
                "source_chr":key[0],"source_pos":key[1],
                "rs671_genome_build_at_position":LOOKUP[key],
                "SAIGE_effect_allele":"Allele2","original_allele2":b,
                "harmonized_effect_allele":"A",
                "KCPS2_beta_A":sign*beta,"KCPS2_se":se,
                "KCPS2_p":p,"KCPS2_neglogP_raw":r["neglogP"],
                "KCPS2_per_variant_N":int(float(r["N"])),
                "KCPS2_af_A":af if sign==1 else 1-af,
                "BP_units":"INVERSE_RANK_NORMAL_TRANSFORM_NOT_MMHG",
                "original_zip_full_md5_tested":False,
                "source_selected_member_crc_pass":True,
                "causal_pathway_estimate":"NOT_ESTABLISHED"})
    if n<minrows or len(observations)!=1:
        raise ValueError(f"KCPS2 incomplete or rs671 ambiguous source rows={n} variants={len(observations)}")
    out.mkdir(parents=True,exist_ok=True)
    p=out/f"G0022_RS671_KCPS2_KOREAN_{trait.upper()}_GWAS.tsv"
    with p.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(observations[0]),delimiter="\t")
        w.writeheader();w.writerows(observations)
    o=observations[0]
    result={"source":"KCPS2_Jee_2025_Zenodo15132424",
        "trait":trait.upper(),"source_member_name":name,"source_genome_build":o["rs671_genome_build_at_position"],
        "source_member_SAIGE_records":n,"selected_member_SHA256_verified":True,
        "selected_member_CRC32_verified":True,"entire_original_ZIP_MD5_verified":False,
        "rs671_effect_A_beta":o["KCPS2_beta_A"],"rs671_effect_se":o["KCPS2_se"],
        "rs671_effect_p":o["KCPS2_p"],"rs671_effect_n":o["KCPS2_per_variant_N"],
        "rs671_A_frequency":o["KCPS2_af_A"],
        "mmHg_interpretation_permitted":False,"alcohol_mediation_inference":False,
        "status":"REAL_KOREAN_BP_RS671_VARIANT_ASSOCIATION_ONLY"}
    (out/f"G0022_RS671_KCPS2_KOREAN_{trait.upper()}_GWAS_SUMMARY.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2));return result
if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--source",type=Path,default=DATA)
    a.add_argument("--out",type=Path,default=RESULT)
    a.add_argument("--trait",required=True,choices=list(FILES))
    a.add_argument("--min-source-rows",type=int,default=100000)
    opt=a.parse_args();extract(opt.source,opt.out,opt.trait,opt.min_source_rows)
