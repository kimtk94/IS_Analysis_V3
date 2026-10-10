#!/usr/bin/env python3
"""Recover KCPS2 Korean rs671 ALCO_AMOUNT association from selectively sourced ZIP member.

The source is one selectively downloaded official ZIP member validated by
original central-directory CRC32 and member SHA256; original 13.2GB ZIP's
full MD5 is NOT available. SAIGE BETA/AF_Allele2 are for Allele2.
Exact GRCh37 and GRCh38 rs671 coordinates inspected; avoid undocumented
assumption that all other loci have identical genome build.
No numeric comparison of transformed phenotypes or causal MR.
"""
import argparse,csv,gzip,hashlib,json,math
from pathlib import Path
DATA=Path("/srv/is-analysis/data/is/gwas_exposure/kcps2_2025_selective")
RESULT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
NAME="SAIGE_ALCO_AMOUNT_INFO.txt.gz"
EXPECTED=("CHR","POS","MarkerID","Allele1","Allele2","AC_Allele2","AF_Allele2",
          "imputationInfo","BETA","SE","Tstat","var","p.value","N","INFO","neglogP","Z")
LOOKUP={("12","112241766"):"GRCh37_hg19",
        ("12","111803962"):"GRCh38_hg38"}
def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1048576),b""):h.update(block)
    return h.hexdigest()
def extract(root,out,minimum_rows=100000):
    infile=root/NAME
    source=json.loads((root/(NAME+".member_crc_verified.json")).read_text())
    if source.get("status")!="SELECTIVE_ZIP_MEMBER_CRC_VERIFIED_ARCHIVE_MD5_NOT_COMPUTED":
        raise ValueError("Unverified KCPS2 selective member provenance")
    if not infile.is_file() or infile.stat().st_size!=source["ZIP_member_size_bytes"]:
        raise ValueError("KCPS2 source archive member incomplete")
    if sha256(infile)!=source["inner_gzip_sha256"]:
        raise ValueError("KCPS2 member SHA256 mismatch")
    selected=[];count=0
    with gzip.open(infile,"rt") as f:
        reader=csv.DictReader(f,delimiter="\t")
        if tuple(reader.fieldnames or ())!=EXPECTED:
            raise ValueError("SAIGE canonical column schema has changed: "+str(reader.fieldnames))
        for r in reader:
            count+=1
            if (r["CHR"].removeprefix("chr"),r["POS"]) not in LOOKUP:
                continue
            ch=r["CHR"].removeprefix("chr");pos=r["POS"]
            a=r["Allele1"].upper();b=r["Allele2"].upper()
            if set((a,b))!={"G","A"}:raise ValueError("rs671 source allele not G/A")
            beta=float(r["BETA"]);se=float(r["SE"]);p=float(r["p.value"])
            af=float(r["AF_Allele2"]);info=float(r["imputationInfo"]);n=int(float(r["N"]))
            if not all(map(math.isfinite,(beta,se,p,af,info))) or se<=0 or not(0<af<1):
                raise ValueError("Invalid rs671 SAIGE GWAS row")
            sign=1 if b=="A" else -1
            selected.append({
               "variant":"ALDH2_rs671_GRCh37_12_112241766_G_A",
               "source_marker_id":r["MarkerID"],"source_genome_build_inferred_from_rs671_coordinate":LOOKUP[(ch,pos)],
               "source_chr":ch,"source_pos":pos,
               "original_allele1":a,"original_effect_allele2":b,
               "effect_allele_harmonized":"A",
               "KCPS2_alcohol_beta_A":sign*beta,"KCPS2_alcohol_se":se,
               "KCPS2_ALCO_AMOUNT_p":p,
               "KCPS2_p_numeric_underflow":int(p==0),
               "KCPS2_neglogP_reported":r["neglogP"],
               "KCPS2_af_A":af if sign==1 else 1-af,
               "KCPS2_imputationInfo":info,"KCPS2_sample_size_at_variant":n,
               "KCPS2_effect_scale":"INVERSE_RANK_NORMAL_TRANSFORM",
               "source_original_Z":r["Z"],
               "source_provenance":"ONE_SELECTED_ORIGINAL_ZIP_MEMBER_CRC_AND_SHA256_PASS",
               "all_archive_MD5_verification":"NOT_DONE",
               "causal_stroke_effect":"NOT_ESTABLISHED"})
    if count<minimum_rows:raise ValueError("KCPS2 source has suspiciously few rows")
    if len(selected)!=1:
        raise RuntimeError(f"Exact rs671 source rows expected 1, found {len(selected)} across {count}; check genome build")
    x=selected[0]
    # Source 2025 paper independently states rs671-A alcohol beta around -0.59,
    # but tolerate version differences while blocking a falsely reversed sign.
    if x["KCPS2_alcohol_beta_A"]>=0:
        raise RuntimeError("KCPS2 source rs671-A effect conflicts published sign")
    out.mkdir(parents=True,exist_ok=True)
    result=out/"G0022_RS671_KCPS2_KOREAN_ALCOHOL_AMOUNT.tsv"
    with result.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(x),delimiter="\t");w.writeheader();w.writerow(x)
    summary={"source":"Jee_et_al_2025_Nature_Communications_KCPS2",
      "official_source_url":"https://zenodo.org/records/15132424",
      "published_cohort_size":153950,"full_KCPS2_archive_MD5_verified":False,
      "selected_ZIP_member_CRC32_verified":True,"selected_member_SHA256_verified":True,
      "selected_member_SAIGE_rows":count,
      "rs671_source_build_at_position":x["source_genome_build_inferred_from_rs671_coordinate"],
      "rs671_effect_allele":"A",
      "rs671_ALCO_AMOUNT_beta_A":x["KCPS2_alcohol_beta_A"],
      "rs671_ALCO_AMOUNT_se":x["KCPS2_alcohol_se"],
      "rs671_ALCO_AMOUNT_p":x["KCPS2_ALCO_AMOUNT_p"],
      "rs671_ALCO_AMOUNT_N":x["KCPS2_sample_size_at_variant"],
      "rs671_ALCO_AMOUNT_A_frequency":x["KCPS2_af_A"],
      "independent_of_Japanese_Koyanagi_BY_STUDY_POPULATION":True,
      "independent_EAS_AIS_STROKE_REPLICATION":False,
      "alcohol_measurement_scale_directly_comparable_to_Japanese":False,
      "causal_mediation_estimated":False,
      "status":"ORIGINAL_KOREAN_SOURCE_VARIANT_RETRIEVED_NONCAUSAL_REPLICATION"}
    (out/"G0022_RS671_KCPS2_KOREAN_SOURCE_AUDIT.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    return summary
if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--source",type=Path,default=DATA)
    a.add_argument("--out",type=Path,default=RESULT)
    a.add_argument("--min-source-rows",type=int,default=100000)
    v=a.parse_args();extract(v.source,v.out,v.min_source_rows)
