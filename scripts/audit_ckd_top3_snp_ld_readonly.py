#!/usr/bin/env python3
"""Read-only Stage2B-to-Stage2C SNP/allele/LD metadata consistency audit.
Does not imply molecular-donor LD, conditional validity or causal identification.
"""
import csv,gzip,json,struct,argparse,math,hashlib
from pathlib import Path
GENES=("SDCCAG8","INHBC","UMOD")
def rows(p):
    op=gzip.open if p.suffix==".gz" else open
    with op(p,"rt",encoding="utf-8-sig",newline="") as f:
        r=csv.DictReader(f,delimiter="\t")
        if not r.fieldnames:raise ValueError("missing table header: "+str(p))
        result=list(r)
    if any(None in x for x in result):raise ValueError("malformed TSV "+str(p))
    return result
def check(gene,stage2b,stage2c,lddir):
    a=stage2b/(gene+".tsv.gz")
    b=stage2c/(gene+".tsv.gz")
    m=lddir/(gene+".ld.meta.tsv")
    v=lddir/(gene+".ld.vars")
    binary=lddir/(gene+".ld.bin")
    files=(a,b,m,v,binary)
    status={"gene":gene,"files":{p.name:{"present":p.is_file(),"bytes":p.stat().st_size if p.is_file() else None} for p in files},"status":"BLOCKED_MISSING_FILE"}
    if not all(x.is_file() for x in files):return status
    A,B,M=rows(a),rows(b),rows(m)
    required={"snp","pos37","allele1_pqtl","allele0_pqtl"}
    if not A or not B or not M or not required.issubset(A[0]):raise ValueError("source input schema invalid")
    for key in ("snp","pos37","beta_pqtl","beta_outcome_aligned"):
        if key not in B[0]:raise ValueError("SuSiE input missing "+key)
    for key in ("snp","ref_id","pos37","effect_allele","other_allele","effect_vs_ld_major_sign"):
        if key not in M[0]:raise ValueError("LD metadata missing "+key)
    def unique(d,label):
        ids=[x["snp"] for x in d]
        if len(set(ids))!=len(ids):raise ValueError(label+" duplicate SNP")
        return {x["snp"]:x for x in d}
    by_a=unique(A,"Stage2B");by_b=unique(B,"Stage2C");by_m=unique(M,"LD")
    order=[x.strip() for x in v.read_text().splitlines() if x.strip()]
    if len(order)!=len(set(order)) or order!=[r["ref_id"] for r in M]:
        raise ValueError("LD variant ordering mismatch")
    if binary.stat().st_size!=4*len(M)**2:raise ValueError("LD binary dimensions invalid")
    ld_not_susie=len(set(by_m)-set(by_b))
    susie_not_ld=len(set(by_b)-set(by_m))
    signs={x["effect_vs_ld_major_sign"] for x in M}
    if not signs.issubset({"1","-1","1.0","-1.0"}):raise ValueError("invalid signed LD orientation")
    mismatch=[]
    for x in M:
        k=x["snp"];a0=by_a.get(k);b0=by_b.get(k)
        if a0 is None:
            mismatch.append("ld_snp_absent_in_abf");continue
        if b0 is None:
            mismatch.append("ld_snp_absent_in_susie_input");continue
        if x["pos37"]!=a0["pos37"] or x["pos37"]!=b0["pos37"]:
            mismatch.append("position_mismatch")
        if x["effect_allele"].upper()!=a0["allele1_pqtl"].upper() or x["other_allele"].upper()!=a0["allele0_pqtl"].upper():
            mismatch.append("allele_mismatch")
    # Random access only: diagonal entries and adjacent samples, not eigen decomposition.
    n=len(M); bad_diag=0;nonfinite=0
    with binary.open("rb") as f:
        for i in range(n):
            f.seek(4*(i*n+i));z=struct.unpack("<f",f.read(4))[0]
            bad_diag+=abs(z-1)>0.001
            nonfinite+=not math.isfinite(z)
    issues=[]
    # Stage2C input intentionally includes variants unavailable in 1KG LD.\n    # Its R runner restricts to meta SNPs using match(meta$snp, dat$snp).\n    if ld_not_susie:issues.append("LD_SNP_ABSENT_IN_SUSIE_INPUT")
    if mismatch:issues.append("ALLELE_OR_POSITION_MISMATCH")
    if bad_diag or nonfinite:issues.append("LD_DIAGONAL_INVALID")
    status.update({"status":"METADATA_PASS_FILTERED_REFERENCE_LD_ONLY" if not issues else "QC_FAIL",
      "stage2b_rows":len(A),"stage2c_rows":len(B),"ld_rows":len(M),
      "ld_not_susie":ld_not_susie,"susie_not_ld":susie_not_ld,
      "abf_ld_shared_snps":len(set(by_a)&set(by_m)),
      "abf_only_snps":len(set(by_a)-set(by_m)),
      "allele_position_errors":len(mismatch),"ld_bad_diagonal":bad_diag,
      "ld_nonfinite_diagonal":nonfinite,"issues":issues,
      "ld_reference_is_molecular_study_matched":False,
      "full_LD_PSD_and_summary_effect_orientation_verified":False,
      "sha256_meta":hashlib.sha256(m.read_bytes()).hexdigest()})
    return status
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--stage2b",type=Path,required=True);p.add_argument("--stage2c",type=Path,required=True)
    p.add_argument("--lddir",type=Path,required=True)
    a=p.parse_args()
    data=[check(g,a.stage2b,a.stage2c,a.lddir) for g in GENES]
    print(json.dumps({"audit":"CKD_STAGE2B_STAGE2C_REFERENCE_LD_READONLY","results":data},indent=2))
    return 0 if all(x["status"]=="METADATA_PASS_FILTERED_REFERENCE_LD_ONLY" for x in data) else 2
if __name__=="__main__":raise SystemExit(main())
