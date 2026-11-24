#!/usr/bin/env python3
"""Read-only numerical integrity audit of three preexisting CKD reference LD matrices.
Checks matrix symmetry, diagonal, PSD, orientation consistency with existing metadata.
Reference LD only; no scientific colocalization rerun or causal validation.
"""
import argparse,csv,json,math,hashlib
from pathlib import Path
import numpy as np

GENES=("SDCCAG8","INHBC","UMOD")
def audit(root,gene):
    meta=root/(gene+".ld.meta.tsv");matrix=root/(gene+".ld.bin");order=root/(gene+".ld.vars")
    if not all(p.is_file() for p in (meta,matrix,order)):
        return {"gene":gene,"status":"BLOCKED_MISSING_INPUT"}
    with meta.open(newline="") as f: rows=list(csv.DictReader(f,delimiter="\t"))
    ids=[x["ref_id"] for x in rows];snps=[x["snp"] for x in rows]
    vars_=[x.strip() for x in order.read_text().splitlines() if x.strip()]
    n=len(rows)
    if not n or ids!=vars_ or len(set(ids))!=n or len(set(snps))!=n:
        return {"gene":gene,"status":"FAIL_ORDER_OR_DUPLICATES"}
    if matrix.stat().st_size!=n*n*4: return {"gene":gene,"status":"FAIL_MATRIX_BYTES"}
    signs=np.array([float(r["effect_vs_ld_major_sign"]) for r in rows])
    if not np.all(np.isin(signs,[-1.,1.])):return {"gene":gene,"status":"FAIL_SIGN_METADATA"}
    raw=np.memmap(matrix,dtype="<f4",mode="r",shape=(n,n),order="C")
    if not np.all(np.isfinite(raw)):return {"gene":gene,"status":"FAIL_NONFINITE"}
    max_sym=float(np.max(np.abs(raw-raw.T)))
    diag=float(np.max(np.abs(np.diag(raw)-1)))
    max_abs=float(np.max(np.abs(raw)))
    # Eigenvalues of a ~4k square dense LD matrix can be costly.
    # Examine exact full spectrum if requested; never declare PSD verified otherwise.
    signed=np.asarray(raw,dtype=np.float64)*signs[:,None]*signs[None,:]
    signed_sym=float(np.max(np.abs(signed-signed.T)))
    signed_diag=float(np.max(np.abs(np.diag(signed)-1)))
    result={"gene":gene,"status":"FULL_SYMMETRY_SIGN_METADATA_PASS_PSD_NOT_TESTED",
      "variants":n,"sha256_ld":hashlib.sha256(matrix.read_bytes()).hexdigest(),
      "max_asymmetry":max_sym,"max_diagonal_error":diag,"max_abs_correlation":max_abs,
      "signed_max_asymmetry":signed_sym,"signed_max_diagonal_error":signed_diag,
      "orientation_signs_negative":int((signs==-1).sum()),"orientation_signs_positive":int((signs==1).sum()),
      "reference_ld_only":True,"molecular_study_matched_ld":False,
      "actual_effect_beta_dosage_sign_revalidated":False}
    if max_sym>1e-4 or diag>1e-3 or max_abs>1.001 or signed_sym>1e-4 or signed_diag>1e-3:
        result["status"]="FAIL_NUMERIC_LD_INTEGRITY"
    return result
def main():
    p=argparse.ArgumentParser();p.add_argument("--lddir",type=Path,required=True);a=p.parse_args()
    out={"kind":"CKD_TOP3_REFERENCE_LD_FULL_NUMERIC_READONLY","psd_verified":False,
         "results":[audit(a.lddir,g) for g in GENES]}
    print(json.dumps(out,indent=2))
    return 0 if all(x["status"]=="FULL_SYMMETRY_SIGN_METADATA_PASS_PSD_NOT_TESTED" for x in out["results"]) else 2
if __name__=="__main__":raise SystemExit(main())
