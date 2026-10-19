#!/usr/bin/env python3
"""Read-only, symmetric eigenvalue/PSD audit of preexisting CKD 1KG EUR LD.

Computes minimum eigenvalue with scipy.sparse.linalg.eigsh, checking full
spectrum only if requested. A negative value beyond numeric tolerance can
identify an LD numerical issue but cannot establish SuSiE nonconvergence cause.
"""
import argparse,csv,json,hashlib,math
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
GENES=("SDCCAG8","INHBC","UMOD")
def check(root,g):
    m=root/(g+".ld.meta.tsv");b=root/(g+".ld.bin")
    if not (m.is_file() and b.is_file()):return {"gene":g,"status":"MISSING"}
    with m.open(newline="") as f:rows=list(csv.DictReader(f,delimiter="\t"))
    n=len(rows)
    if n<2 or b.stat().st_size!=4*n*n:raise ValueError(g+" malformed LD dimension")
    raw=np.memmap(b,mode="r",dtype="<f4",shape=(n,n))
    if not np.all(np.isfinite(raw)):raise ValueError(g+" nonfinite LD")
    if np.max(np.abs(raw-raw.T))>1e-4:raise ValueError(g+" asymmetric LD")
    signs=np.array([float(r["effect_vs_ld_major_sign"]) for r in rows],dtype=np.float64)
    if not np.all(np.isin(signs,[-1.,1.])):raise ValueError(g+" invalid signs")
    # Congruence by diagonal +/-1 preserves the entire eigenspectrum.
    A=np.asarray(raw,dtype=np.float64)
    # Direct symmetric LAPACK solver avoids ARPACK stagnation near repeated zero eigenvalues.
    # Only the smallest eigenpair is requested; this is an exact dense subset solve.
    eigen,vec=eigh(A,subset_by_index=(0,0),driver="evr",check_finite=False,
                   overwrite_a=False)
    eigenvalue=float(eigen[0]);v=vec[:,0]
    residual=float(np.linalg.norm(A@v-eigenvalue*v))
    diagerr=float(np.max(np.abs(np.diag(A)-1.0)))
    threshold=max(1e-5, 10.0*residual)
    result={"gene":g,"variants":n,"min_eigenvalue_estimate":eigenvalue,"eigen_residual_norm":residual,
        "numeric_psd_tolerance":threshold,"max_diagonal_error":diagerr,
        "sign_congruence_preserves_spectrum":True,
        "ld_sha256":hashlib.sha256(b.read_bytes()).hexdigest(),
        "reference_ld_study_matched":False,"susie_nonconvergence_cause_determined":False}
    result["status"]=("PSD_NUMERIC_PASS" if eigenvalue>=-threshold else "PSD_NUMERIC_FAIL")
    return result
def main():
    p=argparse.ArgumentParser();p.add_argument("--lddir",required=True,type=Path);args=p.parse_args()
    entries=[check(args.lddir,g) for g in GENES]
    print(json.dumps({"kind":"REFERENCE_LD_EIGENVALUE_READ_ONLY","results":entries},indent=2))
    return 0 if all(x["status"]=="PSD_NUMERIC_PASS" for x in entries) else 2
if __name__=="__main__":raise SystemExit(main())
