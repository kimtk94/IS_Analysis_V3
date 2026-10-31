#!/usr/bin/env python3
"""Single-transfer, single-decompression multi-label raw UMI pseudobulk.

Fixed metadata cohort is BeamA>=0.75 AND scpred TM_Fibro shared donor set.
Two labels may overlap in the same cell. We explicitly count into *both* sums
in the same sparse read, preserving within-label paired biological donor units.
Never export donor IDs, individual raw counts, or public GEO genotype data.
"""
import argparse,csv,gzip,hashlib,json,re,tempfile
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from ocular_p1_integrity import BASE,file_url,bar,lines
from ocular_p4_fullcohort import pool_name
from ocular_p4_pseudobulk_pilot import GENES,summary_numbers
from ocular_gate_b_050_p4_retry import resilient_geo_download as download

LABELS=("BeamA_075","TM_Fibro")
ARMS=("no-dex","dex")
def get_meta(fp):
    maps={k:defaultdict(dict) for k in LABELS}
    totals={k:defaultdict(Counter) for k in LABELS}
    source_libraries=Counter()
    by_pool=defaultdict(set)
    rows=0
    with gzip.open(fp,"rt",encoding="utf8") as h:
        rd=csv.DictReader(h,delimiter="\t")
        req={"library_id","DEX Treatment","barcode","donor_id_demuxlet","droplet_type",
             "celltypist_chen2025_fine","celltypist_chen2025_fine_score","scpred_vanzyl2022"}
        if not req.issubset(rd.fieldnames or []):raise ValueError("SOURCE_METADATA_SCHEMA_CHANGED")
        for r in rd:
            rows+=1
            lib=r["library_id"].strip();arm=r["DEX Treatment"].strip().lower()
            if arm not in ARMS or pool_name(lib) is None:continue
            source_libraries[(lib,arm)]+=1
            if r["droplet_type"]!="SNG":continue
            donor=r["donor_id_demuxlet"].strip()
            if not donor or donor.lower() in ("na","doublet","unassigned"):continue
            bc=bar(r["barcode"])
            try:score=float(r["celltypist_chen2025_fine_score"])
            except (TypeError,ValueError):score=float("nan")
            flags={"BeamA_075":r["celltypist_chen2025_fine"]=="BeamA" and score>=.75,
                   "TM_Fibro":r["scpred_vanzyl2022"]=="TM_Fibro"}
            if not any(flags.values()):continue
            by_pool[donor].add(pool_name(lib))
            for label,yes in flags.items():
                if yes:
                    key=(lib,arm)
                    if bc in maps[label][key]:raise ValueError("DUPLICATE_SELECTED_BARCODE")
                    maps[label][key][bc]=donor
                    totals[label][(donor,arm)]+=1
    if rows!=169398 or len(source_libraries)!=18 or any(len(s)!=1 for s in by_pool.values()):
        raise ValueError("GEO_METADATA_ROWS_OR_POOL_DRIFT")
    def eligible(label):
        return {d for d,a in totals[label] if
                totals[label][d,"dex"]>=30 and totals[label][d,"no-dex"]>=30}
    beam=eligible("BeamA_075");fibro=eligible("TM_Fibro")
    common=sorted(beam&fibro)
    if len(beam)!=57 or len(fibro)!=57 or len(common)!=55:
        raise ValueError("EXPECTED_SHARED_55_METADATA_DONORS_CHANGED")
    return maps,common,by_pool,totals

def stream_both(path,source_ncol,group_indices):
    """One pass over coordinates; two bincount aggregations per chunk."""
    n=len(group_indices["BeamA_075"])
    if any(len(ix)!=n for ix in group_indices.values()):raise ValueError("BARCODE_IDX_LENGTH_MISMATCH")
    n_donors=max((int(ix.max()) for ix in group_indices.values()),default=-1)+1
    # The joined donor index vectors may not mention all donors in one specific library.
    # Pass a caller-owned n_groups instead of inferring it.
    raise AssertionError("DO_NOT_USE_INFERRED_DIMENSION")
def matrix_onepass(path,source_ncol,mappings,n_genes,n_donors):
    outputs={name:np.zeros((n_genes,n_donors),dtype=np.int64) for name in LABELS}
    selected_umi={name:0 for name in LABELS}
    selected_nnz={name:0 for name in LABELS}
    n=0
    with gzip.open(path,"rb") as f:
        if f.readline().strip()!=b"%%MatrixMarket matrix coordinate integer general":
            raise ValueError("BAD_MATRIX_MARKET_HEADER")
        while True:
            line=f.readline()
            if not line:raise ValueError("NO_MATRIX_DIMENSIONS")
            if not line.startswith(b"%"):break
        rows,cols,nnz=map(int,line.split())
        if rows!=n_genes or cols!=source_ncol or nnz<1:
            raise ValueError("MATRIX_DIMENSION_DRIFT")
        carry=b""
        while True:
            chunk=f.read(4*1024*1024)
            if not chunk:
                if not carry:break
                chunk=carry+b"\n"
                carry=b""
            else:chunk=carry+chunk
            pos=chunk.rfind(b"\n")
            if pos<0:
                if len(chunk)>1024:raise ValueError("UNEXPECTED_LONG_TRIPLE")
                carry=chunk;continue
            data=chunk[:pos+1];carry=chunk[pos+1:]
            triples=np.fromstring(data.decode("ascii"),sep=" ",dtype=np.int64)
            if not len(triples):continue
            if len(triples)%3:raise ValueError("SPLIT_COORDINATE_TRIPLE")
            triples=triples.reshape(-1,3)
            ri=triples[:,0]-1;ci=triples[:,1]-1;v=triples[:,2]
            if (np.any(ri<0)|np.any(ri>=rows)|np.any(ci<0)|np.any(ci>=cols)|np.any(v<=0)):
                raise ValueError("INVALID_SPARSE_COORDINATE")
            n+=len(triples)
            for label,arr in mappings.items():
                inds=arr[ci]
                keep=inds>=0
                if not np.any(keep):continue
                flat=ri[keep]*n_donors+inds[keep]
                weights=np.bincount(flat,weights=v[keep].astype(np.float64),minlength=n_genes*n_donors)
                if weights.max()>=2**53:raise ValueError("COUNTS_FLOAT_PRECISION_LOSS")
                outputs[label]+=np.rint(weights).astype(np.int64).reshape(n_genes,n_donors)
                selected_umi[label]+=int(v[keep].sum())
                selected_nnz[label]+=int(keep.sum())
        if n!=nnz:raise ValueError("MATRIX_NNZ_NOT_MATCHED")
    for label in LABELS:
        if int(outputs[label].sum())!=selected_umi[label]:
            raise ValueError("PSEUDOBULK_UMI_CONSERVATION_FAILED_"+label)
    return outputs,{"nnz":n,"gzip_crc":"PASS",
                   "umi_by_annotation":selected_umi,"selected_nnz_by_annotation":selected_nnz}
def run(args):
    with tempfile.TemporaryDirectory(prefix="ocular-55-") as td:
        root=Path(td)
        meta_path=root/"metadata.tsv.gz"
        meta_src=download(BASE+"series/GSE313nnn/GSE313736/suppl/GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",
                          meta_path,12*1024*1024)
        maps,donors,pools,cellcounts=get_meta(meta_path)
        n=len(donors)
        donor_ix={d:i for i,d in enumerate(donors)}
        umi={lab:{arm:np.zeros((GENES,n),dtype=np.int64) for arm in ARMS} for lab in LABELS}
        cells={lab:{arm:np.zeros(n,dtype=np.int64) for arm in ARMS} for lab in LABELS}
        evidence={"status":"HOLD_RAW_55DUAL","planned_common_donors":55,
            "original_GEO_metadata_sha256":meta_src["sha256"],"per_library":[],
            "no_donor_identity_or_count_matrix_uploaded":True,
            "DEX_library_batch_completely_aliased":True,
            "post_treatment_cell_label_selection":True}
        gene_digest=None;features=None
        for i in range(9):
            pool=chr(65+i)
            for arm,gsm,prefix,lib in (
                    ("no-dex",f"GSM{9373276+i}",f"pool_{pool}",f"HTMC_pool-{pool}_GRCh38"),
                    ("dex",f"GSM{9373285+i}",f"pool_dex{pool}",f"HTMC_pool-dex{pool}_GRCh38")):
                barfp=root/"barcodes.tsv.gz"; featfp=root/"features.tsv.gz";matfp=root/"matrix.mtx.gz"
                x=download(file_url(gsm,prefix,"barcodes.tsv.gz"),barfp,12*1024*1024)
                y=download(file_url(gsm,prefix,"features.tsv.gz"),featfp,12*1024*1024)
                bs=lines(barfp);fs=lines(featfp)
                if len(fs)!=GENES or len(bs)!=len(set(map(bar,bs))):
                    raise ValueError("FEATURE_BARCODE_DIMENSION_ERROR")
                digest=hashlib.sha256(("\n".join(z.split("\t")[0] for z in fs)+"\n").encode()).hexdigest()
                if gene_digest is None:gene_digest=digest;features=fs
                elif gene_digest!=digest:raise ValueError("SOURCE_FEATURE_ORDER_MISMATCH")
                index={}
                rec={"pool":pool,"arm":arm,"GSM":gsm,"feature_hash":y["sha256"],
                     "barcode_hash":x["sha256"]}
                for lab in LABELS:
                    selected=maps[lab].get((lib,arm),{})
                    arr=np.asarray([donor_ix.get(selected.get(bar(z)),-1) for z in bs],dtype=np.int32)
                    index[lab]=arr
                    included=arr>=0
                    joined=np.bincount(arr[included],minlength=n)
                    cells[lab][arm]+=joined
                    expected=sum(1 for d in selected.values() if d in donor_ix)
                    recall=int(included.sum())/max(1,expected)
                    rec[lab+"_matched_barcodes"]=int(included.sum())
                    rec[lab+"_barcode_recall"]=round(recall,5)
                    if recall<.85:raise ValueError("BARCODE_JOIN_RECALL_TOO_LOW "+lab+" "+pool)
                tr=download(file_url(gsm,prefix,"matrix.mtx.gz"),matfp,620*1024*1024)
                out,info=matrix_onepass(matfp,len(bs),index,GENES,n)
                for lab in LABELS:umi[lab][arm]+=out[lab]
                rec["matrix_sha256"]=tr["sha256"]
                rec["matrix_nnz"]=info["nnz"]
                rec["onepass_decompression"]="PASS"
                evidence["per_library"].append(rec)
                matfp.unlink(missing_ok=True)
                barfp.unlink(missing_ok=True)
                featfp.unlink(missing_ok=True)
                print("DUAL_PASS",pool,arm,info["umi_by_annotation"],flush=True)
        mask=np.ones(n,dtype=bool)
        for lab in LABELS:
            for arm in ARMS:mask&=(cells[lab][arm]>=30)
        retained=int(mask.sum())
        evidence["matrix_joined_common_donors_ge30_every_arm_and_annotation"]=retained
        evidence["original_fixed_shared_55"] = 55
        evidence["donors_failed_after_barcode_join"]=55-retained
        evidence["per_arm_annotation"]={
            lab:{arm:{"matched_cells":int(cells[lab][arm][mask].sum()),
                      "umi":int(umi[lab][arm][:,mask].sum())} for arm in ARMS}
            for lab in LABELS}
        # Fail closed instead of silently changing the 55-donor comparison.
        if retained!=55:
            evidence["status"]="HOLD_55_DONOR_JOIN_LOSS_NOT_IDENTICAL_COHORT"
            return evidence
        destination=Path(args.private)
        destination.mkdir(parents=True,exist_ok=True)
        out_features=destination/"gene_features.tsv"
        with out_features.open("w",newline="") as f:
            wr=csv.writer(f,delimiter="\t")
            wr.writerow(["gene_id","gene_symbol"])
            for entry in features:
                a=entry.split("\t")
                wr.writerow([a[0],a[1] if len(a)>1 else a[0]])
        for lab in LABELS:
            directory=destination/lab;directory.mkdir(parents=True,exist_ok=True)
            labels=[]
            with (directory/"samples.tsv").open("w",newline="") as f:
                wr=csv.writer(f,delimiter="\t")
                wr.writerow(["sample_id","pair_id","treatment","pool","selected_cells"])
                for d,idx in sorted(donor_ix.items(),key=lambda z:z[1]):
                    for arm,tr in (("no-dex","Vehicle"),("dex","DEX")):
                        sid=f"p{idx+1:03d}_{tr}"
                        labels.append(sid)
                        wr.writerow([sid,f"p{idx+1:03d}",tr,next(iter(pools[d])),int(cells[lab][arm][idx])])
            matrix=np.empty((GENES,2*n),dtype=np.int64)
            matrix[:,::2]=umi[lab]["no-dex"]
            matrix[:,1::2]=umi[lab]["dex"]
            with gzip.open(directory/"paired_raw_counts.tsv.gz","wt",compresslevel=3,newline="") as f:
                wr=csv.writer(f,delimiter="\t")
                wr.writerow(["gene_id"]+labels)
                for gi,x in zip((z.split("\t")[0] for z in features),matrix):
                    wr.writerow([gi]+[int(v) for v in x])
            # R expects feature file in each directory.
            import shutil
            shutil.copyfile(out_features,directory/"gene_features.tsv")
        evidence["status"]="PASS_55_DUAL_ANNOTATION_RAW_PSEUDOBULK"
        evidence["gene_order_sha256"]=gene_digest
        evidence["features"]=GENES
        evidence["two_labels_in_one_matrix_parse"]=True
        evidence["two_annotation_samples_per_label"]=110
        return evidence

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--private",required=True)
    p.add_argument("--out",required=True)
    a=p.parse_args();out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    report={"status":"HOLD"}
    try:report=run(a)
    except Exception as e:report={"status":"HOLD","error":type(e).__name__+": "+str(e)[:250],
                                  "private_data_cannot_leave_runner":True}
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("GATE_B_55",report["status"],"retained",report.get("matrix_joined_common_donors_ge30_every_arm_and_annotation"),
          "error",report.get("error"),flush=True)
    if report["status"]!="PASS_55_DUAL_ANNOTATION_RAW_PSEUDOBULK":raise SystemExit(2)
if __name__=="__main__":main()
