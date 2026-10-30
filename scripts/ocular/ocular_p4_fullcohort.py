#!/usr/bin/env python3
"""P4 full-cohort, in-memory donor-paired BeamA pseudobulk QC (GSE313736).

All 18 original public 10x matrices are streamed sequentially; sparse raw UMI
counts are summed by biological donor and treatment without cell pseudoreplication.
Artifacts contain aggregate QC and file hashes only; NEVER donor identifiers,
per-donor expression matrices, genotypes, or clinical phenotypes.
"""
import argparse
import csv
import gzip
import hashlib
import json
import re
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from ocular_p1_integrity import BASE, bar, download, file_url, lines
from ocular_p4_pseudobulk_pilot import ANNOTATION, SCORE, GENES, matrix_stream, summary_numbers


def pool_name(lib):
    m = re.search(r"pool[-_](?:dex)?([A-I])(?:_|$)", lib)
    return m.group(1) if m else None


def collect_metadata(path, cutoff):
    # Only high-confidence SNG BeamA records selected. Donor identifiers stay in RAM.
    data = defaultdict(dict)
    n_seen = 0
    all_libraries = Counter()
    donor_pools = defaultdict(set)
    with gzip.open(path, "rt", encoding="utf8") as h:
        reader = csv.DictReader(h, delimiter="\t")
        required = {"library_id", "DEX Treatment", "barcode",
                    "donor_id_demuxlet", "droplet_type", ANNOTATION, SCORE}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("MISSING_COLUMNS " + str(sorted(required-set(reader.fieldnames or []))))
        for r in reader:
            n_seen += 1
            lib, arm, donor = (r["library_id"].strip(),
                               r["DEX Treatment"].strip().lower(),
                               r["donor_id_demuxlet"].strip())
            if arm not in ("dex", "no-dex") or pool_name(lib) is None:
                continue
            all_libraries[(lib, arm)] += 1
            if r["droplet_type"] != "SNG" or r[ANNOTATION] != "BeamA":
                continue
            try:
                confidence = float(r[SCORE])
            except (TypeError, ValueError):
                continue
            if confidence < cutoff or not donor or donor.lower() in ("na", "doublet", "unassigned"):
                continue
            bc = bar(r["barcode"])
            if bc in data[(lib,arm)]:
                raise ValueError("DUPLICATE_METADATA_BARCODE "+lib)
            data[(lib,arm)][bc]=donor
            donor_pools[donor].add(pool_name(lib))
    if n_seen!=169398 or len(all_libraries)!=18:
        raise ValueError("UNEXPECTED_METADATA_COUNT_OR_LIBRARY_COUNT")
    if any(len(x)!=1 for x in donor_pools.values()):
        raise ValueError("DONOR_PRESENT_IN_MULTIPLE_POOLS")
    return data,all_libraries,donor_pools


def evaluate(output,confidence=0.75,min_cells=30):
    result={"dataset":"GSE313736","study_scope":"in-vitro human trabecular meshwork DEX response, not clinical IOP",
            "stage":"P4_full_cohort_raw_UMI_pseudobulk",
            "run_utc":datetime.now(timezone.utc).isoformat(),
            "score_cutoff":confidence,"min_cells_per_donor_per_condition":min_cells,
            "status":"HOLD","per_library":[],"no_donor_ids_saved":True,
            "individual_counts_written":False}
    with tempfile.TemporaryDirectory(prefix="ocular-p4-full-") as td:
        root=Path(td)
        meta_path=root/"cell_metadata.gz"
        result["metadata_transfer"]=download(
            BASE+"series/GSE313nnn/GSE313736/suppl/GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",
            meta_path,12*1024*1024)
        meta,all_libs,donor_pools=collect_metadata(meta_path,confidence)
        per_donor_arm=defaultdict(Counter)
        for (lib,arm), cells in meta.items():
            per_donor_arm[arm].update(cells.values())
        matched=sorted(d for d in (per_donor_arm["no-dex"].keys() & per_donor_arm["dex"].keys())
                       if per_donor_arm["no-dex"][d]>=min_cells and per_donor_arm["dex"][d]>=min_cells)
        if len(matched)<50:
            raise ValueError("FEWER_THAN_50_PAIRED_BEAMA_DONORS")
        group_ix={d:i for i,d in enumerate(matched)}
        result["eligible_paired_donors_by_metadata"]=len(matched)
        result["selected_metadata_cells"]={
          arm:sum(per_donor_arm[arm][d] for d in matched)
          for arm in ("no-dex","dex")
        }
        total={arm:np.zeros((GENES,len(matched)),dtype=np.int64)
               for arm in ("no-dex","dex")}
        joined={arm:np.zeros(len(matched),dtype=np.int64) for arm in total}
        gene_order_digest=None
        for pool_idx in range(9):
            p=chr(ord("A")+pool_idx)
            for arm,condition,gsm,prefix,lib in [
              ("no-dex","vehicle",f"GSM{9373276+pool_idx}",f"pool_{p}",f"HTMC_pool-{p}_GRCh38"),
              ("dex","dex",f"GSM{9373285+pool_idx}",f"pool_dex{p}",f"HTMC_pool-dex{p}_GRCh38")]:
                rec={"pool":p,"condition":condition,"gsm":gsm,"status":"HOLD"}
                lib_data=meta.get((lib,arm))
                if lib_data is None or (lib,arm) not in all_libs:
                    raise ValueError("LIBRARY_KEY_NOT_FOUND "+lib+" "+arm)
                bpath=root/"barcodes.tsv.gz";fpath=root/"features.tsv.gz";mpath=root/"matrix.mtx.gz"
                rec["barcode_file"]=download(file_url(gsm,prefix,"barcodes.tsv.gz"),bpath,12*1024*1024)
                rec["feature_file"]=download(file_url(gsm,prefix,"features.tsv.gz"),fpath,12*1024*1024)
                features=lines(fpath)
                if len(features)!=GENES:
                    raise ValueError("FEATURE_LENGTH_CHANGED "+gsm)
                digest=hashlib.sha256(("\n".join(s.split("\t")[0] for s in features)+"\n").encode()).hexdigest()
                if gene_order_digest is None:
                    gene_order_digest=digest
                elif digest!=gene_order_digest:
                    raise ValueError("FEATURE_ORDER_DIFFERS "+gsm)
                barcodes=lines(bpath)
                if len(barcodes)!=len(set(map(bar,barcodes))):
                    raise ValueError("DUPLICATE_BARCODE "+gsm)
                ix=np.array([group_ix.get(lib_data.get(bar(v)),-1)
                             for v in barcodes],dtype=np.int32)
                good=ix>=0
                count=np.bincount(ix[good],minlength=len(matched))
                joined[arm]+=count
                expected=sum(1 for x in lib_data.values() if x in group_ix)
                recall=int(good.sum())/max(expected,1)
                rec["selected_barcodes"]=int(good.sum())
                rec["selected_barcode_recall"]=round(recall,5)
                rec["source_filtered_barcodes"]=len(barcodes)
                if recall<0.90:raise ValueError("BARCODE_RECALL_LT_90PCT "+gsm)
                if recall<0.95:
                    rec["barcode_recall_warning"]="BETWEEN_90_AND_95_PERCENT; monitor dropout and gene-coverage sensitivity"
                rec["matrix_transfer"]=download(file_url(gsm,prefix,"matrix.mtx.gz"),mpath,620*1024*1024)
                pooled,info=matrix_stream(mpath,len(barcodes),ix,GENES,len(matched))
                total[arm]+=pooled
                rec["matrix_info"]=info
                rec["status"]="PASS_PSEUDOBULK"
                result["per_library"].append(rec)
                mpath.unlink(missing_ok=True);bpath.unlink(missing_ok=True);fpath.unlink(missing_ok=True)
                print("POOL",p,arm,"UMI",info["selected_umis"],"nnz",info["nnz"],flush=True)
        eligible_mask=(joined["dex"]>=min_cells)&(joined["no-dex"]>=min_cells)
        n_kept=int(np.count_nonzero(eligible_mask))
        if n_kept<50:raise ValueError("POST_JOIN_PAIRED_DONORS_LT_50 "+str(n_kept))
        result["metadata_paired_donors_after_actual_barcode_join"]=n_kept
        result["genes"]=GENES
        result["gene_order_sha256"]=gene_order_digest
        result["per_arm"]={}
        filtered={arm:mat[:,eligible_mask] for arm,mat in total.items()}
        for arm in ("no-dex","dex"):
            mat=filtered[arm]
            result["per_arm"][arm]={
              "matched_cells":int(joined[arm][eligible_mask].sum()),
              "total_selected_umi":int(mat.sum()),
              "umi_per_biological_donor":summary_numbers(mat.sum(axis=0)),
              "genes_detected_per_biological_donor":summary_numbers(np.count_nonzero(mat,axis=0))}
        across=filtered["no-dex"]+filtered["dex"]
        result["gene_detected_in_any_donor"]=int(np.count_nonzero(across.sum(axis=1)))
        result["gene_detected_in_at_least_10_donors_each_arm"]=int(np.count_nonzero(
            (np.count_nonzero(filtered["no-dex"],axis=1)>=10)
            &(np.count_nonzero(filtered["dex"],axis=1)>=10)))
        result["status"]="GO_57DONOR_FULL_RAW_PSEUDOBULK_QC" if n_kept==57 else "GO_FULL_RAW_PSEUDOBULK_QC"
        result["restrictions"]=["No clinical steroid-induced IOP phenotype",
                                "No DE testing or independent replication has been performed",
                                "Individual donor count matrices must not be committed to Git"]
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",required=True)
    ap.add_argument("--score",type=float,default=0.75)
    ap.add_argument("--min-cells",type=int,default=30)
    args=ap.parse_args()
    evidence={"status":"HOLD","run_utc":datetime.now(timezone.utc).isoformat()}
    try:
        evidence=evaluate(args.output,args.score,args.min_cells)
    except Exception as e:
        evidence.update(status="HOLD_FULL_PSEUDOBULK",
                        error=type(e).__name__+": "+str(e)[:220])
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(evidence,indent=2,ensure_ascii=False)+"\n")
    print("P4_FULL_RESULT",evidence["status"],
          "paired",evidence.get("metadata_paired_donors_after_actual_barcode_join"),
          "error",evidence.get("error"),flush=True)
    if evidence["status"] not in ("GO_57DONOR_FULL_RAW_PSEUDOBULK_QC",
                                  "GO_FULL_RAW_PSEUDOBULK_QC"):
        raise SystemExit(2)


if __name__=="__main__":
    main()
