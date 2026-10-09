#!/usr/bin/env python3
"""Aggregate broad IS GWS reference readiness without claiming valid fine-mapping."""
import argparse,csv,json
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
REF=Path("/srv/is-analysis/data/is/ld_reference/broad_v1/1kg_eas")
ANCHORS={"IS_XDATA_G0007":"L001","IS_XDATA_G0015":"L002",
         "IS_XDATA_G0022":"L003","IS_XDATA_G0023":"L004"}
def build(root,ref):
    with (root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv").open() as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    output=[]
    for row in rows:
        g=row["group_id"]
        if g in ANCHORS:
            prefix=Path("/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/pgen_gt")/f"BBJ_IS_{ANCHORS[g]}.1KG_EAS.GRCh37"
            stage="LEGACY_REGIONAL_REFERENCE"
            allele={"status":"POSITION_VERIFIED_EARLIER"}
        else:
            prefix=ref/(g+".stable")
            status=ref/(g+".stable.REFERENCE_QC.json")
            allele=json.loads(status.read_text()) if status.is_file() else {}
            stage="NEW_GWS_REFERENCE" if allele.get("status")=="EXTRACTED_AND_CONVERTED" else "MISSING"
        complete=all(Path(str(prefix)+x).is_file() for x in (".pgen",".pvar",".psam"))
        in_scope=row["has_gws"]=="1"
        result={
            "group_id":g,"chrom":row["chr"],"has_gws":int(in_scope),
            "phenotypes":row["phenotypes"],"lead_variant":row["lead_variant"],
            "reference_source":stage if complete else "MISSING",
            "pgen_triplet_complete":int(complete),
            "lead_allele_status":allele.get("lead_allele_status",allele.get("status","NOT_CHECKED")),
            "ld_proxy_status":"READY_TO_COMPUTE_REF_R2" if complete else "BLOCKED_NO_LOCAL_REF",
            "independent_locus_status":"NOT_TESTED",
            "fine_mapping_status":"BLOCKED_ALLELE_HARMONIZATION_AND_FULL_LD_QC"}
        output.append(result)
    with (root/"IS_BROAD_REFERENCE_READINESS.tsv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(output[0]),delimiter="\t")
        writer.writeheader();writer.writerows(output)
    summary={"total_groups":len(rows),"gws_groups":sum(x["has_gws"] for x in output),
        "all_groups_with_local_reference":sum(x["pgen_triplet_complete"] for x in output),
        "gws_with_local_reference":sum(x["has_gws"] and x["pgen_triplet_complete"] for x in output),
        "gws_missing_local_reference":sum(x["has_gws"] and not x["pgen_triplet_complete"] for x in output),
        "all_groups_without_local_reference":sum(not x["pgen_triplet_complete"] for x in output),
        "fine_mapping_passed":0,
        "note":"Regional PGEN is not sufficient for independent signals, fine-mapping, or allele-harmonized MR/coloc"}
    (root/"IS_BROAD_REFERENCE_READINESS_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    return summary
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--ref",type=Path,default=REF)
    args=p.parse_args()
    print(json.dumps(build(args.root,args.ref),indent=2))
