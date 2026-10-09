#!/usr/bin/env python3
"""Broad IS pre-finemapping inventory. No genotype downloads or invasive writes."""
import argparse,csv,json
from pathlib import Path
from collections import defaultdict

def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def pvar_header(line): return line.startswith("##") or line.startswith("#CHROM")
def main():
    a=argparse.ArgumentParser()
    a.add_argument("--root",type=Path,default=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1"))
    a.add_argument("--ref",type=Path,default=Path("/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/pgen_gt"))
    o=a.parse_args()
    groups=rows(o.root/"CROSS_DATASET_INTERVAL_GROUPS.tsv")
    positions=defaultdict(set); byref={}
    for p in sorted(o.ref.glob("*.pvar")):
        stem=p.with_suffix("")
        if not all(Path(str(stem)+x).exists() for x in (".pgen",".psam")):continue
        n=0
        with p.open() as f:
            for line in f:
                if pvar_header(line):continue
                v=line.rstrip("\n").split("\t")
                if len(v)<5:continue
                try:pos=int(v[1])
                except ValueError:continue
                chrom=v[0].removeprefix("chr")
                positions[(chrom,pos)].add(p.name)
                n+=1
        byref[p.name]=n
    output=[]
    for r in groups:
        key=r["lead_variant"].split(":")
        lead_pos=int(key[1])
        exact=sorted(positions.get((r["chr"],lead_pos),set()))
        # Variant-level REF/ALT matching is not asserted here: positional presence only.
        size=int(r["end"])-int(r["start"])+1
        output.append({"group_id":r["group_id"],"chr":r["chr"],"region_start":r["start"],
            "region_end":r["end"],"window_start":max(1,int(r["start"])-500000),
            "window_end":int(r["end"])+500000,"window_bp":size+1000000,
            "lead_variant":r["lead_variant"],"lead_p":r["lead_p"],
            "phenotypes":r["phenotypes"],"has_gws":r["has_gws"],
            "lead_position_in_existing_reference":int(bool(exact)),
            "lead_position_reference_files":";".join(exact),
            "ld_status":"LEAD_POSITION_FOUND_ALLELE_QC_REQUIRED" if exact else "NEEDS_NEW_ANCESTRY_MATCHED_REFERENCE",
            "gene_annotation_status":"UNRESOLVED","fine_mapping_status":"BLOCKED_UNTIL_LD_AND_ALLELE_QC"})
    fields=list(output[0])
    with (o.root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t");w.writeheader();w.writerows(output)
    report={"region_groups":len(output),"reference_pvar_files":byref,
        "lead_position_present":sum(x["lead_position_in_existing_reference"] for x in output),
        "lead_position_absent":sum(not x["lead_position_in_existing_reference"] for x in output),
        "source_warning":"lead position is not allele-verified; no claim of reference LD compatibility",
        "next_gate":"genome-wide EAS phased genotype reference or per-region VCF plus harmonization"}
    (o.root/"IS_EXPANDED_EXECUTION_STATUS.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
