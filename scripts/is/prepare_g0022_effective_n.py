#!/usr/bin/env python3
"""GIGASTROKE EAS case-control approximate effective sample sizes with provenance.

Study-wide counts NOT necessarily per-SNP effective N for meta-analyzed logistic GWAS.
"""
import argparse,csv,json
from pathlib import Path
import yaml

META=Path("/srv/is-analysis/data/is/reference/gigastroke/metadata/yaml")
OUT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/susie_g0022_pilot")
COUNTS=[
 {"accession":"GCST90104544","trait":"AS","cases":27413,"controls":237242,
  "support_url":"https://kp4cd.org/node/1138"},
 {"accession":"GCST90104545","trait":"AIS","cases":19032,"controls":237242,
  "support_url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC13621555/"}
]
def neff(cases,controls,n):
    if cases<=0 or controls<=0 or cases+controls!=n:
        raise ValueError("Case-control count inconsistent with official metadata")
    return 4*cases*controls/n

def build(meta,out):
    table=[]
    for info in COUNTS:
        yamlfile=meta/f'{info["accession"]}_buildGRCh37.tsv.gz-meta.yaml'
        study=yaml.safe_load(yamlfile.read_text())
        entries=study.get("samples") or []
        if len(entries)!=1:raise ValueError("Unexpected sample metadata rows: "+str(yamlfile))
        if entries[0].get("sample_ancestry_category")!=["East Asian"]:
            raise ValueError("Ancestry mismatch: "+str(yamlfile))
        total=int(entries[0]["sample_size"])
        value=neff(info["cases"],info["controls"],total)
        table.append(dict(
            accession=info["accession"],trait=info["trait"],ancestry="EAS",
            cases=info["cases"],controls=info["controls"],
            official_meta_total_n=total,
            approximate_case_control_n_eff=round(value,5),
            working_n_rounded=round(value),
            source_url=info["support_url"],
            method="4*Ncases*Ncontrols/(Ncases+Ncontrols)",
            validation="SAMPLE_TOTAL_MATCHED_OFFICIAL_YAML",
            limitation="STUDY_WIDE_APPROX_NOT_PER_VARIANT_META_EFFECTIVE_N"))
    out.mkdir(parents=True,exist_ok=True)
    with (out/"G0022_GIGASTROKE_EAS_EFFECTIVE_N_APPROX.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=list(table[0]))
        w.writeheader();w.writerows(table)
    result={"status":"SOURCE_COUNT_VALIDATED_APPROXIMATE_NEFF_ONLY",
            "study_records":len(table),"method":"4NcaseNcontrol/(Ncase+Ncontrol)",
            "caution":"Working sensitivity sample size, not true SNP-specific N for GWAS meta analysis.",
            "cases":{r["trait"]:r["cases"] for r in table},
            "controls":{r["trait"]:r["controls"] for r in table},
            "working_n":{r["trait"]:r["working_n_rounded"] for r in table}}
    (out/"G0022_EAS_EFFECTIVE_N_APPROX_SUMMARY.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    return result
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--meta",type=Path,default=META)
    p.add_argument("--out",type=Path,default=OUT)
    args=p.parse_args()
    build(args.meta,args.out)
